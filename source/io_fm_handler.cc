#include "io_fm_handler.h"
#include "layout_info.h"
#include "utils.h"
#include "laph_stdio.h"
#include <unistd.h>


using namespace std;

#define USE_COLLECTIVE_IO


namespace LaphEnv {
   

const int IOFMHandler::ID_string_length = 32;
const IOFMHandler::pos_type IOFMHandler::data_start_pos = 1+IOFMHandler::ID_string_length;
const int IOFMHandler::Nd=4;

ByteHandler IOFMHandler::m_bytehandler;


   // A given site on the 4-d lattice is specified by integers (x,y,z,t)
   // where each of the four integers varies from 
   //
   //     x = 0 .. Nx-1   where   Nx = lattSize()[0],
   //     y = 0 .. Ny-1   where   Ny = lattSize()[1],
   //     z = 0 .. Nz-1   where   Nz = lattSize()[2],
   //     t = 0 .. Nt-1   where   Nt = lattSize()[3]
   //
   // "lexicographic" order is defined such that "x" varies most quickly
   // and "t" varies most slowly.  The lattice is written to disk in
   // lexicographic order.  The lexicographic index for a site is
   // given by
   //           lex_index(x,y,z,t)  =  x + Nx*(y + Ny*(z + Nz*t) )
   //
   // On a parallel machine, the lattice is partitioned up into identical 
   // sublattices.  The "npartitions" parameter specifies how the lattice is 
   // split up.
   //              -npartitions Px Py Pz Pt
   //
   // There are Px*Py*Pz*Pt sublattices.  Each sublattice has extents
   //
   //        Nx/Px by Ny/Py by Nz/Pz by Nt/Pt 
   //
   // where each ratio **must** be an integer.  Let "xinc" = Nx/Px.
   // Start at lex index = 0.  Incrementing by "xinc" always moves you to
   // a different sublattice.  Also, the sites corresponding to the
   // next "xinc" lexico indices must all be on the **same** sublattice.
   //
   // In local host (cpu) memory, the sublattice is laid out
   // lexicographically.  QUDA reorders the sites on the device (gpu)
   // into the preferred layout for the device.

// *************************************************************************

            //  serial code

#ifdef ARCH_SERIAL

void IOFMHandler::delete_file()
{
 check_for_failure(remove(m_filename.c_str()),"Failure deleting file");
}

void IOFMHandler::file_open(IOFMHandler::openmode_type access_mode, 
                            const std::string& errmsg)
{
 fh.open(m_filename.c_str(), std::ios::binary | access_mode);
 check_for_failure(!fh, errmsg);
}
   
void IOFMHandler::file_close()
{
 fh.close();
 check_for_failure(!fh, "Failure during close");
}

const int IOFMHandler::IO_ERR_NO_SUCH_FILE= 1;
const int IOFMHandler::IO_ERR_ACCESS=       1;
const int IOFMHandler::IO_ERR_OTHER=        1;
const int IOFMHandler::IO_SUCCESS=          0;

const IOFMHandler::openmode_type IOFMHandler::IO_MODE_RDONLY= std::ios::in;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_RDWR=   std::ios::in | std::ios::out;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_CREATE= std::ios::trunc;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_EXCL=   std::ios::trunc;

const IOFMHandler::whence_type IOFMHandler::IO_SEEK_BEG  = std::ios::beg;
const IOFMHandler::whence_type IOFMHandler::IO_SEEK_CUR  = std::ios::cur;
const IOFMHandler::whence_type IOFMHandler::IO_SEEK_END  = std::ios::end;


void IOFMHandler::check_for_failure(int errcode, const std::string& mesg)
{
 if (errcode==IOFMHandler::IO_SUCCESS) return;
 errorLaph(make_strf("IOFMHandler error with file %s:\n %s\n",m_filename,mesg));
}

void IOFMHandler::set_striping(int striping_factor, int striping_unit)
{
}


void IOFMHandler::file_seek(off_type offset, whence_type whence)
{
 fh.seekg(offset, whence);
 check_for_failure(fh.fail(),"Failure during seek");
}

IOFMHandler::pos_type IOFMHandler::file_tell()
{
 pos_type current=fh.tellg();
 if (current<0) check_for_failure(1,"Failure during tell");
 return current;
}


void IOFMHandler::printFileHints()
{
 printLaph(make_str("Serial mode: no IOFMHandler file hints\n"));
}


void IOFMHandler::writeCommon(const char *data, size_t nbytes)
{
 pos_type currdisp=fh.tellg();
 fh.seekp(currdisp);  // just to be sure
 fh.write(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common write");
 fh.seekg(currdisp+pos_type(nbytes));
}

 
void IOFMHandler::writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                                   bool datashift)
{
 const char* dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
 writeCommon(dataptr,dinfo.getViewTotalBytes());
}


void IOFMHandler::readCommon(char *data, size_t nbytes, bool broadcast)
{
 pos_type currdisp=fh.tellg();
 fh.read(data,nbytes);
 check_for_failure(fh.fail(),"Failure during common read");
 fh.seekg(currdisp+pos_type(nbytes)); 
}

 
void IOFMHandler::readDistributed(char *data, const DistArrayViewInfo& dinfo,
                                  bool datashift)
{
 char* dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
 readCommon(dataptr,dinfo.getViewTotalBytes());
}

   //  check sum in lexico order in memory, not yet checker boarded
   
void IOFMHandler::compute_lattice_checksum(char* data, size_t bytes_per_site,
                                           size_t lexicostart, size_t nsites)
{
 checksum = m_bytehandler.get_checksum(checksum, data, bytes_per_site*nsites);
}

// *************************************************************************

              //  parallel MPI-IO  code
              
#else


void IOFMHandler::delete_file()
{
 int errcode;
 if ((isPrimaryRank())||(!global_mode)){
    errcode=MPI_File_delete((char*)m_filename.c_str(),MPI_INFO_NULL);}
 if (global_mode) comm_broadcast(&errcode, sizeof(int), 0);
 check_for_failure(errcode,"Failure while deleting file");
}

void IOFMHandler::file_open(IOFMHandler::openmode_type access_mode, 
                            const std::string& errmsg)
{
 if (global_mode)
    check_for_failure(MPI_File_open(MPI_COMM_WORLD, (char*)m_filename.c_str(), 
                      access_mode, finfo, &fh), errmsg);
 else
    check_for_failure(MPI_File_open(MPI_COMM_SELF, (char*)m_filename.c_str(), 
                      access_mode, finfo, &fh), errmsg);
}

void IOFMHandler::file_close()
{
 check_for_failure(MPI_File_close(&fh),"Failure during close");
}

const int IOFMHandler::IO_ERR_NO_SUCH_FILE= MPI_ERR_NO_SUCH_FILE;
const int IOFMHandler::IO_ERR_ACCESS=       MPI_ERR_ACCESS;
const int IOFMHandler::IO_ERR_OTHER=        MPI_ERR_OTHER;
const int IOFMHandler::IO_SUCCESS=          MPI_SUCCESS;

const IOFMHandler::openmode_type IOFMHandler::IO_MODE_RDONLY= MPI_MODE_RDONLY;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_RDWR=   MPI_MODE_RDWR;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_CREATE= MPI_MODE_CREATE;
const IOFMHandler::openmode_type IOFMHandler::IO_MODE_EXCL=   MPI_MODE_EXCL;

const IOFMHandler::whence_type IOFMHandler::IO_SEEK_BEG  = MPI_SEEK_SET;
const IOFMHandler::whence_type IOFMHandler::IO_SEEK_CUR  = MPI_SEEK_CUR;
const IOFMHandler::whence_type IOFMHandler::IO_SEEK_END  = MPI_SEEK_END;


void IOFMHandler::check_for_failure(int errcode, const std::string& mesg)
{
 if (errcode==MPI_SUCCESS) return;
 char errstring[MPI_MAX_ERROR_STRING];
 int len;
 MPI_Error_string(errcode, errstring, &len);
#ifdef USE_COLLECTIVE_IO
 if ((isPrimaryRank())||(!global_mode))
#endif
    {printf("IOFMHandler error with file %s:\n  %s\n%s\n",m_filename.c_str(),
            mesg.c_str(),errstring);}
 errorLaph(make_strf("%d\n",errcode));
}

void IOFMHandler::set_striping(int striping_factor, int striping_unit)
{
 if (striping_factor>1){
    MPI_Info_set(finfo,(char*)"striping_factor",
                (char*) int_to_string(striping_factor).c_str());}
 if (striping_unit>0){
    MPI_Info_set(finfo,(char*)"striping_unit",
                (char*) int_to_string(striping_unit).c_str());}
}

   //  Seek in MPI-IO is relative to the current file view!!
   //  offset 0 refers to the start of the view and is measured
   //  in terms of etypes.  To maintain a similarity with seek
   //  in fstream, it is necessary to return to the default file
   //  view after an file view set in a read/write.

void IOFMHandler::file_seek(off_type offset, whence_type whence)
{
 if (whence==IOFMHandler::IO_SEEK_END){   // MPI_SEEK_END does not seem to work properly
    MPI_Offset fsize;
    check_for_failure(MPI_File_get_size(fh,&fsize),"Failure during seek");
    check_for_failure(MPI_File_seek(fh,fsize+offset,IOFMHandler::IO_SEEK_BEG),
                     "Failure during seek");}
 else{
    check_for_failure(MPI_File_seek(fh,offset,whence),
                      "Failure during seek");}
}

IOFMHandler::pos_type IOFMHandler::file_tell()
{
 pos_type current;
 check_for_failure(MPI_File_get_position(fh,&current),"Failure during tell");
 return current;
}

         // print out MPI-IO file hints
 
void IOFMHandler::printFileHints()
{
 int i,nkeys,flag;
 MPI_Info info_used;
 char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];

 if ((isPrimaryRank())||(!global_mode)){
    std::cout << "File = "<<m_filename<<std::endl<<"File hints:"<<std::endl;
    MPI_File_get_info(fh,&info_used);
    MPI_Info_get_nkeys(info_used, &nkeys);
    for (i=0; i<nkeys; i++){
       MPI_Info_get_nthkey(info_used,i,key);
       MPI_Info_get(info_used,key,MPI_MAX_INFO_VAL,value,&flag);
       std::cout << "    Key = "<<key<<", value = "<<value<<std::endl;}
    MPI_Info_free(&info_used);
    }
}


      // Write at current location in file. Write is done on primary 
      // rank only, but all file pointers updated.  This routine assumes
      // the default file view (entire file, etype = MPI_BYTE).
      // Should be called by all ranks.

void IOFMHandler::writeCommon(const char *data, size_t nbytes)
{
 pos_type currdisp=file_tell();
 MPI_Status status;
 int errcode,count;
 if ((isPrimaryRank())||(!global_mode)){
    errcode=MPI_File_write_at(fh,currdisp,(char*)data,nbytes,MPI_BYTE,&status);
    if (errcode==MPI_SUCCESS){
       errcode=MPI_Get_count(&status,MPI_BYTE,&count);
       if (uint(count)!=nbytes) errcode=MPI_ERR_COUNT;}}
 if (global_mode) comm_broadcast(&errcode, sizeof(int), 0);
 check_for_failure(errcode,"Failure during common write");
 check_for_failure(MPI_File_seek(fh,currdisp+nbytes,MPI_SEEK_SET),
                   "Failure during common write");
}


void IOFMHandler::writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                                   bool datashift)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (!global_mode)
    check_for_failure(IO_ERR_ACCESS,"Distributed write failure: not global mode");
 if (read_mode){
    read_mode=false; checksum=0;}
 pos_type currdisp=file_tell();
              // set the file view
 check_for_failure(MPI_File_set_view(fh,currdisp,
                   dinfo.getFileElemType(),dinfo.getFileViewType(),
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 MPI_Status status;
 int errcode,count;
 char *dataptr=const_cast<char *>(data);
 if (datashift) dataptr+=dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
       // do the collective write!
#ifdef USE_COLLECTIVE_IO
 errcode=MPI_File_write_all(fh,dataptr,dinfo.getViewElementsThisRank(),
                            dinfo.getFileElemType(),&status);
#else
 bool sync=true;
 errcode=MPI_File_write(fh,dataptr,dinfo.getViewElementsThisRank(),
                        dinfo.getFileElemType(),&status);
 comm_broadcast(&sync, sizeof(bool), 0);
#endif
 if (errcode==MPI_SUCCESS){
    errcode=MPI_Get_count(&status,dinfo.getFileElemType(),&count);
    if (count!=dinfo.getViewElementsThisRank()) errcode=MPI_ERR_COUNT;}
 check_for_failure(errcode,"Failure during distributed write");
   // reset view to entire file and etype = bytes so seeks will work like fstream
 check_for_failure(MPI_File_set_view(fh,0,MPI_BYTE,MPI_BYTE,
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 check_for_failure(MPI_File_seek(fh,currdisp+dinfo.getViewTotalBytes(),
                   MPI_SEEK_SET),"Failure during distributed write");
} 


         // Read from current location in file. Memory for "data" 
         // must already be allocated.  File pointer updated.
         // Read is done on primary rank. If "broadcast" is true,
         // data is broadcast to all other ranks. This routine assumes
         // the default file view (entire file, etype = MPI_BYTE).

void IOFMHandler::readCommon(char *data, size_t nbytes, bool broadcast)
{
 pos_type currdisp=file_tell();
 MPI_Status status;
 int errcode,count;
 if ((isPrimaryRank())||(!global_mode)){
    errcode=MPI_File_read_at(fh,currdisp,data,nbytes,MPI_BYTE,&status);
    if (errcode==MPI_SUCCESS){
       errcode=MPI_Get_count(&status,MPI_BYTE,&count);
       if (uint(count)!=nbytes) errcode=MPI_ERR_COUNT;}}
 if (global_mode) comm_broadcast(&errcode, sizeof(int), 0);
 check_for_failure(errcode,"Failure during primary read");
 if ((global_mode)&&(broadcast)) comm_broadcast(data,nbytes,0);
 check_for_failure(MPI_File_seek(fh,currdisp+nbytes,MPI_SEEK_SET),
                   "Failure during primary read");
}


void IOFMHandler::readDistributed(char *data, const DistArrayViewInfo& dinfo,
                                  bool datashift)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!global_mode)
    check_for_failure(IO_ERR_ACCESS,"Distributed read failure: not global mode");
 if (!read_mode){
    read_mode=true; checksum=0;}
 pos_type currdisp=file_tell();
 check_for_failure(MPI_File_set_view(fh,currdisp,
                   dinfo.getFileElemType(),dinfo.getFileViewType(),
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 MPI_Status status;
 int errcode,count;
 char *dataptr=data;
 if (datashift) dataptr+=dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
      // do the collective read!
#ifdef USE_COLLECTIVE_IO
 errcode=MPI_File_read_all(fh,dataptr,dinfo.getViewElementsThisRank(),
                           dinfo.getFileElemType(),&status);
#else
 bool sync=true;
 errcode=MPI_File_read(fh,dataptr,dinfo.getViewElementsThisRank(),
                       dinfo.getFileElemType(),&status);
 comm_broadcast(&sync, sizeof(bool), 0);
#endif
 if (errcode==MPI_SUCCESS){
    errcode=MPI_Get_count(&status,dinfo.getFileElemType(),&count);
    if (count!=dinfo.getViewElementsThisRank()) errcode=MPI_ERR_COUNT;}
 check_for_failure(errcode,"Failure during distributed read");
   // reset view to entire file and etype = bytes so seeks will work like fstream
 check_for_failure(MPI_File_set_view(fh,0,MPI_BYTE,MPI_BYTE,
                   (char*)"native",finfo),"Failure during MPI_File_set_view");
 check_for_failure(MPI_File_seek(fh,currdisp+dinfo.getViewTotalBytes(),
                   MPI_SEEK_SET),"Failure during distributed read");
} 

    //  check sum after distribution to ranks (assume lexicographic ordering)
    
void IOFMHandler::compute_lattice_checksum(char* data, size_t bytes_per_site,
                                           size_t lexicostart, size_t nsites)
{
 uint lastrank=0;
 uint thisrank=LayoutInfo::getMyRank();
 uint rank=0;
 uint count=0;
 const uint xinc = LayoutInfo::getRankLattExtents()[0];
 size_t chunksize=bytes_per_site*xinc;
 size_t lexicostop=lexicostart+nsites;
      // must loop in order of sites
 for (uint site=lexicostart;site<lexicostop;site += xinc){
    rank=LayoutInfo::getRankFromLexicoLinearIndex(site);
    if (lastrank!=rank){
            // send checksum from last rank to new rank
       commSendRecv(&checksum,sizeof(ByteHandler::n_uint32_t),lastrank,rank);
       lastrank=rank;}
    if (thisrank==rank){
       checksum = m_bytehandler.get_checksum(checksum, data+count*bytes_per_site, chunksize);
       count+=xinc;}
    comm_barrier();
    }
 comm_broadcast(&checksum, sizeof(ByteHandler::n_uint32_t), rank);
}


#endif


// *************************************************************************





IOFMHandler::IOFMHandler(bool global) : read_only(true), openflag(false), read_mode(true),
                                    global_mode(global), endian_format('U'), 
                                    endian_convert(false), checksum_on(false), 
                                    is_new_file(false), checksum(0), lattinfo(1),
                                    lattcmplxbytes(LattField::get_cpu_prec_bytes())
{
#ifdef ARCH_PARALLEL
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#endif
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
}


IOFMHandler::IOFMHandler(const std::string& filename, OpenMode mode,
                         const std::string& filetype_id, char endianness,
                         int striping_factor, int striping_unit,
                         bool turn_on_checksum, bool global) : global_mode(global), lattinfo(1),
                         lattcmplxbytes(LattField::get_cpu_prec_bytes())
{
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
#ifdef ARCH_PARALLEL
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}

IOFMHandler::IOFMHandler(const std::string& filename, std::ios_base::openmode mode,
                         const std::string& filetype_id, char endianness,
                         int striping_factor, int striping_unit,
                         bool turn_on_checksum, bool global) : global_mode(global), lattinfo(1),
                         lattcmplxbytes(LattField::get_cpu_prec_bytes())
{
 static_assert(sizeof(int)==4,"sizeof(int) must be 4");
#ifdef ARCH_PARALLEL
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info_create");
#endif
 openflag=false;
 open(filename,mode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOFMHandler::open(const std::string& filename, IOFMHandler::OpenMode mode,
                       const std::string& filetype_id, char endianness,
                       int striping_factor, int striping_unit,
                       bool turn_on_checksum)
{
 close();
 read_only=(mode==ReadOnly)?true:false;

 m_filename=tidyString(filename);
 if (m_filename.empty())
     check_for_failure(IO_ERR_NO_SUCH_FILE,"Empty file name");

 bool exists=file_exists();
 if ((mode==ReadOnly)&&(!exists))
    check_for_failure(IO_ERR_NO_SUCH_FILE,
           "Failure during ReadOnly open: file does not exist");
 else if ((mode==ReadWriteFailIfExists)&&(exists))
    check_for_failure(IO_ERR_ACCESS,
           "Failure during ReadWrite open: file exists and FailIfExists mode");
 else if ((mode==ReadWriteEraseIfExists)&&(exists)){
    delete_file();
    exists=false;}
 if (exists){
    IOFMHandler::openmode_type access=(mode==ReadOnly) ? 
                IO_MODE_RDONLY : IO_MODE_RDWR;
    is_new_file=false;
    open_existing_file(filetype_id,access);}
 else{
    is_new_file=true;
    open_new_file(filetype_id,endianness,striping_factor,striping_unit);}

 checksum_on=turn_on_checksum;
 checksum=0;
 read_mode=true;
}


void IOFMHandler::open(const std::string& filename, std::ios_base::openmode mode,
                       const std::string& filetype_id, char endianness,
                       int striping_factor, int striping_unit,
                       bool turn_on_checksum)
{
 IOFMHandler::OpenMode iomode;
 if (!(mode & std::ios_base::out)) iomode=ReadOnly;
 else{
    if (mode & std::ios_base::trunc) iomode=ReadWriteEraseIfExists;
    else iomode=ReadWriteFailIfExists;}
 open(filename,iomode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOFMHandler::openReadOnly(const std::string& filename, 
                               const std::string& filetype_id,
                               bool turn_on_checksum)
{
 open(filename,ReadOnly,filetype_id,'N',1,0,turn_on_checksum);
}


void IOFMHandler::openNew(const std::string& filename, bool fail_if_exists,
                          const std::string& filetype_id, char endianness,
                          int striping_factor, int striping_unit,
                          bool turn_on_checksum)
{
 OpenMode iomode=(fail_if_exists) ? ReadWriteFailIfExists : ReadWriteEraseIfExists;
 open(filename,iomode,filetype_id,endianness,striping_factor,striping_unit,
      turn_on_checksum);
}


void IOFMHandler::openUpdate(const std::string& filename,
                             const std::string& filetype_id, char endianness,
                             int striping_factor, int striping_unit,
                             bool turn_on_checksum)
{
 open(filename,ReadWriteUpdateIfExists,filetype_id,endianness,
      striping_factor,striping_unit,turn_on_checksum);
}



void IOFMHandler::open_existing_file(const std::string& filetype_id,
                                     IOFMHandler::openmode_type access_mode)
{
 file_open(access_mode,"Could not open existing file: "+m_filename);
     // get endian info, and check file id
 openflag=true;
 string ID_string;
 readIDstring(endian_format,ID_string);  // on primary rank only
 bool flag=false;
 if ((isPrimaryRank())||(!global_mode)){
    flag=true;
    if ((endian_format!='B')&&(endian_format!='L')){
       flag=false; std::cerr << "Invalid endian format"<<std::endl;}
    if (tidyString(filetype_id)!=ID_string){
       flag=false; 
       std::cerr << "File = "<<m_filename<<std::endl;
       std::cerr << "File ID mismatch:"<<std::endl;
       std::cerr << "File contains ID: <"<<ID_string<<">"<<std::endl;
       std::cerr << "ID requested was: <"<<tidyString(filetype_id)<<">"<<std::endl;}}
#ifdef ARCH_PARALLEL
 if (global_mode) comm_broadcast(&flag, sizeof(bool), 0);
#endif
 if (!flag)
    check_for_failure(IO_ERR_OTHER,"Error during open");
#ifdef ARCH_PARALLEL
 if (global_mode) comm_broadcast(&endian_format, sizeof(char), 0);
#endif
 if (m_bytehandler.big_endian())
    endian_convert=(endian_format=='L')?true:false;
 else 
    endian_convert=(endian_format=='B')?true:false;
}


void IOFMHandler::open_new_file(const std::string& filetype_id,
                                char endianness,
                                int striping_factor, int striping_unit)
{
 if (endianness=='N'){  // native 
    endian_format=m_bytehandler.big_endian() ? 'B':'L';
    endian_convert=false;}
 else if ((endianness!='B')&&(endianness!='L'))
    check_for_failure(IO_ERR_OTHER,"Invalid endian format");
 else{
    endian_format=endianness;
    if (m_bytehandler.big_endian()) 
       endian_convert=(endian_format=='L')?true:false;
    else 
       endian_convert=(endian_format=='B')?true:false;}

 set_striping(striping_factor,striping_unit);
 IOFMHandler::openmode_type amode = IO_MODE_RDWR | IO_MODE_CREATE;
 
 file_open(amode,"Could not open file "+m_filename);

 openflag=true;
 writeIDstring(filetype_id);
}


       // Destructor

IOFMHandler::~IOFMHandler() 
{
 if (openflag) clear(); 
#ifdef ARCH_PARALLEL
 MPI_Info_free(&finfo);
#endif
}


     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
bool IOFMHandler::peekString(std::string& stringvalue, size_t byte_offset,
                             const std::string& filename, const std::string& filetype_id)
{
 stringvalue.clear();
 string fname=tidyString(filename);
 if (fname.empty()) return false;
 bool flag=false;
 if ((isPrimaryRank())||(!global_mode)){ 
    flag=peeker(stringvalue,byte_offset,fname,filetype_id);}
#ifdef ARCH_PARALLEL
 if (global_mode){
    comm_broadcast(&flag, sizeof(bool), 0);
    broadcastString(stringvalue);}
#endif
 return flag;
}



bool IOFMHandler::peeker(std::string& stringvalue, size_t byte_offset,
                         const std::string& fname, const std::string& filetype_id)
{
 ifstream in(fname.c_str(), std::ios::binary | std::ios::in);
 if (!in) return false;
 string ID_string(ID_string_length+1,' ');
 in.read((char*)&ID_string[0],ID_string_length+1);
 if (in.fail()) return false;
 char endian=ID_string[0];
 if ((endian!='B')&&(endian!='L')) return false;
 ID_string.erase(0,1);
 ID_string=tidyString(ID_string);
 if (tidyString(filetype_id)!=ID_string) return false;
 in.seekg(pos_type(byte_offset),std::ios::cur);
 if (in.fail()) return false;
 unsigned int n;
 in.read((char*)&n,sizeof(int));
 if (in.fail()) return false;
 if (m_bytehandler.big_endian())
    endian_convert=(endian=='L')?true:false;
 else 
    endian_convert=(endian=='B')?true:false;
 if (endian_convert) m_bytehandler.byte_swap(&n,sizeof(int),1);
 if (n>16777216) return false;  // too large for string...must be corrupt data
 stringvalue.resize(n);
 in.read((char*)&stringvalue[0],sizeof(char)*n);
 if (in.fail()){ stringvalue.clear(); return false;}
 return true;       // in destructor will close the file
}



      // Clear, then reset finfo

void IOFMHandler::close()
{ 
 if (!openflag) return;
 clear();
#ifdef ARCH_PARALLEL
 MPI_Info_free(&finfo);
 check_for_failure(MPI_Info_create(&finfo),"Failure during MPI_Info reset");
#endif
}

     // Close the file, reset all data members except finfo

void IOFMHandler::clear()
{ 
 file_close();
 m_filename.clear();
 openflag=false;
 endian_convert=false;
 endian_format='U';  // undefined
 checksum=0;
 read_only=true;
 read_mode=true;
 checksum_on=false;
 is_new_file=false;
}

         //  Set the file pointer relative to start of file,
         //  end of file (positive is backward), or current location

void IOFMHandler::seekFromStart(off_type offset)
{
 file_seek(offset+data_start_pos,IOFMHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOFMHandler::seekFromCurr(off_type offset)
{ 
 file_seek(offset,IOFMHandler::IO_SEEK_CUR);
 checksum=0;
}
 
void IOFMHandler::seekFromEnd(off_type offset)
{ 
 file_seek(offset,IOFMHandler::IO_SEEK_END);
 checksum=0;
}

void IOFMHandler::seek(pos_type offset)
{ 
 file_seek(offset+data_start_pos,IOFMHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOFMHandler::seekBegin(off_type offset)
{ 
 file_seek(offset+data_start_pos,IOFMHandler::IO_SEEK_BEG);
 checksum=0;
}

void IOFMHandler::seekRelative(off_type offset)
{ 
 file_seek(offset,IOFMHandler::IO_SEEK_CUR);
 checksum=0;
}

void IOFMHandler::seekEnd(off_type offset)
{ 
 file_seek(offset,IOFMHandler::IO_SEEK_END);
 checksum=0;
}

void IOFMHandler::rewind()
{ 
 file_seek(data_start_pos,IOFMHandler::IO_SEEK_BEG);
 checksum=0;
}


      //  Get the file pointer location in bytes from start of data

IOFMHandler::pos_type IOFMHandler::currentPosition()
{
 return file_tell()-data_start_pos;
}

IOFMHandler::pos_type IOFMHandler::tell()
{
 return file_tell()-data_start_pos;
}


void IOFMHandler::turnOnChecksum()
{
 if (checksum_on) return;
 checksum_on=true;
 checksum=0;
}

void IOFMHandler::turnOffChecksum()
{
 checksum_on=false;
}

void IOFMHandler::resetChecksum()
{
 checksum=0;
}

ByteHandler::n_uint32_t IOFMHandler::getChecksum()
{
 if (checksum_on) return checksum;
 check_for_failure(true,"Invalid call to getChecksum since checksums not turned on");
 return checksum;  // to avoid compiler warnings
}



        // Print out file ID on primary rank. Does not 
        // change file pointers.

void IOFMHandler::printFileID()
{
 char endian;
 string ID_string;
 pos_type curr=file_tell();
 readIDstring(endian,ID_string);
 file_seek(curr,IOFMHandler::IO_SEEK_BEG);
 if ((isPrimaryRank())||(!global_mode)){
    std::cout << "File = "<<m_filename<<std::endl<<"ID string = <"
              <<tidyString(ID_string)<<">"<<std::endl;}
}

   // Write ID string at start of file. Length of string
   // is always "ID_string_length".  First character will
   // be 'B' or 'L' to indicate endian format to use.

void IOFMHandler::writeIDstring(const std::string& ID_string)
{
 check_for_failure(int(ID_string.length())>ID_string_length,
                  "IOFMHandler file ID string too long: cannot exceed "
                  +int_to_string(ID_string_length)+" characters");
 string buf(1,endian_format);
 buf+=ID_string;
 int nblanks=ID_string_length-ID_string.length();
 if (nblanks>0)
    buf+=string(nblanks,' ');
 file_seek(0,IOFMHandler::IO_SEEK_BEG);
 writeCommon(&buf[0],ID_string_length+1);
}


   // Reads endian format and ID string from current file.
   // Read is done only on primary rank and is NOT broadcast.
   // Results are returned in "endianness" and "ID_string"
   // on primary rank.  All file pointers updated.

void IOFMHandler::readIDstring(char &endianness, std::string& ID_string)
{
 if ((isPrimaryRank())||(!global_mode)){
    ID_string.resize(ID_string_length+1);}
 file_seek(0,IOFMHandler::IO_SEEK_BEG);
 readCommon(&ID_string[0],ID_string_length+1,false);
 if ((isPrimaryRank())||(!global_mode)){
    endianness=ID_string[0];
    ID_string.erase(0,1);
    ID_string=tidyString(ID_string);}
}


      // Checks if a file exists.
      
bool IOFMHandler::file_exists()
{
 bool result=false;
 if ((isPrimaryRank())||(!global_mode)){
    result = (access(m_filename.c_str(),F_OK) == 0) ? true : false;}
#ifdef ARCH_PARALLEL
 if (global_mode){
    comm_broadcast(&result, sizeof(bool), 0);}
#endif
 return result;
}


      // Converts an integer to a string
      
string IOFMHandler::int_to_string(int intval)
{ 
 std::ostringstream oss;
 oss << intval;
 return oss.str();
}

      // Removes leading and trailing blanks in a string

string IOFMHandler::tidyString(const string& str)   
{
 string tmp;
 for (unsigned int i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 size_t start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 size_t len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}


// *******************************************************************

      //   Main input/output routines that do the byte-swapping (if needed),
      //   update the check sum, and do the read/write.

void IOFMHandler::write_common(const char* output, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
 if (endian_convert){
    m_bytehandler.byte_swap(const_cast<char *>(output), elementbytes, nelements);
    if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);
    m_bytehandler.byte_swap(const_cast<char *>(output), elementbytes, nelements);}
 else{
    if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, output, elementbytes*nelements);
    writeCommon(output, elementbytes*nelements);}
}


void IOFMHandler::read_common(char* input, size_t elementbytes, size_t nelements)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 readCommon(input, elementbytes*nelements);
 if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, input, elementbytes*nelements);
 if (endian_convert) m_bytehandler.byte_swap(input, elementbytes, nelements);
}


      //   Main input/output routines that do the byte-swapping (if needed),
      //   update the check sum, and do the read/write.  "bytes_per_site" is the total number
      //   of bytes needed to define the entire quantity on one lattice site (includes
      //   color matrix, spin vectors, etc.).  "bytes_per_word" 
      //   is the number of bytes in the underlying numbers or words (int, double, complex, 
      //   etc) of each component of the quantity at each site.  This information is needed
      //   for byte-swapping.

void IOFMHandler::write_lattice(const char* output, size_t bytes_per_word, 
                                const DistArrayViewInfo& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Write failure--no open file");
 if (read_only)
    check_for_failure(IO_ERR_ACCESS,"Write failure--read only file");
 if (read_mode){
    read_mode=false; checksum=0;}
      // convert from even-odd checkerboard to lexico
 vector<char> lexbuf(LayoutInfo::getRankLatticeNumSites()*dinfo.getBytesPerSite());
 LayoutInfo::evenodd_to_lexico(lexbuf.data(),output,dinfo.getBytesPerSite());
      // set up info for writing
 IOH_int bytes_this_rank=dinfo.getViewBytesThisRank();
 char* buf=lexbuf.data();
 bool datashift=true;
 char* b=buf+dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
 if (endian_convert){
    size_t words_this_rank=bytes_this_rank/bytes_per_word;
    m_bytehandler.byte_swap(b, bytes_per_word, words_this_rank);
    if (checksum_on) 
       compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                                dinfo.getViewLexicoSpan());
    writeDistributed(buf, dinfo, datashift);
    if (datashift) 
       m_bytehandler.byte_swap(b, bytes_per_word, words_this_rank);
    }
 else{ 
    if (checksum_on){
       compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                                dinfo.getViewLexicoSpan());}
    writeDistributed(buf,dinfo,datashift);}
}



void IOFMHandler::read_lattice(char* input, size_t bytes_per_word, 
                               const DistArrayViewInfo& dinfo)
{
 if (!openflag)
    check_for_failure(IO_ERR_OTHER,"Read failure--no open file");
 if (!read_mode){
    read_mode=true; checksum=0;}
 vector<char> lexbuf(LayoutInfo::getRankLatticeNumSites()*dinfo.getBytesPerSite());
 IOH_int bytes_this_rank=dinfo.getViewBytesThisRank();
 char* buf=lexbuf.data();
 bool datashift=true;
 char* b=buf+dinfo.getElementOffsetThisRank()*dinfo.getBytesPerElement();
 readDistributed(buf, dinfo, datashift);
 if (checksum_on)
    compute_lattice_checksum(b,dinfo.getBytesPerSite(),dinfo.getViewLexicoStart(),
                             dinfo.getViewLexicoSpan());
 if (endian_convert){
    size_t words_this_rank=bytes_this_rank/bytes_per_word;
    m_bytehandler.byte_swap(b, bytes_per_word, words_this_rank);}
    // convert from lexico to even-odd checkerboard
 LayoutInfo::lexico_to_evenodd(input,buf,dinfo.getBytesPerSite());
}



void IOFMHandler::write(const LattField& output)
{
 int sizeT=output.bytesPerSite();
 lattinfo.resetBytes(sizeT);
 unsigned int sz[IOFMHandler::Nd+2];
 sz[0]=LayoutInfo::Ndim;
 for (int k=0;k<LayoutInfo::Ndim;k++)
    sz[k+1]=LayoutInfo::getLattExtents()[k];
 sz[LayoutInfo::Ndim+1]=sizeT;
 multi_write(sz,LayoutInfo::Ndim+2);
 write_lattice((const char *)(output.getDataConstPtr()),output.bytesPerWord(), lattinfo);
}

void IOFMHandler::write(const std::vector<LattField>& output)
{
 int n=output.size();
 write(n);
 if (n==0) return;
 for (int k=0;k<n;k++) write(output[k]);
}

void IOFMHandler::read(LattField& input)
{
 unsigned int sz[IOFMHandler::Nd+2];
 multi_read(sz,LayoutInfo::Ndim+2);
 bool errflag=(sz[0]!=(unsigned int)(LayoutInfo::Ndim));
 for (int k=0;k<LayoutInfo::Ndim;k++)
    errflag=errflag || (sz[k+1]!=(unsigned int)(LayoutInfo::getLattExtents()[k]));
 check_for_failure(errflag, "Invalid lattice read");
 int sizeT=sz[IOFMHandler::Nd+1];
 lattinfo.resetBytes(sizeT);
 input.reset_by_bytes_per_site(sizeT);  // read using precision in file
 int cmplxbytes=input.bytesPerSite()/input.elemsPerSite();
 read_lattice((char *)(input.getDataPtr()),cmplxbytes, lattinfo);
 if (cmplxbytes!=lattcmplxbytes){
    lattcmplxbytes=cmplxbytes;}   // done for IOMap to help with numbytes
 input.to_quda_precision();   // convert to quda precision if needed
}


void IOFMHandler::read(std::vector<LattField>& input)
{
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>1024), "std::vector too large...bad file location?");
 input.resize(n);
 for (int k=0;k<n;k++) read(input[k]);
}

size_t IOFMHandler::numbytes(const LattField& data)
{
 return (LayoutInfo::Ndim+2)*sizeof(int)+lattcmplxbytes*data.elemsPerSite()
         *LayoutInfo::getLatticeNumSites();
}

size_t IOFMHandler::numbytes(const std::vector<LattField>& data)
{
 if (data.size()>0){
    size_t bytes=sizeof(int);
    for (uint k=0;k<data.size();++k){
       bytes+=numbytes(data[k]);}
    return bytes;}
 return sizeof(int);
}

// *******************************************************************


void IOFMHandler::write(const std::string& output)
{
 int n=output.length();
 write_common((const char*)&n, sizeof(int), 1);
 write_common(output.data(), sizeof(char), n);
}


void IOFMHandler::read(std::string& input, bool broadcast)
{
 int n;
 read_common((char*)&n, sizeof(int), 1);
    // if reading in wrong location, could get nonsense here,
    // so limit to a 16MB string
 check_for_failure(n>16777216, "string for read too large...bad file location?");
 char* str = new(nothrow) char[n];
 check_for_failure((str==0),"IOFMHandler::read---unable to allocate memory");
 readCommon(str, sizeof(char)*n, broadcast);
 if ((isPrimaryRank())||(!global_mode)||(global_mode && broadcast)){
    if (checksum_on) checksum = m_bytehandler.get_checksum(checksum, str, sizeof(char)*n);
    input.assign(str, n);}
#ifdef ARCH_PARALLEL
 if (global_mode && (!broadcast) && checksum_on) comm_broadcast(&checksum, sizeof(checksum), 0);
#endif
 delete[] str;
}


    // the number of bytes that these quantities occupy in an
    // IOFMHandler file

size_t numbytes(IOFMHandler& ioh, const std::string& output)
{
 return sizeof(int)+output.length();
}


// *************************************************


#ifdef ARCH_PARALLEL

 // **************************************************************************

  //  Objects of this class contain info about multi-dimensional arrays 
  //  distributed across the MPI ranks.  Blocking assumed, column major site ordering.   
  //  Assumed distributed in a regular way...same size on each rank.  
  //  No distribution of the MPI ranks is assumed.  This class is mainly
  //  used by BinaryFileHandler for creating file views for MPI-IO. 
      
  //  The main purpose of this class is to create an appropriate "filetype"
  //  to pass to MPI-IO.  The routine  MPI_Type_create_subarray is used.


         // copy constructors
        
DistArrayViewInfo::DistArrayViewInfo(const DistArrayViewInfo& in) 
                     :    ndim(in.ndim), global_sizes(in.global_sizes),
                          local_sizes(in.local_sizes),
                          start_indices(in.start_indices),
                          nelem_this_rank(in.nelem_this_rank),
                          gvvol(in.gvvol), lexico_start(in.lexico_start),
                          el_offset_this_rank(in.el_offset_this_rank),
                          nbytes_per_element(in.nbytes_per_element)
{
 MPI_Type_dup(in.ftype,&ftype);
 MPI_Type_dup(in.etype,&etype);
}
    

DistArrayViewInfo& DistArrayViewInfo::operator=(const DistArrayViewInfo& in)
{
 ndim=in.ndim; 
 global_sizes=in.global_sizes;
 local_sizes=in.local_sizes;
 start_indices=in.start_indices;
 nelem_this_rank=in.nelem_this_rank;
 gvvol=in.gvvol; 
 lexico_start=in.lexico_start;
 el_offset_this_rank=in.el_offset_this_rank;
 nbytes_per_element=in.nbytes_per_element;
 MPI_Type_dup(in.ftype,&ftype);
 MPI_Type_dup(in.etype,&etype);
 return *this;
}

         // destructor

DistArrayViewInfo::~DistArrayViewInfo()
{
 MPI_Type_free(&ftype);
 MPI_Type_free(&etype);
}

   //  Sets up the MPI file view for this MPI rank assuming a full lattice 
   //  quantity.  "gsizes" are the global lattice sizes (the entire lattice).
   //  "numranks" should contains the number of partitions in each of the
   //  space-time dimensions.  "comm_coord" contains the communications
   //  coordinate of this rank. 

void DistArrayViewInfo::setup_full(const vector<int>& gsizes, const vector<int>& numranks,
                                   const vector<int>& comm_coord)
{
 ndim=gsizes.size();
 global_sizes=gsizes;
 local_sizes.resize(ndim);
 start_indices.resize(ndim);
 IOH_int lavol=1;
 gvvol=1;
 for (int k=0;k<ndim;k++){
    local_sizes[k]=global_sizes[k]/numranks[k];
    start_indices[k]=comm_coord[k]*local_sizes[k];
    gvvol*=global_sizes[k];
    lavol*=local_sizes[k];}
 nelem_this_rank=lavol;
 lexico_start=0;
 el_offset_this_rank=0;
}                      


void DistArrayViewInfo::set_bytes(IOH_int bytes_per_element)
{
 error_return(bytes_per_element<=0,"invalid input set_bytes");
 nbytes_per_element=bytes_per_element;
}



bool DistArrayViewInfo::checker(const vector<int>& gsizes, 
                                const vector<int>& numranks)
{
 if ((gsizes.size()==0)||(gsizes.size()!=numranks.size()))
    return false;
 for (int k=0;k<int(gsizes.size());k++)
    if ((gsizes[k]<=0)||(numranks[k]<=0))
       return false;
 for (int k=0;k<int(gsizes.size());k++){
    int ll=gsizes[k]/numranks[k];
    if (numranks[k]*ll != gsizes[k]) return false;}
 int check=1;
 for (int k=0;k<int(gsizes.size());k++) check*=numranks[k];
 if (check!=int(comm_size())) return false;
 return true;
}


void DistArrayViewInfo::error_return(bool cond, const string& msg)
{
 if (cond){
    errorLaph(make_strf("Error in DistArrayViewInfo: %s\n",msg));}
}

         // constructor for full lattice quantity
         
DistArrayViewInfo::DistArrayViewInfo(IOH_int bytes_per_site)
{
 setup_full(LayoutInfo::getLattExtents(),LayoutInfo::getCommNumPartitions(),LayoutInfo::getMyCommCoords());
 set_bytes(bytes_per_site);
    // make and store the MPI filetype
 create_filetype();
} 

         // constructor for distributed array (non-lattice), let MPI
         // distribute over the ranks
         

DistArrayViewInfo::DistArrayViewInfo(const vector<int>& gsizes, 
                                     IOH_int bytes_per_element)
{
 int nranks=comm_size();
 int nd=gsizes.size();
 error_return(nd==0,"bad global sizes");
 vector<int> numranks(nd,0);
 error_return(MPI_Dims_create(nranks,nd,&numranks[0]),
          "Could not distribute array among processors");
 error_return(!checker(gsizes,numranks),"invalid input");
 vector<int> comm_coord(LayoutInfo::getMyCommCoords()); 
 setup_full(gsizes,numranks,comm_coord);
 set_bytes(bytes_per_element);
    // make and store the MPI filetype
 create_filetype();
}


DistArrayViewInfo::DistArrayViewInfo(const vector<int>& gsizes, 
                                     const vector<int>& numranks,
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes,numranks),"invalid input");
 vector<int> comm_coord(LayoutInfo::getMyCommCoords()); 
 setup_full(gsizes,numranks,comm_coord);
 set_bytes(bytes_per_element);
    // make and store the MPI filetype
 create_filetype();
}


void DistArrayViewInfo::create_filetype()
{
 MPI_Type_contiguous(nbytes_per_element,MPI_BYTE,&etype);
 MPI_Type_commit(&etype);
 error_return(MPI_Type_create_subarray(ndim,&global_sizes[0],&local_sizes[0],
                  &start_indices[0],MPI_ORDER_FORTRAN,etype,&ftype),
                  "Could not create filetype");            
 MPI_Type_commit(&ftype);
}


void DistArrayViewInfo::resetBytes(IOH_int bytes_per_element)
{
 if (bytes_per_element==nbytes_per_element) return;
 set_bytes(bytes_per_element);
 MPI_Type_free(&ftype);
 MPI_Type_free(&etype);
 create_filetype();
}


#else
 // ***************************************************************
               // serial code

DistArrayViewInfo::DistArrayViewInfo(const DistArrayViewInfo& in) 
                     :    ndim(in.ndim), local_sizes(in.local_sizes),
                          gvvol(in.gvvol), lexico_start(in.lexico_start),
                          nbytes_per_element(in.nbytes_per_element) {}
    

DistArrayViewInfo& DistArrayViewInfo::operator=(const DistArrayViewInfo& in)
{
 ndim=in.ndim; 
 local_sizes=in.local_sizes;
 gvvol=in.gvvol; 
 lexico_start=in.lexico_start;
 nbytes_per_element=in.nbytes_per_element;
 return *this;
}

         // destructor

DistArrayViewInfo::~DistArrayViewInfo() {}

 

void DistArrayViewInfo::setup_full(const vector<int>& gsizes)
{
 ndim=gsizes.size();
 local_sizes=gsizes;
 gvvol=1;
 for (int k=0;k<ndim;k++)
    gvvol*=gsizes[k];
 lexico_start=0;
}                      


void DistArrayViewInfo::set_bytes(IOH_int bytes_per_element)
{
 error_return(bytes_per_element<=0,"invalid input set_bytes");
 nbytes_per_element=bytes_per_element;
}



bool DistArrayViewInfo::checker(const vector<int>& gsizes)
{
 if (gsizes.size()==0)
    return false;
 for (int k=0;k<int(gsizes.size());k++)
    if (gsizes[k]<=0) return false;
 return true;
}


void DistArrayViewInfo::error_return(bool cond, const string& msg)
{
 if (cond){
    errorLaph(make_strf("Error in DistArrayViewInfo: %s\n",msg));}
}

         // constructor for full lattice quantity
         
DistArrayViewInfo::DistArrayViewInfo(IOH_int bytes_per_site)
{
 setup_full(LayoutInfo::getLattExtents());
 set_bytes(bytes_per_site);
} 


         // constructor for distributed array (non-lattice), let MPI
         // distribute over the ranks
         

DistArrayViewInfo::DistArrayViewInfo(const vector<int>& gsizes, 
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes),"invalid input");
 setup_full(gsizes);
 set_bytes(bytes_per_element);
}

DistArrayViewInfo::DistArrayViewInfo(const vector<int>& gsizes, 
                                     const vector<int>& numranks,
                                     IOH_int bytes_per_element)
{
 error_return(!checker(gsizes),"invalid input");
 setup_full(gsizes);
 set_bytes(bytes_per_element);
}


void DistArrayViewInfo::resetBytes(IOH_int bytes_per_element)
{
 if (bytes_per_element==nbytes_per_element) return;
 set_bytes(bytes_per_element);
}



#endif


// ***************************************************************
}
