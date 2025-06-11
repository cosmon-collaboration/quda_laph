#ifndef IO_FM_HANDLER_H
#define IO_FM_HANDLER_H

#include "array.h"
#include "byte_handler.h"
#include "latt_field.h"
#include <fstream>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

typedef std::complex<double> cmplx;
typedef std::complex<float> fcmplx;


namespace LaphEnv {



 // *********************************************************************************
 // *                                                                               *
 // *       class IOFMHandler:    parallel random access input/output               *
 // *                                                                               *
 // *   Author: Colin Morningstar (Carnegie Mellon University)                      *
 // *                                                                               *
 // *   This class attempts to make serial/parallel input/output in local and       *
 // *   global mode easier to achieve.  This is a low-level class which has         *
 // *   the functionality of fstreams except it uses MPI-IO in a parallel setting.  *
 // *   The FM stands for "fstreams or mpi-io" since fstreams is used when          *
 // *   compiled for serial runs and MPI-IO is used when compiled for parallel      *
 // *   runs.                                                                       *
 // *                                                                               *
 // *   Some features of this class are as follows: (a) in a parallel architecture, *
 // *   it uses MPI-IO so that input/output is parallel and efficient; (b) it       *
 // *   allows the user to choose between big-endian and little-endian format;      *
 // *   (c) it allows the user to turn check sums on or off; (d) MPI-IO striping    *
 // *   can be used in a parallel architecture; (e) objects are constructed and     *
 // *   subsequently operate in either global or local mode.  Currently, only       *
 // *   blocking I/O operations are supported, but it would be simple to implement  *
 // *   some of MPI-IO's non-blocking operations.                                   *
 // *                                                                               *

 // *   NOTE: this is a utility class used by IOFMMap, which is the class meant     *
 // *   for the end user.  This class uses primitive seek functions to read/write,  *
 // *   just like fstreams.  IOFMMap acts like a C++ map but with storage in a      *
 // *   file, and IOFMMap allows access to the data through a key such as in a C++  *
 // *   map.  In this way, IOFMMap is meant to mimic the features of HDF5, except   *
 // *   that the key need not be just a string.  IOFMHandler can maintain a         *
 // *   checksum during I/O operations, but it does not write checksums into the    *
 // *   files, whereas IOFMMap does.                                                *
 // *                                                                               *
 // *   An object of class IOFMHandler must be defined as global or local in the    *
 // *   constructor, and this property cannot be changed.  If global, then it is    *
 // *   assumed that an object is created on every rank.  In global mode, lattice   *
 // *   and other distributed data can be read and written, and common data is      *
 // *   written and read according to how MPI-IO does common read/writes.  Only     *
 // *   one copy is written to file because the data is assumed to be common or     *
 // *   the same on all ranks.  After a common read, the data will be available     *
 // *   on all ranks.  In local mode, the object is assumed to exist only on the    *
 // *   current rank. Different ranks will write to different files.  Reads can     *
 // *   occur from different files. Only non-distributed data can be read and       *
 // *   written.  Read/writes are done through the current rank, and no data is     *
 // *   broadcast to other ranks.                                                   *
 // *                                                                               *
 // *   Files written by this class will start with a single character 'B' or 'L'   *
 // *   to indicate endian-ness.  Then an ID string is output, which always has     *
 // *   the length "ID_string_length"=32 (spaces are padded if needed).  Files read *
 // *   by this class expect this starting behavior.  An error occurs if the        *
 // *   ID string given in the open command on an existing file does not match      *
 // *   that in the file.                                                           *
 // *                                                                               *
 // *   All errors are considered fatal in this class, so the end user need not     *
 // *   perform any error checking.  The ID string is useful for ensuring that      *
 // *   the file contains the kind of information that the user is expecting.       *
 // *                                                                               *
 // *   Files can be opened in one of four modes:                                   *
 // *     (1) ReadOnly  -- fails if the file does not exist or read not allowed     *
 // *     (2) ReadWriteFailIfExists --  fails if the file exists                    *
 // *     (3) ReadWriteEraseIfExists  -- erases the file if it exists               *
 // *     (4) ReadWriteUpdateIfExists  -- updates existing file or creates new      *
 // *   Modes 2,3,4 will create a new file if it does not exist.  When creating     *
 // *   a new file, the endian format choices are 'B', 'L', or 'N' (native).        *
 // *   Random access to the file is employed.                                      *
 // *                                                                               *
 // *   There are separate open routines for the different modes, or you can        *
 // *   use the OpenMode enum in the general open routine.  Alternatively, you      *
 // *   can use iostream openmodes as follows:                                      *
 // *      ReadOnly                  =   ios::in                                    *
 // *      ReadWriteFailIfExists     =   ios::in | ios::out                         *
 // *      ReadWriteEraseIfExists    =   ios::in | ios::out | ios::trunc            *
 // *                                                                               *
 // *   There are seek and tell members for moving around in the file.  It is an    *
 // *   error to seek before the start of the data, but you can seek past the end   *
 // *   of the file.  You can seek relative to the start of the data (33 characters *
 // *   past the start of the file), the current location, or the end of file.      *
 // *   Use negative offsets relative to the end of file to go backwards.           *
 // *                                                                               *
 // *   Input is done with read(..) commands, and output is done with write(..)     *
 // *   commands. For lattice objects, collective I/O is used (in global mode). For *
 // *   non-lattice objects, the read() and write() commands are either done by     *
 // *   the primary rank (global mode) or the current rank (local mode).            *
 // *   The data types currently supported for read/write are                       *
 // *                                                                               *
 // *       - basic data types (int, bool, float, string, etc.)                     *
 // *       - vector and Array of basic data types                                  *
 // *       - LattField objects                                                     *
 // *       - vector of LattField objects                                           *
 // *                                                                               *
 // *   The quantities above are contiguous in memory and regular (each element is  *
 // *   the same size) and these facts have been used to speed up I/O.              *
 // *                                                                               *
 // *   (a) For lattice objects, the read() and write() commands do collective I/O  *
 // *   in MPI-IO. All ranks participate so all ranks must call the read/write      *
 // *   command.  Before the read/write call, ensure that all file pointers are set *
 // *   to the same point in the file (such as by a seek call by all ranks).  After *
 // *   each read/write, all file pointers are advanced by the size of the          *
 // *   lattice object.  Only applies to global mode.                               *
 // *                                                                               *
 // *   (b) For reading/writing non-lattice objects in global mode, all ranks must  *
 // *   call the read/write subroutine, but I/O is only done on the primary rank    *
 // *   (or however MPI-IO does such input/output). For reads, the results are      *
 // *   broadcast to all ranks.  All file pointers are advanced by the size of the  *
 // *   object after the I/O operation.  Before the read/write, it is safest to     *
 // *   have all ranks call a seek operation so that all file pointers are          *
 // *   initially pointing to the same location in the file.                        *
 // *                                                                               *
 // *   (c) Reading/writing non-lattice objects in "local" mode must be done        *
 // *   carefully.  Typically, a seek is called only on one rank, and the local     *
 // *   read/write is only called by the one rank.  The file pointer for that       *
 // *   one rank will be advanced afterwards by the size of the object. For local   *
 // *   reads, the results are NOT broadcast to all other ranks.                    *
 // *   Simultaneous local writes to the same file can cause major problems with    *
 // *   consistency!!  The "local" mode is more meant to be used during             *
 // *   simultaneous writes/reads to **different** files.                           *
 // *   Use local mode cautiously, especially when writing data.                    *
 // *                                                                               *
 // *   Examples of writing and reading:                                            *
 // *                                                                               *
 // *      IOFMHandler io;     <-- default is global mode                           *
 // *      io.open(....);                                                           *
 // *      int k=5;     io.seekFromStart(...); write(io,k);                         *
 // *      float x=5.4; io.seekFromStart(...); write(io,x);                         *
 // *      int j; float y;  io.seekFromStart(...); read(io,j,y); // in sequence     *
 // *      string str("ABC"); io.seekFromStart(...); write(io,str);                 *
 // *      multi1d<float> g(1024); io.seekFromStart(...); write(io,g);              *
 // *      LattField zn(...);  io.seekFromStart(...); write(io,z);                  *
 // *                                                                               *
 // *   When reading and writing distributed arrays, a column-major order is        *
 // *   assumed in the file (MPI_ORDER_FORTRAN), and column-major order of the      *
 // *   sublattice on each rank is assumed.  No ordering of the ranks is assumed.   *
 // *                                                                               *
 // *   If check sums are turned on, objects of this class maintain a check sum.    *
 // *   Doing an explicit seek resets the checksum; changing from a read to a       *
 // *   write or vice versa resets the checksum; or an explicit reset can also      *
 // *   be done.  The checksum is updated as bytes are read or written.             *
 // *   Successive reads or successive writes update the checksum. Checksums are    *
 // *   NOT stored in the file.  IOMap uses IOFMHandler, and IOMap maintains        *
 // *   checksums in the files if desired.                                          *
 // *                                                                               *
 // *   Precision conversion takes place automatically for lattice field reads.     *
 // *   Lattice field writes can only take place using the current quda precision.  *
 // *                                                                               *
 // *   For local read/writes, create an IOFMHandler object on one processor in     *
 // *   local mode:                                                                 *
 // *                                                                               *
 // *      bool global_mode=false;                                                  *
 // *      IOFMHandler localhandler(global_mode);    <-  local write                *
 // *                                                                               *
 // *   Use local I/O **very carefully**.                                           *
 // *                                                                               *
 // *********************************************************************************


 //  IOFMHandler class is defined later in this file, after the 
 //  helper class DistArrayViewInfo.


#if defined(ARCH_SERIAL)
   typedef size_t       IOH_int; 
#elif defined(ARCH_PARALLEL)
   typedef MPI_Offset   IOH_int;
#else
   #error "Either ARCH_PARALLEL or ARCH_SERIAL must be defined"
#endif

 
  /*  DistArrayViewInfo objects contain info about multi-dimensional arrays 
      distributed across the MPI ranks.  Blocking assumed, column major.   
      Assumed distributed in a regular way...same size on each rank.  
      Distribution on ranks is column major.  This class is mainly
      used by IOFMHandler for creating file views for MPI-IO. It supports
      entire global arrays, and single slices of the most major index
      (the right-most index).  The end user does not need to know
      anything about this class.
      
      The main purpose of this class is to create an appropriate "filetype"
      to pass to MPI-IO.  The routine  MPI_Type_create_subarray is used.  
      
      Warning:  This class stores an MPI quantity, so make sure all
      destructors get called before MPI_finalize().   */

#ifdef ARCH_PARALLEL
                     //   parallel version

class DistArrayViewInfo
{

    int ndim;
    std::vector<int> global_sizes, local_sizes, start_indices;
    IOH_int nelem_this_rank, gvvol, lexico_start, el_offset_this_rank;
    IOH_int nbytes_per_element;
    MPI_Datatype ftype;
    MPI_Datatype etype;
    
  public:
  
    DistArrayViewInfo(IOH_int bytes_per_site);   // assumes a full lattice
    
    DistArrayViewInfo(int time_slice, IOH_int bytes_per_site);   // for a lattice time slice

                 // this lets MPI determine how to distribute arrays on ranks
    DistArrayViewInfo(const std::vector<int>& gsizes, IOH_int bytes_per_element);
    
    DistArrayViewInfo(const std::vector<int>& gsizes, const std::vector<int>& numranks,
                      IOH_int bytes_per_element);

    DistArrayViewInfo(const std::vector<int>& gsizes, const std::vector<int>& numranks,
                      int major_index, IOH_int bytes_per_element);

    void resetBytes(IOH_int bytes_per_site);

    DistArrayViewInfo(const DistArrayViewInfo& in);

    DistArrayViewInfo& operator=(const DistArrayViewInfo& in);

    ~DistArrayViewInfo();

   

    int getArrayDimensions() const
     { return ndim; }
     
    const std::vector<int>& getGlobalArraySizes() const
     { return global_sizes; }

    const std::vector<int>& getLocalArraySizes() const
     { return local_sizes; }

    const std::vector<int>& getViewLocalStarts() const
     { return start_indices; }
     
    IOH_int getBytesPerSite() const
     { return nbytes_per_element; }
     
    IOH_int getBytesPerElement() const
     { return nbytes_per_element; }
     
    IOH_int getViewElementsThisRank() const
     { return nelem_this_rank; }

    IOH_int getViewBytesThisRank() const
     { return nbytes_per_element*nelem_this_rank; }

    IOH_int getElementOffsetThisRank() const
     { return el_offset_this_rank; }

    IOH_int getViewTotalBytes() const
     { return gvvol*nbytes_per_element; }
     
    IOH_int getViewLexicoStart() const   // in terms of elements
     { return lexico_start; }            // helps with check sums
    
    IOH_int getViewLexicoSpan() const    // in terms of elements
     { return gvvol; }                   // helps with check sums


    const MPI_Datatype& getFileViewType() const
     { return ftype; }

    const MPI_Datatype& getFileElemType() const
     { return etype; }
    

 
 private:
 
    void setup_full(const std::vector<int>& gsizes, const std::vector<int>& numranks,
                    const std::vector<int>& rank_coord);

    void set_bytes(IOH_int bytes_per_element);
    
    bool checker(const std::vector<int>& gsizes, const std::vector<int>& numranks);
    
    void error_return(bool errcond, const std::string& mesg);


    void create_filetype();
     
};

#else
                     //   serial version
 
class DistArrayViewInfo
{

    int ndim;
    std::vector<int> local_sizes;
    IOH_int gvvol, lexico_start;
    IOH_int nbytes_per_element;
    
  public:
  
    DistArrayViewInfo(IOH_int bytes_per_site);   // assumes a full lattice
    
                 // this lets MPI determine how to distributes array on ranks
    DistArrayViewInfo(const std::vector<int>& gsizes, IOH_int bytes_per_element);
    
             // numranks ignored below...for compatibility with parallel compiles
    DistArrayViewInfo(const std::vector<int>& gsizes, const std::vector<int>& numranks,
                      IOH_int bytes_per_element);

    void resetBytes(IOH_int bytes_per_site);

    DistArrayViewInfo(const DistArrayViewInfo& in);

    DistArrayViewInfo& operator=(const DistArrayViewInfo& in);

    ~DistArrayViewInfo();

   
   
    int getArrayDimensions() const
     { return ndim; }
     
    const std::vector<int>& getGlobalArraySizes() const
     { return local_sizes; }

    const std::vector<int>& getLocalArraySizes() const
     { return local_sizes; }

    const std::vector<int> getViewLocalStarts() const
     { std::vector<int> tmp(ndim,0); return tmp; }

    IOH_int getBytesPerSite() const
     { return nbytes_per_element; }
     
    IOH_int getBytesPerElement() const
     { return nbytes_per_element; }
     
    IOH_int getViewElementsThisRank() const
     { return gvvol; }

    IOH_int getViewBytesThisRank() const
     { return nbytes_per_element*gvvol; }

    IOH_int getElementOffsetThisRank() const
     { return lexico_start; }

    IOH_int getViewTotalBytes() const
     { return gvvol*nbytes_per_element; }
     
    IOH_int getViewLexicoStart() const   // in terms of elements
     { return lexico_start; }        // helps with check sums
    
    IOH_int getViewLexicoSpan() const    // in terms of elements
     { return gvvol; }               // helps with check sums



 
 private:
 
    void setup_full(const std::vector<int>& gsizes);

    void set_bytes(IOH_int bytes_per_element);
    
    bool checker(const std::vector<int>& gsizes);
    
    void error_return(bool errcond, const std::string& mesg);

     
};

#endif


 // ***********************************************************
 // *                                                         *
 // *          Now for the main event:  IOFMHandler           *
 // *                                                         *
 // ***********************************************************


class IOFMHandler
{

#if defined(ARCH_SERIAL)
   std::fstream fh;
#elif defined(ARCH_PARALLEL)
   MPI_File fh;        // the MPI-IO file handler
   MPI_Info finfo;     // MPI-IO file hints
#else
   #error "Either ARCH_PARALLEL or ARCH_SERIAL must be defined"
#endif

   bool read_only;
   bool openflag;
   bool read_mode;
   bool global_mode;
   
   char endian_format;            // 'B' for big-endian, 'L' for little-endian
   bool endian_convert;

   bool checksum_on;
   std::string m_filename;
   bool is_new_file;
 
   ByteHandler::n_uint32_t checksum;
   DistArrayViewInfo lattinfo;
   int lattcmplxbytes;    // needed to handle precision conversion in IOMap

      // disallow copying
   IOFMHandler(const IOFMHandler&);
   IOFMHandler(IOFMHandler&);
   IOFMHandler& operator=(const IOFMHandler&);
   IOFMHandler& operator=(IOFMHandler&);


 public:
 
   enum OpenMode { ReadOnly, ReadWriteFailIfExists, ReadWriteEraseIfExists, 
                   ReadWriteUpdateIfExists };

   explicit IOFMHandler(bool global=true);

   IOFMHandler(const std::string& filename, OpenMode mode=ReadOnly,
               const std::string& filetype_id="", bool global=true, 
               char endianness='N', bool turn_on_checksum=false,
               int striping_factor=1, int striping_unit=0);

   IOFMHandler(const std::string& filename, std::ios_base::openmode mode,
               const std::string& filetype_id="", bool global=true,
               char endianness='N', bool turn_on_checksum=false,
               int striping_factor=1, int striping_unit=0);

   void open(const std::string& filename, OpenMode mode=ReadOnly,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false, int striping_factor=1,
             int striping_unit=0);

   void open(const std::string& filename, std::ios_base::openmode mode,
             const std::string& filetype_id="", char endianness='N',
             bool turn_on_checksum=false, int striping_factor=1,
             int striping_unit=0);

   void openReadOnly(const std::string& filename, const std::string& filetype_id="",
                     bool turn_on_checksum=false);

   void openNew(const std::string& filename, bool fail_if_exists=true,
                const std::string& filetype_id="", char endianness='N',
                bool turn_on_checksum=false, int striping_factor=1,
                int striping_unit=0);

   void openUpdate(const std::string& filename, const std::string& filetype_id="", 
                   char endianness='N', bool turn_on_checksum=false,
                   int striping_factor=1, int striping_unit=0);

   ~IOFMHandler();

          // closes current file if open, otherwise no action taken
   void close();


          // informational routines
          
   bool isOpen() const { return openflag; }
   
   bool isNewFile() const { return is_new_file; }

   std::string getFileName() const { return m_filename; }
   
   bool isChecksumOn() const { return checksum_on; }
   
   bool isEndianConversionOn() const { return endian_convert; }
   
   bool isFileLittleEndian() const { return (endian_format=='L'); }
   
   bool isFileBigEndian() const { return (endian_format=='B'); }

   bool isReadOnly() const { return read_only; }

   bool isGlobal() const { return global_mode; }
   
   bool isLocal() const { return !global_mode; }
   

#if defined(ARCH_SERIAL)
   typedef std::iostream::pos_type   pos_type;  // position in buffer
   typedef std::iostream::off_type   off_type;  // offset in buffer
   typedef std::iostream::seekdir    whence_type;
   typedef std::iostream::openmode   openmode_type;
#elif defined(ARCH_PARALLEL)
   typedef MPI_Offset   pos_type;    // position in buffer
   typedef MPI_Offset   off_type;    // offset in buffer
   typedef int          whence_type;
   typedef int          openmode_type;
#endif

         //  Set the file pointer relative to start of data in file,
         //  end of file (use negative to go backward), or current 
         //  location.  Checksum is reset if in use.

         //  Caution: make sure to convert sizeof(...) quantities
         //  to an IOFMHandler::off_type(sizeof(...)) if you wish to use 
         //  negative offsets.  sizeof(...) returns an unsigned 
         //  integer type.

   void seekFromStart(off_type offset);  // start means start of data (not file)
   void seekFromCurr(off_type offset);
   void seekFromEnd(off_type offset);
    
   void seek(pos_type offset);          // from start of data
   void seekBegin(off_type offset);     // from start of data
   void seekRelative(off_type offset);
   void seekEnd(off_type offset);
   void rewind();   // puts pointer at start of data in file

         //  Getting the file pointer location in bytes from start of data

   pos_type tell();
   pos_type currentPosition();


         //  Check summing
   
   void turnOnChecksum();
   void turnOffChecksum();
   void resetChecksum();   
   ByteHandler::n_uint32_t getChecksum();


         // print out MPI-IO file hints or ID on primary rank

   void printFileHints();
   void printFileID();

     // Open file, read string at a particular position, then close. Returns
     // true if the file exists and can be opened and its file type matches 
     // "filetype_id", returns false otherwise.  The string is returned in 
     // "stringvalue".  This routine is used by objects in data_io_handler.h 
     // when the multi-file handlers build up their maps of file keys.
           
   bool peekString(std::string& stringvalue, size_t position, 
                   const std::string& filename,
                   const std::string& filetype_id="");

         // write routines  

   void write(const std::string& output);

   template<typename T>
   void write(const T& output);    // basic types

   template<typename T>
   void multi_write(const T* output, int n);    // basic types

   template <typename T>
   void write(const std::vector<T>& output);   // T is basic type

   template <typename T>
   void write(const Array<T>& output);

   void write(const LattField& output);

   void write(const std::vector<LattField>& output);

          // read routines
 
   void read(std::string& input, bool broadcast=true); 

   template<typename T>
   void read(T& input);    // basic types

   template<typename T>
   void multi_read(T* input, int n);    // basic types

   template <typename T>
   void read(std::vector<T>& input);   // basic types

   template <typename T>
   void read(Array<T>& input);   // basic types

   void read(LattField& input);

   void read(std::vector<LattField>& input);

         //  number of bytes routines (will be useful by IOMap for
         //  verifying read/writes and determining whether overwrites
         //  should occur)
  
   size_t numbytes(const std::string& output);

   template<typename T>
   size_t numbytes(const T& output);    // basic types

   template <typename T>
   size_t numbytes(const std::vector<T>& output);  // T is basic type

   template <typename T>
   size_t numbytes(const Array<T>& output);   // basic types

   size_t numbytes(const LattField& output);

   size_t numbytes(const std::vector<LattField>& output);

   void reset_lattcmplxbytes() 
    {lattcmplxbytes=LattField::get_cpu_prec_bytes();}


 private:

      // private utility routines
       
   void open_existing_file(const std::string& filetype_id,
                           openmode_type access_mode);

   void open_new_file(const std::string& filetype_id,
                      char endianness,
                      int striping_factor, int striping_unit);

   void clear();

   bool file_exists();

   void writeCommon(const char *data, size_t nbytes);

   void writeDistributed(const char *data, const DistArrayViewInfo& dinfo,
                         bool datashift=false);

   void readCommon(char *data, size_t nbytes, bool broadcast=true);

   void readDistributed(char *data, const DistArrayViewInfo& dinfo,
                        bool datashift=false);

   void check_for_failure(int errcode, const std::string& mesg="");

   std::string int_to_string(int intval);

   std::string tidyString(const std::string& str);  

   void readIDstring(char &endianness, std::string& ID_string);

   void writeIDstring(const std::string& ID_string);

   void write_common(const char* output, size_t element_size, size_t nelements);

   void write_lattice(const char* output, size_t bytes_per_word, 
                      const DistArrayViewInfo& dinfo);

   void read_common(char* output, size_t element_size, size_t nelements);
  
   void read_lattice(char* input, size_t bytes_per_word,
                     const DistArrayViewInfo& dinfo);

   void compute_lattice_checksum(char* data, size_t bytes_per_site, 
                                 size_t lexicostart, size_t nsites);

   void delete_file();

   void file_open(openmode_type access_mode, const std::string& errmsg);

   void file_close();

   void file_seek(off_type offset, whence_type whence);

   pos_type file_tell();

   void set_striping(int striping_factor, int striping_unit);

   bool peeker(std::string& stringvalue, size_t position, 
               const std::string& filename,
               const std::string& filetype_id="");
   
       // constants
       
   static const int ID_string_length;
   static const pos_type data_start_pos;

   static const int IO_ERR_NO_SUCH_FILE;
   static const int IO_ERR_ACCESS;
   static const int IO_ERR_OTHER;
   static const int IO_SUCCESS;

   static const openmode_type IO_MODE_RDONLY;
   static const openmode_type IO_MODE_RDWR;
   static const openmode_type IO_MODE_CREATE;
   static const openmode_type IO_MODE_EXCL;

   static const whence_type IO_SEEK_BEG;
   static const whence_type IO_SEEK_CUR;
   static const whence_type IO_SEEK_END;

   static const int Nd;
   static ByteHandler m_bytehandler;

         // for static (compile time) assertion
//   template <bool b>
//   void static_assert()
//   { typedef char asserter[b?1:-1]; }

       // purpose of this class is to cause compiler errors if
       // constructed with a non basic type 

    class BasicType
     {
      public:
        explicit BasicType( char s ) {}
        explicit BasicType( int s ) {}
        explicit BasicType( unsigned int s ) {}
        explicit BasicType( short int s ) {}
        explicit BasicType( unsigned short int s ) {}
        explicit BasicType( long int s ) {}
        explicit BasicType( unsigned long int s ) {}
        explicit BasicType( long long int s ) {}
        explicit BasicType( unsigned long long int s ) {}
        explicit BasicType( float s ) {}
        explicit BasicType( double s ) {}
        explicit BasicType( bool s ) {}
        explicit BasicType( cmplx s) {}
        explicit BasicType( fcmplx s) {}
     };


};


 // ***************************************************************


template <typename T>
void IOFMHandler::write(const T& output)
{ 
 IOFMHandler::BasicType dummy(output);  // cause compiler error if T is not a valid basic type
 write_common((const char*)&output, sizeof(T), 1);
}


template <typename T>
void IOFMHandler::multi_write(const T* output, int n)
{ 
 IOFMHandler::BasicType dummy(*output);  // cause compiler error if T is not a valid basic type
 write_common((const char*)output, sizeof(T), n);
}


template <typename T>
void IOFMHandler::write(const std::vector<T>& output)
{
 IOFMHandler::BasicType dummy(output[0]);  // cause compiler error if T is not a valid basic type
 int n=output.size();
 write(n);
 if (n==0) return;
 write_common((const char*)&output[0],sizeof(T),size_t(n));
}

template <typename T>
void IOFMHandler::write(const Array<T>& output)
{
 unsigned int n=output.numDimensions();
 write(n);
 if (n==0) return;
 write(output.m_sizes);
 write(output.m_store);
}

  // *************************************************

  
template<typename T>
void write(IOFMHandler& io, const T& output)   // basic types
 { io.write(output); }

template <typename T>
void write(IOFMHandler& io, const std::vector<T>& output)
 { io.write(output); }

template <typename T>
void write(IOFMHandler& io, const Array<T>& output)  // T is basic type
 { io.write(output); }

inline void write(IOFMHandler& io, const LattField& output)
 { io.write(output); }

inline void write(IOFMHandler& io, const std::vector<LattField>& output)
 { io.write(output); }

// ************************************************************  
  

template <typename T>
void IOFMHandler::read(T& input)
{ 
 IOFMHandler::BasicType dummy(input);  // cause compiler error if T is not a valid basic type
 read_common((char*)&input, sizeof(T), 1);
}


template <typename T>
void IOFMHandler::multi_read(T* input, int n)
{ 
 IOFMHandler::BasicType dummy(*input);  // cause compiler error if T is not a valid basic type
 read_common((char*)input, sizeof(T), n);
}


template <typename T>
void IOFMHandler::read(std::vector<T>& input)
{
 IOFMHandler::BasicType dummy(input[0]);  // cause compiler error if T is not a valid basic type
 int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure((n<0)||(n>16777216), "std::vector too large...bad file location?");
 input.resize(n);
 read_common((char*)&input[0],sizeof(T),size_t(n));
}

template <typename T>
void IOFMHandler::read(Array<T>& input)
{
 unsigned int n;
 read(n);
 if (n==0) return;
    // if reading in wrong location, could get nonsense here,
    // so place a reasonable limit
 check_for_failure(n>128, "Array number of dimensions too large...bad file location?");
 read(input.m_sizes);
 read(input.m_store);
}

// **********************************************

template<typename T>
void read(IOFMHandler& io, T& input)   // basic types
 { io.read(input); }

template <typename T>
void read(IOFMHandler& io, std::vector<T>& input)
 { io.read(input); }

template <typename T>
void read(IOFMHandler& io, Array<T>& input)  // T is basic type
 { io.read(input); }

inline void read(IOFMHandler& io, LattField& input)
 { io.read(input); }

inline void read(IOFMHandler& io, std::vector<LattField>& input)
 { io.read(input); }


// ************************************************************  


template <typename T>
size_t IOFMHandler::numbytes(const T& data)
{ 
 IOFMHandler::BasicType dummy(data);  // cause compiler error if T is not a valid basic type
 return sizeof(T);
}

template <typename T>
size_t IOFMHandler::numbytes(const std::vector<T>& data)
{
 IOFMHandler::BasicType dummy(data[0]);  // cause compiler error if T is not a valid basic type
 return sizeof(int)+data.size()*sizeof(T);
}

template <typename T>
size_t IOFMHandler::numbytes(const Array<T>& data)
{
 return sizeof(unsigned int)+numbytes(data.m_sizes)+numbytes(data.m_store);
}

// ************************************************************  


template<typename T>
size_t numbytes(IOFMHandler& io, const T& data)   // generic routine
 { return data.numbytes(); }

template<> inline size_t numbytes(IOFMHandler& io, const char& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const unsigned int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const short int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const unsigned short int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const long int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const unsigned long int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const long long int& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const unsigned long long int& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const float& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const double& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const bool& data)  
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const cmplx& data) 
{ return io.numbytes(data); }

template<> inline size_t numbytes(IOFMHandler& io, const fcmplx& data) 
{ return io.numbytes(data); }

template <typename T>
size_t numbytes(IOFMHandler& io, const std::vector<T>& data)
 { return io.numbytes(data); }

template <typename T>
size_t numbytes(IOFMHandler& io, const Array<T>& data)
 { return io.numbytes(data); }

inline size_t numbytes(IOFMHandler& io, const LattField& data)
 { return io.numbytes(data); }

inline size_t numbytes(IOFMHandler& io, const std::vector<LattField>& data)
 { return io.numbytes(data); }


// **************************************************************
}
#endif
