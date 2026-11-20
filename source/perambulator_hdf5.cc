#include "perambulator_hdf5.h"
#include "laph_stdio.h"
#include "utils.h"
#include <hdf5.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

using namespace std;

namespace LaphEnv {

namespace {

void h5check(long long status, const string& filename, const char* what)
{
 if (status<0){
    errorLaph(make_strf("HDF5 failure in PerambulatorHDF5Writer (%s) for file %s",
                        what,filename));}
}

   // the same value on all ranks as on the primary rank
void broadcastBool(bool& flag)
{
#ifdef ARCH_PARALLEL
 int iflag=flag?1:0;
 comm_broadcast(&iflag,sizeof(int),0);
 flag=(iflag!=0);
#endif
}

string matrixPath(int src_time, int src_spin, int snk_time, int snk_spin)
{
 return make_strf("/PerambulatorData/srcTime%d_snkTime%d/srcSpin%d/snkSpin%d",
                  src_time,snk_time,src_spin,snk_spin);
}

   // true if every link along the absolute "path" exists
bool pathExists(hid_t file, const string& path)
{
 size_t pos=0;
 while (pos!=string::npos){
    pos=path.find('/',pos+1);
    string prefix=path.substr(0,pos);
    if (H5Lexists(file,prefix.c_str(),H5P_DEFAULT)<=0) return false;}
 return true;
}

}


string PerambulatorHDF5Writer::open(const string& filename, const string& header_xml,
                                    int nev, bool overwrite)
{
 close();
 m_filename=filename;
 m_nev=nev;
 string existing_header;
 if (isPrimaryRank()){
       // a sieve buffer holding a whole matrix, so that writing or reading
       // one (strided) column is one read-modify-write of the file
    hid_t fapl=H5Pcreate(H5P_FILE_ACCESS);
    h5check(fapl,filename,"create file access property list");
    h5check(H5Pset_sieve_buf_size(fapl,max(size_t(1)<<20,size_t(nev)*size_t(nev)*sizeof(double))),
            filename,"set sieve buffer size");
    bool create=overwrite||(!fileExists(filename,false));
    if (create){
       m_file=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,fapl);
       h5check(m_file,filename,"create file");
       hid_t group=H5Gcreate2(m_file,"/PerambulatorData",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
       h5check(group,filename,"create group /PerambulatorData");
       hid_t strtype=H5Tcopy(H5T_C_S1);
       h5check(strtype,filename,"copy string type");
       h5check(H5Tset_size(strtype,H5T_VARIABLE),filename,"set string type size");
       hid_t attspace=H5Screate(H5S_SCALAR);
       h5check(attspace,filename,"create attribute dataspace");
       hid_t attr=H5Acreate2(group,"Header",strtype,attspace,H5P_DEFAULT,H5P_DEFAULT);
       h5check(attr,filename,"create Header attribute");
       const char* hdr=header_xml.c_str();
       h5check(H5Awrite(attr,strtype,&hdr),filename,"write Header attribute");
       H5Aclose(attr); H5Sclose(attspace); H5Tclose(strtype);
       h5check(H5Gclose(group),filename,"close group");}
    else{
       m_file=H5Fopen(filename.c_str(),H5F_ACC_RDWR,fapl);
       h5check(m_file,filename,"open existing file");
       hid_t attr=H5Aopen_by_name(m_file,"/PerambulatorData","Header",H5P_DEFAULT,H5P_DEFAULT);
       h5check(attr,filename,"open Header attribute of existing file");
       hid_t strtype=H5Tcopy(H5T_C_S1);
       h5check(strtype,filename,"copy string type");
       h5check(H5Tset_size(strtype,H5T_VARIABLE),filename,"set string type size");
       char* hdr=nullptr;
       h5check(H5Aread(attr,strtype,&hdr),filename,"read Header attribute");
       existing_header=(hdr)?string(hdr):string(" ");   // non-empty: file existed
       H5free_memory(hdr);
       H5Tclose(strtype); H5Aclose(attr);}
    H5Pclose(fapl);}
#ifdef ARCH_PARALLEL
 broadcastString(existing_header);
#endif
 return existing_header;
}


bool PerambulatorHDF5Writer::queryColumn(int src_time, int src_spin, const vector<int>& snk_times,
                                         int nspin, int srcev)
{
 if (!isOpen()) errorLaph("PerambulatorHDF5Writer::queryColumn: no file open");
 bool complete=true;
 if (isPrimaryRank()){
    const hsize_t nev=m_nev;
    const hsize_t total=nev*nev;
    hid_t fspace=H5Screate_simple(1,&total,nullptr);
    h5check(fspace,m_filename,"create file dataspace");
    hsize_t start=srcev, stride=nev, count=nev;
    h5check(H5Sselect_hyperslab(fspace,H5S_SELECT_SET,&start,&stride,&count,nullptr),
            m_filename,"select column");
    hid_t mspace=H5Screate_simple(1,&count,nullptr);
    h5check(mspace,m_filename,"create memory dataspace");
    vector<double> col(nev);
    for (size_t it=0;(it<snk_times.size())&&complete;++it)
    for (int snk_spin=1;(snk_spin<=nspin)&&complete;++snk_spin){
       string path=matrixPath(src_time,src_spin,snk_times[it],snk_spin);
       for (const char* part : {"/re","/im"}){
          if (!complete) break;
          string dpath=path+part;
          if (!pathExists(m_file,dpath)){
             complete=false; break;}
          hid_t dset=H5Dopen2(m_file,dpath.c_str(),H5P_DEFAULT);
          h5check(dset,m_filename,"open dataset");
          h5check(H5Dread(dset,H5T_NATIVE_DOUBLE,mspace,fspace,H5P_DEFAULT,col.data()),
                  m_filename,"read column");
          H5Dclose(dset);
          for (double v : col){
             if (std::isnan(v)){
                complete=false; break;}}}}
    H5Sclose(mspace); H5Sclose(fspace);}
 broadcastBool(complete);
 return complete;
}


void PerambulatorHDF5Writer::putColumns(int src_time, int src_spin, int snk_time, int snk_spin,
                                        const vector<int>& srcevs,
                                        const vector<complex<double>>& data)
{
 if (!isOpen()) errorLaph("PerambulatorHDF5Writer::putColumns: no file open");
 if (!isPrimaryRank()) return;
 const hsize_t nev=m_nev;
 const hsize_t ncol=srcevs.size();
 if (ncol==0) return;
 if (data.size()!=nev*ncol){
    errorLaph("PerambulatorHDF5Writer::putColumns: data size does not match");}

    // HDF5 maps the buffer onto the selection in increasing file offset
    // n*nev+m, so the columns are put in increasing order
 vector<size_t> perm(ncol);
 std::iota(perm.begin(),perm.end(),0);
 std::sort(perm.begin(),perm.end(),[&srcevs](size_t a,size_t b){return srcevs[a]<srcevs[b];});
 const hsize_t total=nev*nev;
 hid_t fspace=H5Screate_simple(1,&total,nullptr);
 h5check(fspace,m_filename,"create file dataspace");
 for (hsize_t k=0;k<ncol;++k){
    int m=srcevs[perm[k]];
    if ((m<0)||(hsize_t(m)>=nev)||((k>0)&&(m==srcevs[perm[k-1]]))){
       errorLaph("PerambulatorHDF5Writer::putColumns: invalid source eigvec indices");}
    hsize_t start=m, stride=nev, count=nev;
    h5check(H5Sselect_hyperslab(fspace,(k==0)?H5S_SELECT_SET:H5S_SELECT_OR,
                                &start,&stride,&count,nullptr),m_filename,"select columns");}
 const hsize_t nbuf=nev*ncol;
 hid_t mspace=H5Screate_simple(1,&nbuf,nullptr);
 h5check(mspace,m_filename,"create memory dataspace");
 vector<double> re(nbuf), im(nbuf);
 for (hsize_t n=0;n<nev;++n)
 for (hsize_t k=0;k<ncol;++k){
    const complex<double>& z=data[n*ncol+perm[k]];
    re[n*ncol+k]=z.real(); im[n*ncol+k]=z.imag();}

    // new datasets: contiguous, allocated at creation and filled with NaN
 hid_t dcpl=H5Pcreate(H5P_DATASET_CREATE);
 h5check(dcpl,m_filename,"create dataset creation property list");
 const double nan=std::numeric_limits<double>::quiet_NaN();
 h5check(H5Pset_fill_value(dcpl,H5T_NATIVE_DOUBLE,&nan),m_filename,"set fill value");
 h5check(H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY),m_filename,"set allocation time");
 h5check(H5Pset_fill_time(dcpl,H5D_FILL_TIME_ALLOC),m_filename,"set fill time");
 hid_t lcpl=H5Pcreate(H5P_LINK_CREATE);
 h5check(lcpl,m_filename,"create link creation property list");
 h5check(H5Pset_create_intermediate_group(lcpl,1),m_filename,"set intermediate group creation");

 string path=matrixPath(src_time,src_spin,snk_time,snk_spin);
 const char* parts[2]={"/re","/im"};
 const double* bufs[2]={re.data(),im.data()};
 for (int p=0;p<2;++p){
    string dpath=path+parts[p];
    hid_t dset;
    if (pathExists(m_file,dpath)){
       dset=H5Dopen2(m_file,dpath.c_str(),H5P_DEFAULT);
       h5check(dset,m_filename,"open dataset");}
    else{
       hid_t space=H5Screate_simple(1,&total,nullptr);
       h5check(space,m_filename,"create dataspace");
       dset=H5Dcreate2(m_file,dpath.c_str(),H5T_NATIVE_DOUBLE,space,lcpl,dcpl,H5P_DEFAULT);
       h5check(dset,m_filename,"create dataset");
       H5Sclose(space);}
    h5check(H5Dwrite(dset,H5T_NATIVE_DOUBLE,mspace,fspace,H5P_DEFAULT,bufs[p]),
            m_filename,"write columns");
    h5check(H5Dclose(dset),m_filename,"close dataset");}
 H5Pclose(lcpl); H5Pclose(dcpl);
 H5Sclose(mspace); H5Sclose(fspace);
}


void PerambulatorHDF5Writer::flush()
{
 if (isPrimaryRank()&&(m_file>=0)){
    h5check(H5Fflush(m_file,H5F_SCOPE_GLOBAL),m_filename,"flush file");}
}


void PerambulatorHDF5Writer::close()
{
 if (isPrimaryRank()&&(m_file>=0)){
    h5check(H5Fclose(m_file),m_filename,"close file");}
 m_file=-1;
 m_filename.clear();
}

}
