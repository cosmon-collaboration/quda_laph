#include "gauge_field_hdf5.h"
#include "layout_info.h"
#include "laph_stdio.h"
#include <hdf5.h>
#include <algorithm>
#include <cctype>
#include <cstring>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#ifndef H5_HAVE_PARALLEL
#error "A parallel build (ARCH_PARALLEL) needs an HDF5 library built with parallel (MPI) support"
#endif
#endif

using namespace std;

namespace LaphEnv {


bool isHDF5FileName(const string& filename)
{
 string lower(filename);
 std::transform(lower.begin(),lower.end(),lower.begin(),
                [](unsigned char c){return std::tolower(c);});
 auto ends_with=[&lower](const string& suffix){
    return (lower.size()>suffix.size())
        && (lower.compare(lower.size()-suffix.size(),suffix.size(),suffix)==0);};
 return ends_with(".hdf5") || ends_with(".h5");
}


namespace {

    // Aborts (on all ranks) with a message if an HDF5 call failed
    // (HDF5 returns a negative hid_t/herr_t on failure).

void h5check(long long status, const string& filename, const char* what)
{
 if (status<0){
    errorLaph(make_strf("HDF5 failure in writeGaugeFieldHDF5 (%s) for file %s",
                        what,filename));}
}

}


void writeGaugeFieldHDF5(const string& filename, const vector<LattField>& gauge,
                         const string& header_xml)
{
 const int Ndim=LayoutInfo::Ndim;
 const hsize_t ncm=FieldNcolor*FieldNcolor;   // complex numbers per color matrix
 if (int(gauge.size())!=Ndim){
    errorLaph("writeGaugeFieldHDF5: gauge field must have one LattField per direction");}
 for (int mu=0;mu<Ndim;++mu){
    if (gauge[mu].getFieldSiteType()!=FieldSiteType::ColorMatrix){
       errorLaph("writeGaugeFieldHDF5: gauge field LattFields must be ColorMatrix fields");}}

    // global extents, local extents, and global coordinates of the local origin
 const vector<int>& gext=LayoutInfo::getLattExtents();
 const vector<int>& lext=LayoutInfo::getRankLattExtents();
 const vector<int>& comm=LayoutInfo::getMyCommCoords();
 vector<hsize_t> N(Ndim),n(Ndim),start(Ndim);
 for (int k=0;k<Ndim;++k){
    N[k]=gext[k]; n[k]=lext[k]; start[k]=hsize_t(comm[k])*hsize_t(lext[k]);}
 const hsize_t loc_nsites=LayoutInfo::getRankLatticeNumSites();
 const int start_parity=LayoutInfo::getMyStartParity();

    // host precision of each direction's field
 vector<bool> dprec(Ndim);
 for (int mu=0;mu<Ndim;++mu){
    const size_t cmplx_bytes=gauge[mu].bytesPerSite()/gauge[mu].elemsPerSite();
    dprec[mu]=(cmplx_bytes==2*sizeof(double));
    if ((!dprec[mu])&&(cmplx_bytes!=2*sizeof(float))){
       errorLaph("writeGaugeFieldHDF5: unsupported precision of gauge field");}}

    // create the file (collectively in a parallel build)
 hid_t fapl=H5Pcreate(H5P_FILE_ACCESS);
 h5check(fapl,filename,"create file access property list");
#ifdef ARCH_PARALLEL
 h5check(H5Pset_fapl_mpio(fapl,MPI_COMM_WORLD,MPI_INFO_NULL),filename,"set MPI-IO file access");
 h5check(H5Pset_all_coll_metadata_ops(fapl,true),filename,"set collective metadata reads");
 h5check(H5Pset_coll_metadata_write(fapl,true),filename,"set collective metadata writes");
#endif
 hid_t file=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,fapl);
 h5check(file,filename,"create file");
 hid_t group=H5Gcreate2(file,"/GaugeField",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 h5check(group,filename,"create group /GaugeField");

    // XML header as a scalar variable-length string attribute
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

    // Per timeslice: file dataspace of all (mu,x,y,z,c1,c2) entries, and the
    // part of it this rank writes on the timeslices it holds: for each local
    // (mu,x), a regular hyperslab over the local y values, each block being
    // the local z values times the color matrix entries.  The local data are
    // packed in lexicographic (mu,x,y,z,c1,c2) order (local coordinates); since
    // the file offset of each entry increases monotonically in this order, it
    // is also the order in which HDF5 maps the buffer onto the selection.
 const hsize_t slice_total=hsize_t(Ndim)*N[0]*N[1]*N[2]*ncm;
 const hsize_t slice_local=hsize_t(Ndim)*n[0]*n[1]*n[2]*ncm;
 hid_t fspace_mine=H5Screate_simple(1,&slice_total,nullptr);
 h5check(fspace_mine,filename,"create file dataspace");
 bool first=true;
 for (hsize_t mu=0;mu<hsize_t(Ndim);++mu)
 for (hsize_t x=0;x<n[0];++x){
    hsize_t offset=(((mu*N[0]+start[0]+x)*N[1]+start[1])*N[2]+start[2])*ncm;
    hsize_t stride=N[2]*ncm, count=n[1], block=n[2]*ncm;
    h5check(H5Sselect_hyperslab(fspace_mine,first?H5S_SELECT_SET:H5S_SELECT_OR,
                                &offset,&stride,&count,&block),filename,"select hyperslab");
    first=false;}
 if (hsize_t(H5Sget_select_npoints(fspace_mine))!=slice_local){
    errorLaph("writeGaugeFieldHDF5: file selection does not match local data size");}
 hid_t mspace_mine=H5Screate_simple(1,&slice_local,nullptr);
 h5check(mspace_mine,filename,"create memory dataspace");
    // empty selections, for timeslices held by other ranks (H5Dwrite is collective)
 hid_t fspace_none=H5Scopy(fspace_mine);
 h5check(fspace_none,filename,"copy file dataspace");
 h5check(H5Sselect_none(fspace_none),filename,"select none in file");
 hid_t mspace_none=H5Scopy(mspace_mine);
 h5check(mspace_none,filename,"copy memory dataspace");
 h5check(H5Sselect_none(mspace_none),filename,"select none in memory");

 hid_t dxpl=H5Pcreate(H5P_DATASET_XFER);
 h5check(dxpl,filename,"create transfer property list");
#ifdef ARCH_PARALLEL
 h5check(H5Pset_dxpl_mpio(dxpl,H5FD_MPIO_COLLECTIVE),filename,"set collective transfer");
#endif
 hid_t dcpl=H5Pcreate(H5P_DATASET_CREATE);
 h5check(dcpl,filename,"create dataset creation property list");
 h5check(H5Pset_fill_time(dcpl,H5D_FILL_TIME_NEVER),filename,"set fill time");

 vector<double> d_re(slice_local), d_im(slice_local);
 for (hsize_t tg=0;tg<N[3];++tg){
    const bool mine=(tg>=start[3])&&(tg<start[3]+n[3]);
    if (mine){
          // reorder from the even-odd checkerboard layout of LattField,
          // converting to double precision
       const hsize_t t=tg-start[3];
       hsize_t k=0;
       for (int mu=0;mu<Ndim;++mu){
          const char* data=gauge[mu].getDataConstPtr();
          for (hsize_t x=0;x<n[0];++x)
          for (hsize_t y=0;y<n[1];++y)
          for (hsize_t z=0;z<n[2];++z){
             const hsize_t site=(x+n[0]*(y+n[1]*(z+n[2]*t)))/2
                               +(loc_nsites/2)*((start_parity+x+y+z+t)%2);
             if (dprec[mu]){
                const double* src=reinterpret_cast<const double*>(data)+2*ncm*site;
                for (hsize_t c=0;c<ncm;++c,++k){
                   d_re[k]=src[2*c]; d_im[k]=src[2*c+1];}}
             else{
                const float* src=reinterpret_cast<const float*>(data)+2*ncm*site;
                for (hsize_t c=0;c<ncm;++c,++k){
                   d_re[k]=src[2*c]; d_im[k]=src[2*c+1];}}}}}

    const string tname=make_strf("t%d",int(tg));
    hid_t tgroup=H5Gcreate2(group,tname.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    h5check(tgroup,filename,"create timeslice group");
    const char* names[2]={"re","im"};
    const double* bufs[2]={d_re.data(),d_im.data()};
    for (int part=0;part<2;++part){
       hid_t dset=H5Dcreate2(tgroup,names[part],H5T_NATIVE_DOUBLE,fspace_mine,
                             H5P_DEFAULT,dcpl,H5P_DEFAULT);
       h5check(dset,filename,"create dataset");
       h5check(H5Dwrite(dset,H5T_NATIVE_DOUBLE,mine?mspace_mine:mspace_none,
                        mine?fspace_mine:fspace_none,dxpl,bufs[part]),
               filename,"write dataset");
       h5check(H5Dclose(dset),filename,"close dataset");}
    h5check(H5Gclose(tgroup),filename,"close timeslice group");}

 H5Pclose(dcpl); H5Pclose(dxpl);
 H5Sclose(mspace_none); H5Sclose(fspace_none);
 H5Sclose(mspace_mine); H5Sclose(fspace_mine);
 h5check(H5Gclose(group),filename,"close group");
 h5check(H5Fclose(file),filename,"close file");
 H5Pclose(fapl);
}

}
