#ifdef TESTING
#include <string>
#include "task_tests.h"
#include "xml_handler.h"
#include "laph_stdio.h"
#include "latt_field.h"
#include "layout_info.h"
#include "io_handler_fm.h"
#include "io_map.h"
#include "utils.h"
#include "gauge_configuration_info.h"
#include "filelist_info.h"
#include "data_io_handler.h"
#include "laph_noise_info.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdarg.h>

using namespace std;
using namespace LaphEnv;

namespace QLTestEnv {

// ************************************************

bool do_write_check(IOHandlerFM& ioh, int expected_pos)
{
 bool flag=true;
 int pos=ioh.currentPosition();
 if (!isIntegerSameAllRanks(pos)){
    printLaph("Write error: not all ranks see same position afterwards");
    flag=false;}
 if (pos!=expected_pos){
    printLaph(make_strf("Write error: position afterwards is not as expected: got %d expected %d",pos,expected_pos));
    flag=false;}
 if (ioh.isChecksumOn()){
    ByteHandler::n_uint32_t csum=ioh.getChecksum();
    if (!isIntegerSameAllRanks(int(csum))){
      printLaph("Write error: not all ranks see same checksum");
      flag=false;}
    printLaph(make_strf("Checksum after write is %u",int(csum)));}
 return flag;
}

bool do_read_check(IOHandlerFM& ioh, int expected_pos, bool value_correct)
{
 bool flag=true;
 int pos=ioh.currentPosition();
 if (!isIntegerSameAllRanks(pos)){
    printLaph("Read error: not all ranks see same position afterwards");
    flag=false;}
 if (pos!=expected_pos){
    printLaph(make_strf("Read error: position afterwards is not as expected: got %d expected %d",pos,expected_pos));
    flag=false;}
 if (!value_correct){
    printLaph("Read error: incorrect value(s) read");
    flag=false;}
 if (ioh.isChecksumOn()){
    ByteHandler::n_uint32_t csum=ioh.getChecksum();
    if (!isIntegerSameAllRanks(int(csum))){
      printLaph("Read error: not all ranks see same checksum");
      flag=false;}
    printLaph(make_strf("Checksum after read is %u",int(csum)));}
 return flag;
}

std::vector<int> make_int_vector(int a, int b, int c, int d)
{
 vector<int> res(4); res[0]=a; res[1]=b; res[2]=c; res[3]=d;
 return res;
}



void testIOHandlerFM(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestIOHandlerFM")==0)
 return;

 bool success=true;
 printLaph("Running TestIOHandlerFM");
 XMLHandler xmlr(xml_in,"TestIOHandlerFM");
 
 string testfilename;
 xmlread(xmlr,"TestFileName",testfilename,"testIOHandlerFM");
 testfilename=tidyString(testfilename);
 if (testfilename.empty()){
    errorLaph("Empty test file name");}
 printLaph(make_strf("Test file name is %s",testfilename));

 bool turn_on_checksum=false;
 if (xml_tag_count(xmlr,"ChecksumsOn")>0)
    turn_on_checksum=true;
 char endianness='N';
 if (xml_tag_count(xmlr,"LittleEndian")>0)
    endianness='L';
 else if (xml_tag_count(xmlr,"BigEndian")>0)
    endianness='B';
 uint MTseed1=0;
 xmlreadif(xmlr,"MTSeed1",MTseed1,"TestIOHandlerFM");
 uint MTseed2=0;
 xmlreadif(xmlr,"MTSeed2",MTseed2,"TestIOHandlerFM");
 uint MTseed3=0;
 xmlreadif(xmlr,"MTSeed3",MTseed3,"TestIOHandlerFM");
 uint MTseed4=0;
 xmlreadif(xmlr,"MTSeed4",MTseed4,"TestIOHandlerFM");
 double tol=-0.1;
 xmlreadif(xmlr,"Tol",tol,"TestIOHandlerFM");

 if (turn_on_checksum) printLaph("Checksums are ON");
 else printLaph("Checksums are OFF");
 printLaph(make_strf("Endianness = %c\n",endianness));
 
 int k=5;
 float x=7.4;
 string str("ABC");
 vector<float> g(3);
 for (int q=0;q<3;++q) g[q]=q;
 std::complex<double> z(-2.2,1.1);
 Array<float> h(2,3);
 float xx=0.34;
 for (int jj=0;jj<3;++jj)
 for (int ii=0;ii<2;++ii){
    h(ii,jj)=xx; xx+=1.0;}

 vector<std::complex<float>> zvec(45);
 for (uint k=0;k<45;k++) zvec[k]=complex<float>(k*3.2,(k-1)*0.452);

     // Make a lattice field of complex numbers, color vectors, etc.
 LatticeAssigner QLA(MTseed1);
 LattField LC(FieldSiteType::Complex);
 QLA.assign_field(LC,"LatCmplx");
 QLA.reSeed(MTseed2);
 LattField LV(FieldSiteType::ColorVector);
 QLA.assign_field(LV,"LatColorVec");
 QLA.reSeed(MTseed4);
 LattField LSV(FieldSiteType::ColorSpinVector);
 QLA.assign_field(LSV,"LatColorSpinVec");

 GaugeFieldAssigner GC(MTseed1,MTseed2,MTseed3,MTseed4);
 vector<LattField> Config(4,FieldSiteType::ColorMatrix);
 GC.assign_gauge_field(Config,"GaugeField");

 vector<vector<int>> checkSites;
 checkSites.push_back(vector<int>{0,0,0,0});
 const vector<int> N(LayoutInfo::getLattSizes());
 checkSites.push_back(make_int_vector(N[0]/3,N[1]/2,N[2]/4,N[3]/2));
 checkSites.push_back(make_int_vector(N[0]/2,2*N[1]/3,N[2]/6,N[3]/5));
 checkSites.push_back(make_int_vector(N[0]/4,3*N[1]/5,3*N[2]/4,N[3]/4));
 checkSites.push_back(make_int_vector(N[0]-1,N[1]-1,N[2]-1,N[3]-1));

 if (!fileExists(testfilename)){
    printLaph(make_strf("Writing to IOHandlerFM using file %s\n",testfilename));
    IOHandlerFM ioh(testfilename,IOHandlerFM::ReadWriteEraseIfExists,"IOTester",
                    endianness,1,0,turn_on_checksum);
    bool check1=ioh.isChecksumOn();
    if (check1) printLaph("Write IOHandlerFM has checksums turned on");
    else printLaph("Write IOHandlerFM has checksums turned off");
    if (!isBooleanSameAllRanks(check1)){ success=false; printLaph("ERROR not all checksum states agree");}
    printLaph(make_strf("Is endian conversion on? %d",ioh.isEndianConversionOn()));
    printLaph(make_strf("Is file little endian? %d",ioh.isFileLittleEndian()));
    printLaph(make_strf("Is file big endian? %d",ioh.isFileBigEndian()));
    printLaph(make_strf("Is file read only? %d",ioh.isReadOnly()));
    printLaph(make_strf("Is IOHandler global? %d",ioh.isGlobal()));
    printLaph(make_strf("Is IOHandler local? %d",ioh.isLocal()));

    size_t loc=0;
    printLaph("Writing integer");
    write(ioh,k); loc+=sizeof(int);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing float");
    write(ioh,x); loc+=sizeof(float);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing string of 3 char");
    write(ioh,str); loc+=sizeof(int)+str.length();
    success=success && do_write_check(ioh,loc);
    printLaph("Writing vector of 3 floats");
    write(ioh,g); loc+=sizeof(int)+g.size()*sizeof(float);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing complex<double>");
    write(ioh,z); loc+=sizeof(z);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing Array(2,3) of floats");
    write(ioh,h); loc+=sizeof(int)*(h.numDimensions()+3)+sizeof(float)*h.size();
    success=success && do_write_check(ioh,loc);
    printLaph("Writing vector of 45 complex<float>");
    write(ioh,zvec); loc+=zvec.size()*sizeof(complex<float>)+sizeof(int);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing lattice of complex");
    printField(LC,"LatticeComplex",checkSites); 
    write(ioh,LC); loc+=sizeof(int)*(LayoutInfo::Ndim+2)+LayoutInfo::getLatticeNumSites()*LC.bytesPerSite();
    success=success && do_write_check(ioh,loc);
    printLaph("Writing an integer");
    write(ioh,k+9); loc+=sizeof(int);
    success=success && do_write_check(ioh,loc);
    printLaph("Writing lattice of color vector");
    printField(LV,"LatticeColorVector",checkSites); 
    write(ioh,LV); loc+=sizeof(int)*(LayoutInfo::Ndim+2)+LV.bytesPerSite()*LayoutInfo::getLatticeNumSites();
    success=success && do_write_check(ioh,loc);
    printLaph("Writing lattice of color-spin vector");
    printField(LSV,"LatticeColorSpin",checkSites); 
    write(ioh,LSV); loc+=sizeof(int)*(LayoutInfo::Ndim+2)+LayoutInfo::getLatticeNumSites()*LSV.bytesPerSite();
    success=success && do_write_check(ioh,loc);
    printLaph("Writing vector of lattice of color matrix");
    printField(Config[0],"GaugeField[0]",checkSites); 
    printField(Config[2],"GaugeField[2]",checkSites); 
    write(ioh,Config); loc+=sizeof(int)*(1+LayoutInfo::Ndim*(LayoutInfo::Ndim+2))
                          +LayoutInfo::Ndim*Config[0].bytesPerSite()*LayoutInfo::getLatticeNumSites();
    success=success && do_write_check(ioh,loc);

    if (ioh.isChecksumOn()){
       printLaph(make_strf("checksum is %u",ioh.getChecksum()));}
    }
 else{
    printLaph("Writing tests skipped since file already exists");}

#ifdef ARCH_PARALLEL
 comm_barrier();
#endif

 printLaph(make_strf("Reading from IOHandlerFM using file %s\n",testfilename));

 {IOHandlerFM Q(testfilename,IOHandlerFM::ReadOnly,"IOTester",
                endianness,1,0,turn_on_checksum);
 bool check1=Q.isChecksumOn();
 if (check1) printLaph("Read IOHandlerFM has checksums turned on");
 else printLaph("Read IOHandlerFM has checksums turned off");
 if (!isBooleanSameAllRanks(check1)){ success=false; printLaph("ERROR not all checksum states agree");}
 printLaph(make_strf("Is endian conversion on? %d",Q.isEndianConversionOn()));
 printLaph(make_strf("Is file little endian? %d",Q.isFileLittleEndian()));
 printLaph(make_strf("Is file big endian? %d",Q.isFileBigEndian()));
 printLaph(make_strf("Is file read only? %d",Q.isReadOnly()));
 printLaph(make_strf("Is IOHandler global? %d",Q.isGlobal()));
 printLaph(make_strf("Is IOHandler local? %d",Q.isLocal()));

 size_t loc=0;
 printLaph("Reading integer");
 int kread;  read(Q,kread); loc+=sizeof(int);
 bool valflag=(kread==k);
 success=success && do_read_check(Q,loc,valflag);
 printLaph("Reading float");
 float xread; read(Q,xread); loc+=sizeof(float);
 valflag=(std::abs(xread-x)<1e-10);
 success=success && do_read_check(Q,loc,valflag);
 printLaph("Reading string");
 string strread; read(Q,strread); loc+=sizeof(int)+str.length();
 valflag=(strread==str);
 success=success && do_read_check(Q,loc,valflag);
 printLaph("Reading vector of float");
 vector<float> gread; read(Q,gread); loc+=sizeof(int)+g.size()*sizeof(float);
 valflag=true;
 for (uint q=0;q<g.size();++q) 
    if (g[q]!=gread[q]) {valflag=false;}
 success=success && do_read_check(Q,loc,valflag);
 complex<double> zread; read(Q,zread);  loc+=sizeof(z);
 valflag=(std::abs(zread-z)<1e-10);
 success=success && do_read_check(Q,loc,valflag);
 Array<float> hh;
 printLaph("Reading Array of float");
 read(Q,hh); loc+=sizeof(int)*(h.numDimensions()+3)+sizeof(float)*h.size();
 valflag=true;
 uint hdim=hh.numDimensions();
 if (hdim!=h.numDimensions()) { valflag=false;}
 else{
    for (int jj=0;jj<3;++jj)
    for (int ii=0;ii<2;++ii){
       if (h(ii,jj)!=hh(ii,jj)){ valflag=false;}}}
 success=success && do_read_check(Q,loc,valflag);
 vector<std::complex<float>> zvecread;
 printLaph("Reading vector of complex<float>");
 read(Q,zvecread); loc+=zvec.size()*sizeof(complex<float>)+sizeof(int);
 valflag=true;
 if (zvecread.size()!=45) { valflag=false;}
 for (uint k=0;k<45;k++) 
    if (zvec[k]!=zvecread[k]) { valflag=false;}
 success=success && do_read_check(Q,loc,valflag);
 LattField LCread; 
 printLaph("Reading lattice complex quantity");
 read(Q,LCread);
 if (LCread.getFieldSiteType()!=FieldSiteType::Complex){
    success=false; printLaph("Lattice complex has wrong site type");}

 size_t cmplx_filebytes=LCread.bytesPerSite();
 size_t colvec_filebytes=3*cmplx_filebytes;
 size_t colspinvec_filebytes=12*cmplx_filebytes;
 size_t colmat_filebytes=9*cmplx_filebytes;
 size_t loctry=sizeof(int)*(LayoutInfo::Ndim+2)+LayoutInfo::getLatticeNumSites()*cmplx_filebytes;
 printField(LCread,"LatticeComplex",checkSites); 
 LatticeChecker QLC(QLA);
 QLA.reSeed(MTseed1);
 valflag=true;
 if (!QLC.check_field(LCread,"LatCmplx",true,tol)) { valflag=false;}
 if (size_t(Q.currentPosition())<(loc+loctry)){
    cmplx_filebytes/=2; colvec_filebytes/=2; colspinvec_filebytes/=2; colmat_filebytes/=2;}
 else if (size_t(Q.currentPosition())>(loc+loctry)){
    cmplx_filebytes*=2; colvec_filebytes*=2; colspinvec_filebytes*=2; colmat_filebytes*=2;}
 loc+=sizeof(int)*(LayoutInfo::Ndim+2)+LayoutInfo::getLatticeNumSites()*cmplx_filebytes;
 success=success && do_read_check(Q,loc,valflag);

 printLaph("Reading integer");
 read(Q,kread); loc+=sizeof(int);
 valflag=(kread==(k+9));
 success=success && do_read_check(Q,loc,valflag);

 printLaph("Reading lattice color vector quantity"); 
 LattField LVread;  
 read(Q,LVread);
 if (LVread.getFieldSiteType()!=FieldSiteType::ColorVector){ 
    success=false; printLaph("Lattice color-vector has wrong site type");}
 loc+=sizeof(int)*(LayoutInfo::Ndim+2)+colvec_filebytes*LayoutInfo::getLatticeNumSites();
 printField(LVread,"LatticeColorVector",checkSites); 
 valflag=true;
 QLA.reSeed(MTseed2);
 if (!QLC.check_field(LVread,"LatColorVector",true,tol)) { valflag=false;}
 success=success && do_read_check(Q,loc,valflag);

 printLaph("Reading lattice color-spin vector quantity"); 
 LattField LSVread;  
 read(Q,LSVread); 
 if (LSVread.getFieldSiteType()!=FieldSiteType::ColorSpinVector){
    success=false; printLaph("Lattice color-spin has wrong site type");}
 loc+=sizeof(int)*(LayoutInfo::Ndim+2)+colspinvec_filebytes*LayoutInfo::getLatticeNumSites();
 printField(LSVread,"LatticeColorSpinVector",checkSites); 
 valflag=true;
 QLA.reSeed(MTseed4);
 if (!QLC.check_field(LSVread,"LatSpinColorVector",true,tol)) { valflag=false;}
 success=success && do_read_check(Q,loc,valflag);

 printLaph("Reading lattice gauge field quantity"); 
 vector<LattField> Gread;  
 read(Q,Gread); 
 for (int k=0;k<int(Gread.size());++k){
    if (Gread[k].getFieldSiteType()!=FieldSiteType::ColorMatrix){
       success=false; printLaph("Lattice color matrix has wrong site type");}}
 loc+=sizeof(int)*(1+LayoutInfo::Ndim*(LayoutInfo::Ndim+2))
        +LayoutInfo::Ndim*colmat_filebytes*LayoutInfo::getLatticeNumSites();
 printField(Gread[0],"GaugeField[0]",checkSites); 
 printField(Gread[2],"GaugeField[2]",checkSites); 
 valflag=true;
 GaugeFieldChecker GFC(GC);
 GC.reSeed(MTseed1,MTseed2,MTseed3,MTseed4);
 if (!GFC.check_gauge_field(Gread,"GaugeField",true,tol)) { valflag=false;}
 success=success && do_read_check(Q,loc,valflag);

 if (Q.isChecksumOn()){
    printLaph(make_strf("Checksum is %u",Q.getChecksum()));}
 
 }

 success=globalAnd(success);
 if (success) printLaph("ALL GLOBAL TESTS PASSED!!");
 else         printLaph("Some global tests FAILED");

   //   Now for tests in local mode

 int myrank=LayoutInfo::getMyRank();
 testfilename+="_"+make_string(myrank);
 string logfile("testiohandler_"); logfile+=make_string(myrank)+".log";
 ofstream olog(logfile);

 if (!fileExists(testfilename)){

    olog << "Doing local write to file "<<testfilename<<endl<<endl;
    IOHandlerFM ioh(testfilename,IOHandlerFM::ReadWriteEraseIfExists,"IOTester",
                    endianness,1,0,turn_on_checksum,false);
    bool check1=ioh.isChecksumOn();
    if (check1) olog <<"Write IOHandlerFM has checksums turned on"<<endl;
    else olog <<"Write IOHandlerFM has checksums turned off"<<endl;
    olog <<"Is endian conversion on? "<<ioh.isEndianConversionOn()<<endl;
    olog <<"Is file little endian? "<<ioh.isFileLittleEndian()<<endl;
    olog <<"Is file big endian? "<<ioh.isFileBigEndian()<<endl;
    olog <<"Is file read only? "<<ioh.isReadOnly()<<endl;
    olog <<"Is IOHandler global? "<<ioh.isGlobal()<<endl;
    olog <<"Is IOHandler local? "<<ioh.isLocal()<<endl;

    size_t loc=0;
    olog <<"Writing integer"<<endl;
    write(ioh,k); loc+=sizeof(int);
    success=success && do_write_check(ioh,loc);
    olog <<"Writing float"<<endl;
    write(ioh,x); loc+=sizeof(float);
    success=success && do_write_check(ioh,loc);
    olog <<"Writing string of 3 char"<<endl;
    write(ioh,str); loc+=sizeof(int)+str.length();
    success=success && do_write_check(ioh,loc);
    olog <<"Writing vector of 3 floats"<<endl;
    write(ioh,g); loc+=sizeof(int)+g.size()*sizeof(float);
    success=success && do_write_check(ioh,loc);
    olog <<"Writing complex<double>"<<endl;
    write(ioh,z); loc+=sizeof(z);
    success=success && do_write_check(ioh,loc);
    olog <<"Writing Array(2,3) of floats"<<endl;
    write(ioh,h); loc+=sizeof(int)*(h.numDimensions()+3)+sizeof(float)*h.size();
    success=success && do_write_check(ioh,loc);
    olog <<"Writing vector of 45 complex<float>"<<endl;
    write(ioh,zvec); loc+=zvec.size()*sizeof(complex<float>)+sizeof(int);
    success=success && do_write_check(ioh,loc);
    olog <<"Writing an integer"<<endl;
    write(ioh,k+9); loc+=sizeof(int);
    success=success && do_write_check(ioh,loc);
    if (ioh.isChecksumOn()){
       olog <<"checksum is "<<ioh.getChecksum()<<endl;}
    }
 else{
    olog <<"Writing tests skipped since file already exists"<<endl;}

 olog << "Reading from IOHandlerFM using file "<<testfilename<<endl;

 {IOHandlerFM Q(testfilename,IOHandlerFM::ReadOnly,"IOTester",
                endianness,1,0,turn_on_checksum,false);
 bool check1=Q.isChecksumOn();
 if (check1) olog << "Read IOHandlerFM has checksums turned on"<<endl;
 else olog << "Read IOHandlerFM has checksums turned off"<<endl;
 olog << "Is endian conversion on? "<<Q.isEndianConversionOn()<<endl;
 olog << "Is file little endian? "<<Q.isFileLittleEndian()<<endl;
 olog << "Is file big endian? "<<Q.isFileBigEndian()<<endl;
 olog << "Is file read only? "<<Q.isReadOnly()<<endl;
 olog << "Is IOHandler global? "<<Q.isGlobal()<<endl;
 olog << "Is IOHandler local? "<<Q.isLocal()<<endl;

 size_t loc=0;
 olog << "Reading integer"<<endl;
 int kread;  read(Q,kread); loc+=sizeof(int);
 bool valflag=(kread==k);
 success=success && do_read_check(Q,loc,valflag);
 olog << "Reading float"<<endl;
 float xread; read(Q,xread); loc+=sizeof(float);
 valflag=(std::abs(xread-x)<1e-10);
 success=success && do_read_check(Q,loc,valflag);
 olog << "Reading string"<<endl;
 string strread; read(Q,strread); loc+=sizeof(int)+str.length();
 valflag=(strread==str);
 success=success && do_read_check(Q,loc,valflag);
 olog << "Reading vector of float"<<endl;
 vector<float> gread; read(Q,gread); loc+=sizeof(int)+g.size()*sizeof(float);
 valflag=true;
 for (uint q=0;q<g.size();++q) 
    if (g[q]!=gread[q]) {valflag=false;}
 success=success && do_read_check(Q,loc,valflag);
 complex<double> zread; read(Q,zread);  loc+=sizeof(z);
 valflag=(std::abs(zread-z)<1e-10);
 success=success && do_read_check(Q,loc,valflag);
 Array<float> hh;
 olog << "Reading Array of float"<<endl;
 read(Q,hh); loc+=sizeof(int)*(h.numDimensions()+3)+sizeof(float)*h.size();
 valflag=true;
 uint hdim=hh.numDimensions();
 if (hdim!=h.numDimensions()) { valflag=false;}
 else{
    for (int jj=0;jj<3;++jj)
    for (int ii=0;ii<2;++ii){
       if (h(ii,jj)!=hh(ii,jj)){ valflag=false;}}}
 success=success && do_read_check(Q,loc,valflag);
 vector<std::complex<float>> zvecread;
 olog << "Reading vector of complex<float>"<<endl;
 read(Q,zvecread); loc+=zvec.size()*sizeof(complex<float>)+sizeof(int);
 valflag=true;
 if (zvecread.size()!=45) { valflag=false;}
 for (uint k=0;k<45;k++) 
    if (zvec[k]!=zvecread[k]) { valflag=false;}
 success=success && do_read_check(Q,loc,valflag);
 olog << "Reading integer"<<endl;
 read(Q,kread); loc+=sizeof(int);
 valflag=(kread==(k+9));
 success=success && do_read_check(Q,loc,valflag);

 if (Q.isChecksumOn()){
    olog << "Checksum is "<<Q.getChecksum()<<endl;}}

 olog.close();
 success=globalAnd(success);
 if (success) printLaph("ALL LOCAL TESTS PASSED!!");
 else         printLaph("Some local tests FAILED");
}

// ***********************************************


class TwoSpin
{

    unsigned int s1,s2;    // each value between 0 and 3 (say)
    
  public:
  
    TwoSpin(int in1, int in2) {assign(in1,in2);} // no default constructor by design !!
    
    TwoSpin(const unsigned int* buf) {assign(buf[0],buf[1]);}   

    TwoSpin(const TwoSpin& rhs) : s1(rhs.s1), s2(rhs.s2) {}
    
    TwoSpin& operator=(const TwoSpin& rhs)
     {s1=rhs.s1; s2=rhs.s2; return *this;}
    
    std::string output() const;
    
    bool operator<(const TwoSpin& rhs) const
    { return ((s1<rhs.s1) || ( (s1==rhs.s1)&&(s2<rhs.s2) ) ); }
    
    static int numints() {return 2;}
    
    size_t numbytes() const {return 2*sizeof(int);}
    
    void copyTo(unsigned int* buf) const { buf[0]=s1; buf[1]=s2;}
  
  private:
  
    void assign(int in1, int in2);

};

   // no default constructor is needed
   
inline void TwoSpin::assign(int in1, int in2)
{
 if ((in1<0)||(in1>3)||(in2<0)||(in2>3)){
    errorLaph("Invalid TwoSpin initialization");}
 s1=static_cast<unsigned int>(in1);
 s2=static_cast<unsigned int>(in2);
}

inline std::string  TwoSpin::output() const
{
 std::stringstream oss;
 oss << "("<< s1 <<", "<< s2 <<")";
 return oss.str();
}

void pprintf(ofstream& fout, bool global, const char* format, ...)
{
 va_list args;
 va_start(args, format);
 char buffer[2048];
 vsprintf (buffer, format, args);
 va_end(args);  
 if (global){
    printLaph(buffer);
    }
 else{
    fout << buffer;}
}

// **********************************************************************


void testIOMap(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestIOMap")==0) return;

 printLaph("\nNow starting IOMap tests");
 XMLHandler xmlr(xml_in,"TestIOMap");

 string filename,mode,endianness,whichtest,Rcksum,Wcksum;
 xmlread(xmlr,"FileName",filename,"TestIOMap");
 xmlread(xmlr,"Mode",mode,"TestIOMap");
 xmlread(xmlr,"Endianness",endianness,"TestIOMap");
 xmlread(xmlr,"WhichTest",whichtest,"TestIOMap");
 xmlread(xmlr,"WriteChecksumTest",Wcksum,"TestIOMap");
 xmlread(xmlr,"ReadChecksumTest",Rcksum,"TestIOMap");
 double tol=-0.1;
 xmlreadif(xmlr,"Tol",tol,"TestIOMap");
 
 char endian;
 if (endianness=="little") endian='L';
 else if (endianness=="big") endian='B';
 else endian='N';
 bool writetest,readtest;
 if (whichtest=="write"){ writetest=true; readtest=false;}
 else if (whichtest=="read"){ writetest=false; readtest=true;}
 else { writetest=true; readtest=true;}
 bool testreadchecksum;
 if (Rcksum=="yes") testreadchecksum=true; else testreadchecksum=false;
 bool testwritechecksum;
 if (Wcksum=="yes") testwritechecksum=true; else testwritechecksum=false;
 
 bool global;
 if (mode=="local") global=false; else global=true;
 
 //bool checkflag=true;

 printLaph("\n\nStarting IOMap tests");
 printLaph(make_strf("Lattice size is %d %d %d %d",LayoutInfo::getLattSizes()[0],
             LayoutInfo::getLattSizes()[1],LayoutInfo::getLattSizes()[2],
             LayoutInfo::getLattSizes()[3]));
 printLaph(make_strf("Number of MPI ranks = %d",LayoutInfo::getNumRanks()));
 printLaph(make_strf("Number of MPI ranks in each lattice direction %d %d %d %d",
             LayoutInfo::getCommNumPartitions()[0],LayoutInfo::getCommNumPartitions()[1],
             LayoutInfo::getCommNumPartitions()[2],LayoutInfo::getCommNumPartitions()[3]));
 printLaph(make_strf("Testing endianness = %c",endian));
 printLaph(make_strf("Main output file name stub = %s",filename));
 printLaph(make_strf("IOMap mode = %s",mode.c_str()));
 printLaph(make_strf("size of pos_type = %d",sizeof(IOHandlerFM::pos_type)));
 printLaph(make_strf("size of off_type = %d",sizeof(IOHandlerFM::off_type)));
 printLaph(make_strf("size of int = %d",sizeof(int)));
 printLaph(make_strf("size of long = %d",sizeof(long)));
 {ByteHandler mb;
 printLaph(make_strf("Is this machine bigendian? %d",mb.big_endian()));}
 int myrank=LayoutInfo::getMyRank();
 if (!global) filename+=int_to_string(myrank);
 string logfile("testiomap_"); logfile+=make_string(myrank)+".log";
 ofstream olog;
 if (!global) olog.open(logfile);

 IOMap<TwoSpin,vector<int> > iomap(global);
 iomap.setHighVerbosity();
 string header("This is the header string");

 if (writetest){

 printLaph("starting WRITE tests");

 iomap.openNew(filename,"TestIOMap ID String",header,false,endian,1,0,
               testwritechecksum,true);

 vector<int> g(3); g[0]=1; g[1]=2; g[2]=3;
 if (!global) for (int k=0;k<3;++k) g[k]+=10*myrank;
 iomap.put( TwoSpin(1,2), g);

 vector<int> h(2); h[0]=1; h[1]=2;
 if (!global) for (int k=0;k<2;++k) h[k]+=10*myrank;
 iomap.put( TwoSpin(3,2), h);

 vector<int> hh(2); hh[0]=5; hh[1]=8;
 if (!global) for (int k=0;k<2;++k) hh[k]+=10*myrank;
 iomap.put( TwoSpin(1,1), hh);

 vector<int> hhh(6); hhh[0]=55; hhh[1]=88; hhh[2]=66;
 hhh[3]=77; hhh[4]=44; hhh[5]=99;
 if (!global) for (int k=0;k<6;++k) hhh[k]+=10*myrank;
 iomap.put( TwoSpin(2,2), hhh);
 vector<int> hhhkeep(hhh);

 pprintf(olog,global,"\nHeader in file is %s\n",iomap.getHeader().c_str());
 pprintf(olog,global,"Header should be  This is the header string\n\n");

 hhh[0]=555; hhh[1]=888; hhh[2]=666;
 hhh[3]=777; hhh[4]=444; hhh[5]=999;
 if (!global) for (int k=0;k<6;++k) hhh[k]+=10*myrank;
 iomap.put( TwoSpin(2,0), hhh);
 
 iomap.flush();
 sleep(6);  // for testing flush and file corruption during abort

 vector<int> p(1); p[0]=8;
 if (!global) for (int k=0;k<1;++k) p[k]+=10*myrank;
 iomap.put( TwoSpin(1,1), p );

 vector<TwoSpin> keys;
 iomap.getKeys(keys);
 for (uint i=0;i<keys.size();i++){
    pprintf(olog,global,"key[%d] = %s\n",i,keys[i].output().c_str());}

 vector<int> val;
 iomap.get(TwoSpin(1,2),val);
 pprintf(olog,global,"val = ");
 for (uint j=0;j<val.size();j++) pprintf(olog,global," %d",val[j]);


 string results("\n results above should be ");
 for (uint k=0;k<g.size();++k) results+=" "+make_string(g[k]);
 results+="\n";
 pprintf(olog,global,results.c_str());

 hh[0]=5050; hh[1]=8080;
 if (!global) for (int k=0;k<2;++k) hh[k]+=10*myrank;
 iomap.put( TwoSpin(0,0), hh);

 iomap.put( TwoSpin(3,3), hhh);

 iomap.getKeys(keys);
 for (uint i=0;i<keys.size();i++){
    pprintf(olog,global,"key[%d] = %s\n",i,keys[i].output().c_str());}

 iomap.get(TwoSpin(2,2),val);
 pprintf(olog,global,"val = ");
 for (uint j=0;j<val.size();j++) pprintf(olog,global," %d",val[j]);
 results="\n results above should be ";
 for (uint k=0;k<hhhkeep.size();++k) results+=" "+make_string(hhhkeep[k]);
 results+="\n";
 pprintf(olog,global,results.c_str());

 iomap.close();
 }

 if (global){

 IOMap<TwoSpin,vector<LattField> > giom;
 giom.setHighVerbosity();
 string header("This is the header string for vector<LattField>");

 if (writetest){

 pprintf(olog,global,"starting LattField for ColorMatrix WRITE tests\n\n");
 
 giom.openNew(filename+"latt","TestIOMap ID String",header,false,endian,1,0,
              testwritechecksum,true);

 GaugeFieldAssigner GFA(3321294511,827120843,1942865517,4019325174);
 vector<LattField> GG;
 GFA.assign_gauge_field(GG,"GaugeField1");
 giom.put( TwoSpin(1,2), GG);

 GFA.reSeed(1196140740,809094426,2348838239,4264392720);
 GFA.assign_gauge_field(GG,"GaugeField2");
 giom.put( TwoSpin(3,2), GG);

 pprintf(olog,global,"\nHeader in file is %s\n",giom.getHeader().c_str());
 pprintf(olog,global," vector<LattField> of ColorMatrix write done\n\n");

 GG.clear();
 GG.resize(3);
 LatticeAssigner LFA(827120843,-7.0,7.0);
 GG[0].reset(FieldSiteType::Complex);
 LFA.assign_field(GG[0],"LatticeComplex");
 LFA.reSeed(4264392720);
 GG[1].reset(FieldSiteType::ColorVector);
 LFA.assign_field(GG[1],"LatticeColorVector");
 LFA.reSeed(2348838239);
 GG[2].reset(FieldSiteType::ColorSpinVector);
 LFA.assign_field(GG[2],"LatticeColorSpinVector");
 giom.put( TwoSpin(1,3), GG); 
 pprintf(olog,global," vector<LattField> of various site types write done\n\n");
 giom.close();
 }}


 if (readtest){
 
 IOMap<TwoSpin,vector<int> > ff(global);
 ff.setHighVerbosity();
 ff.openUpdate(filename,"TestIOMap ID String",header,endian,1,0,
               testreadchecksum,true);
 
 pprintf(olog,global,"%s\n",ff.getHeader().c_str());
 
 vector<int> qq(5); qq[0]=100; qq[1]=200; qq[2]=300; qq[3]=400; qq[4]=500;
 ff.put(TwoSpin(2,1),qq);
 
 qq[0]=1000; qq[1]=2000; qq[2]=3000; qq[3]=4000; qq[4]=5000;
 ff.put(TwoSpin(2,2),qq);
 
 
 vector<TwoSpin> newkeys;
 ff.getKeys(newkeys);
 for (uint i=0;i<newkeys.size();i++){
    pprintf(olog,global,"newkey[%d] = %s\n",i,newkeys[i].output().c_str());}
    
 ff.close();

 IOMap<TwoSpin,vector<int> > fff;
 fff.setHighVerbosity();
 fff.openReadOnly(filename,"TestIOMap ID String",header,true);
 pprintf(olog,global,"header = %s\n",header.c_str());

 vector<int> vv;
 fff.get(TwoSpin(2,1),vv);
 pprintf(olog,global,"vv = ");
 for (uint j=0;j<vv.size();j++) pprintf(olog,global," %d",vv[j]); 
 pprintf(olog,global,"\n");
 }


 if ((global)&&(readtest)){

 IOMap<TwoSpin,vector<LattField>> giom;
 giom.setHighVerbosity();

 printLaph("starting LatticeColorMatrix READ tests");
 
 string header;
 giom.openReadOnly(filename+"latt","TestIOMap ID String",header,
                   testwritechecksum);

 vector<LattField> inGG;
 giom.get( TwoSpin(1,2), inGG);

 GaugeFieldAssigner GFA(3321294511,827120843,1942865517,4019325174);
 GaugeFieldChecker GFC(GFA);
 bool check=GFC.check_gauge_field(inGG,"GaugeField1",true,tol);
 if (check) printLaph("GaugeField1 agrees");
 else printLaph("GaugeField1 MISMATCH");

 giom.get( TwoSpin(3,2), inGG);
 GFC.reSeed(1196140740,809094426,2348838239,4264392720);
 check=GFC.check_gauge_field(inGG,"GaugeField2",true,tol);
 if (check) printLaph("GaugeField2 agrees");
 else printLaph("GaugeField2 MISMATCH");

 printLaph(make_strf("\nHeader in file is %s",giom.getHeader()));
 printLaph(" vector<LattField> with ColorMatrix read done\n");

 giom.get( TwoSpin(1,3), inGG);
 LatticeAssigner LFA(827120843,-7.0,7.0);
 LatticeChecker LFC(LFA);
 check=LFC.check_field(inGG[0],"LatticeComplex",true,tol);
 LFA.reSeed(4264392720);
 check=check && LFC.check_field(inGG[1],"LatticeColorVector",true,tol);
 LFA.reSeed(2348838239);
 check=check && LFC.check_field(inGG[2],"LatticeColorSpinVector",true,tol);
 if (check) printLaph("vector<LattField> of various site types agrees");
 else printLaph("vector<LattField> of various site types MISMATCH");
 giom.close();
 } 

}


// *****************************************************************************



class SimpleHandler
{

 public:

    struct FileKey
    {
       int fkey;
       FileKey(int val) : fkey(val) {}
       FileKey(XMLHandler& xmlr){xmlread(xmlr,"FileKey",fkey,"");}
       FileKey(const FileKey& rhs) : fkey(rhs.fkey) {}
       FileKey& operator=(const FileKey& rhs) {fkey=rhs.fkey; return *this;}
       ~FileKey(){}
       void output(XMLHandler& xmlw) const {xmlw.set_root("FileKey",make_string(fkey));}
       bool operator<(const FileKey& rhs) const {return (fkey<rhs.fkey);}
       bool operator==(const FileKey& rhs) const {return (fkey==rhs.fkey);}
    };

    struct RecordKey
    {
       int rkey;
       RecordKey() : rkey(0) {}
       RecordKey(int val) {rkey=val;}
       RecordKey(const RecordKey& in) {rkey=in.rkey;}
       ~RecordKey() {}
       void output(XMLHandler& xmlw) const {xmlw.set_root("RecordKey",make_string(rkey));}
       bool operator<(const RecordKey& rhs) const {return (rkey<rhs.rkey);}
       bool operator==(const RecordKey& rhs) const {return (rkey==rhs.rkey);}
       explicit RecordKey(const unsigned int* buf){ rkey=*buf;}
       void copyTo(unsigned int* buf) const {*buf=rkey;}
       static int numints() {return 1;}
       size_t numbytes() const {return sizeof(int);}
       
    };


    struct DataType
    {
       int result;
       DataType() {}
       DataType(int in) : result(in) {}
       DataType(const DataType& in) : result(in.result) {}
       DataType& operator=(const DataType& in)
        {result=in.result; return *this;}
       ~DataType() {}
      // void output(TextFileWriter& tout) const         // HERE
      // {tout << result;}
       size_t numbytes() const {return sizeof(int);}
    };

 private:

    GaugeConfigurationInfo m_gauge;
    FileListInfo m_files;
    bool write_mode;
    
 //   typedef DataPutHandlerMF<SimpleHandler,FileKey,RecordKey,DataType> DPutType;
 //   typedef DataGetHandlerMF<SimpleHandler,FileKey,RecordKey,DataType> DGetType;
    
 //   typedef DataPutHandlerMFNOM<SimpleHandler,FileKey,RecordKey,DataType> DPutType;
 //   typedef DataGetHandlerMFNOM<SimpleHandler,FileKey,RecordKey,DataType> DGetType;

    typedef DataPutHandlerMFO<SimpleHandler,FileKey,RecordKey,DataType> DPutType;
    typedef DataGetHandlerMFO<SimpleHandler,FileKey,RecordKey,DataType> DGetType;

    DPutType *DHinsert;
    DGetType *DHread;  

 public:

    SimpleHandler(const GaugeConfigurationInfo& ingauge,
                  const FileListInfo& files, bool writemode=false);

    void setWriteMode(bool over_write=false);

    void setReadMode();

    ~SimpleHandler() {delete DHinsert; delete DHread;}

    bool checkHeader(XMLHandler& xmlr, int suffix);

    void writeHeader(XMLHandler& xmlout, const FileKey& fkey, int suffix);

    void putData(int val1, int val2, int data);

    const DataType& getData(int val1, int val2);

    bool queryData(int val1, int val2);

    bool queryFile(int val1);

    void removeData(int val1, int val2);

    void getFileMap(XMLHandler& xmlout) const;

    void getKeys(XMLHandler& xmlout);

    void outputKeys(XMLHandler& xmlout);

    void clearData() {if (DHread!=0) DHread->clearData();}

    void merge(const FileListInfo& mergefiles);
    
    const FileListInfo& getFileListInfo() const;




    void open(int val1);

    void putData(int val2, int data);  // insert into current file

    void flush();   // flush current file

    void close();   // close current file

    void flush(int val1);

    void close(int val1);
    
    void flushAll();

    void closeAll();


    bool queryData(int val2);  // query in current open file



};


SimpleHandler::SimpleHandler(const GaugeConfigurationInfo& ingauge,
                             const FileListInfo& files, bool writemode) 
         : m_gauge(ingauge), m_files(files), write_mode(writemode), 
           DHinsert(0), DHread(0)
{ printLaph("Starting SimpleHandler constructor");
 if (writemode){ printLaph("write mode firing up DHinsert");
    DHinsert=new DPutType(*this,m_files,"FileID","SimpleHeader"); printLaph("done DHinsert");}
 else{ printLaph("read mode firing up DHread");
    DHread=new DGetType(*this,m_files,"FileID","SimpleHeader");printLaph("done DHread");}
}


bool SimpleHandler::checkHeader(XMLHandler& xmlr, int suffix)
{
printLaph("starting checkHeader");
 GaugeConfigurationInfo gauge_info(xmlr); printLaph("read gaugeinfo");
 try{printLaph("trying checkequal gaugeinfo");
     m_gauge.checkEqual(gauge_info);}
 catch(const std::exception& xp){ printLaph("result is false"); return false;}
 printLaph("result is true");
 return true;
}

void SimpleHandler::writeHeader(XMLHandler& xmlout, 
                               const SimpleHandler::FileKey& fkey,
                               int suffix)
{
printLaph("about to write header");
 xmlout.set_root("SimpleHeader"); printLaph(make_strf("%s",xmlout.str()));
 XMLHandler xmlg,xmlk;
 m_gauge.output(xmlg);printLaph(make_strf("%s",xmlg.str()));
 fkey.output(xmlk);printLaph(make_strf("%s",xmlk.str()));
 xmlout.put_child(xmlg);
 xmlout.put_child(xmlk);
              // HERE
      printLaph("done writing header");
} 

void SimpleHandler::setWriteMode(bool over_write)
{
 if (!write_mode){
    delete DHread; DHread=0;
    write_mode=true;
    DHinsert=new DPutType(*this,m_files,"FileID","SimpleHeader");}
 if (over_write) DHinsert->setOverWrite();
 else DHinsert->setNoOverWrite();
}

void SimpleHandler::setReadMode()
{
 if (!write_mode) return;
 delete DHinsert; DHinsert=0;
 write_mode=false;
 DHread=new DGetType(*this,m_files,"FileID","SimpleHeader");
}

void SimpleHandler::putData(int val1, int val2, int data)
{
 if (!write_mode){
    errorLaph("attempt to insert data before setWriteMode called");}
 DHinsert->putData(FileKey(val1),RecordKey(val2),DataType(data));
}

bool SimpleHandler::queryData(int val1, int val2)
{
 if (write_mode) return DHinsert->queryData(FileKey(val1),RecordKey(val2));
 else return DHread->queryData(FileKey(val1),RecordKey(val2));
}

const FileListInfo& SimpleHandler::getFileListInfo() const
{
 if (write_mode) return DHinsert->getFileListInfo();
 else return DHread->getFileListInfo();
}

bool SimpleHandler::queryFile(int val1)
{
 if (write_mode) return DHinsert->queryFile(FileKey(val1));
 else return DHread->queryFile(FileKey(val1));
}

const SimpleHandler::DataType& SimpleHandler::getData(int val1, int val2)
{
 if (write_mode){
    errorLaph("attempt to get data before setReadMode called");}
 return DHread->getData(FileKey(val1),RecordKey(val2));
}

void SimpleHandler::removeData(int val1, int val2)
{
 if (write_mode){
    errorLaph("attempt to get data before setReadMode called");}
 DHread->removeData(FileKey(val1),RecordKey(val2));
}

void SimpleHandler::getFileMap(XMLHandler& xmlout) const
{
 if (!write_mode) DHread->getFileMap(xmlout);
}

void SimpleHandler::getKeys(XMLHandler& xmlout)
{
 if (!write_mode){ /*
    push(xmlout,"SimpleHandlerKeys");
    multi1d<int> key(2);
    list<pair<FileKey,list<RecordKey> > > keys(DHread->getKeys());
    for (list<pair<FileKey,list<RecordKey> > >::const_iterator 
          it=keys.begin();it!=keys.end();it++){
       key[0]=it->first.fkey;
       for (list<RecordKey>::const_iterator 
               jt=it->second.begin();jt!=it->second.end();jt++){
          key[1]=jt->rkey;
          write(xmlout,"Key",key);}}
    pop(xmlout); */}   // HERE
}

void SimpleHandler::outputKeys(XMLHandler& xmlout)
{
  if (!write_mode) DHread->outputKeys(xmlout);
}

void SimpleHandler::merge(const FileListInfo& mergefiles)
{
 if (!write_mode){
    errorLaph("attempt to merge before setWriteMode called");}
// DHinsert->merge(mergefiles,"SimpleHeader");
}


void SimpleHandler::open(int val1)
{
  if (write_mode) DHinsert->open(FileKey(val1));
}

void SimpleHandler::putData(int val2, int data)  // insert into current file
{
  if (write_mode) DHinsert->putData(RecordKey(val2),DataType(data));
}

void SimpleHandler::flush()   // flush current file
{
  if (write_mode) DHinsert->flush();
}

void SimpleHandler::close()   // close current file
{
  if (write_mode) DHinsert->close();
}

void SimpleHandler::flush(int val1)
{
  if (write_mode) DHinsert->flush(FileKey(val1));
}

void SimpleHandler::close(int val1)
{
  if (write_mode) DHinsert->close(FileKey(val1));
}
    
void SimpleHandler::flushAll()
{
  if (write_mode) DHinsert->flushAll();
}

void SimpleHandler::closeAll()
{
  if (write_mode) DHinsert->closeAll();
}


bool SimpleHandler::queryData(int val2)  // query in current open file
{
  if (write_mode) return DHinsert->queryData(RecordKey(val2));
  return false;
}


/*                 // HERE
void write(BinaryWriter& xmlout, const SimpleHandler::RecordKey& rkey)
{
 rkey.binaryWrite(xmlout);
}

void readPrimaryNode(BinaryReader& xmlin, SimpleHandler::RecordKey& rkey)
{
 rkey.binaryReadPrimaryNode(xmlin);
}

void read(BinaryReader& xmlin, SimpleHandler::RecordKey& rkey)
{
 rkey.binaryRead(xmlin);
}

void write(BinaryWriter& xmlout, const SimpleHandler::DataType& data)
{
 write(xmlout,data.result);
}

void read(BinaryReader& xmlin, SimpleHandler::DataType& data)
{
 read(xmlin,data.result);
}

void output(const SimpleHandler::DataType& res, TextFileWriter& tout)
{
 res.output(tout);
}
*/
void write(IOHandlerFM& ioh, const SimpleHandler::DataType& data)
{
 ioh.write(data.result);
}

void read(IOHandlerFM& ioh, SimpleHandler::DataType& data)
{
 ioh.read(data.result);
}


// ***********************************************************************************


class SimpHandler
{

 public:

    struct DataType
    {
       int result;
       DataType() {}
       DataType(int in) : result(in) {}
       DataType(const DataType& in) : result(in.result) {}
       DataType& operator=(const DataType& in)
        {result=in.result; return *this;}
       ~DataType() {}
       //void output(TextFileWriter& tout) const
       //{tout << result;}
       size_t numbytes() const {return sizeof(int);}
    };

    struct RecordKey
    {
       int rkey;
       RecordKey() : rkey(0) {}
       RecordKey(int val) {rkey=val;}
       RecordKey(const RecordKey& in) {rkey=in.rkey;}
       ~RecordKey() {}
       void output(XMLHandler& xmlw) const {/*write(xmlw,"RecordKey",rkey);*/}   // HERE
       bool operator<(const RecordKey& rhs) const {return (rkey<rhs.rkey);}
       bool operator==(const RecordKey& rhs) const {return (rkey==rhs.rkey);}
       //void binaryReadPrimaryNode(BinaryReader& bin)
       // { bin.readArrayPrimaryNode((char*)&rkey, sizeof(rkey), 1);}
       //void binaryRead(BinaryReader& bin) {read(bin,rkey);}
       //void binaryWrite(BinaryWriter& bout) const {write(bout,rkey);}
       explicit RecordKey(const unsigned int* buf){ rkey=*buf;}
       void copyTo(unsigned int* buf) const {*buf=rkey;}
       static int numints() {return 1;}
       size_t numbytes() const {return sizeof(int);}
    };


 private:

    GaugeConfigurationInfo m_gauge;
    string m_file;
    int m_mode;
    
   // typedef DataPutHandlerSF<SimpHandler,RecordKey,DataType> DPutType;
   // typedef DataGetHandlerSF<SimpHandler,RecordKey,DataType> DGetType;

   // typedef DataPutHandlerSFNOM<SimpHandler,RecordKey,DataType> DPutType;
   // typedef DataGetHandlerSFNOM<SimpHandler,RecordKey,DataType> DGetType;

    typedef DataPutHandlerSFO<SimpHandler,RecordKey,DataType> DPutType;
    typedef DataGetHandlerSFO<SimpHandler,RecordKey,DataType> DGetType;

    DPutType *DHinsert;
    DGetType *DHread;  

 public:

    SimpHandler(const GaugeConfigurationInfo& ingauge,
                const string& file, int mode=0);  // 0=read,1=write,2=overwrite

    void setMode(int mode);

    ~SimpHandler() {delete DHinsert; delete DHread;}

    void putData(int val, int data);

    const DataType& getData(int val);

    bool queryData(int val);
    
    string getFileName() const;

    void removeData(int val);

    void outputKeys(XMLHandler& xmlout);

    void getKeys(XMLHandler& xmlout);

    void clearData() {if (DHread!=0) DHread->clearData();}

    bool checkHeader(XMLHandler& xmlr);

    void writeHeader(XMLHandler& xmlout);

 private:
 

    friend class DataPutHandlerSF<SimpHandler,RecordKey,DataType>;
    friend class DataGetHandlerSF<SimpHandler,RecordKey,DataType>;  
    //friend void write(BinaryWriter& xmlout, const RecordKey& rkey);
    //friend void readPrimaryNode(BinaryReader& xmlin, RecordKey& rkey);
    //friend void read(BinaryReader& xmlin, RecordKey& rkey);

};


SimpHandler::SimpHandler(const GaugeConfigurationInfo& ingauge,
                         const string& file, int mode) 
         : m_gauge(ingauge), m_file(file), DHinsert(0), DHread(0)
{
 m_mode=-1;
 setMode(mode);
}


bool SimpHandler::checkHeader(XMLHandler& xmlr)
{
 GaugeConfigurationInfo gauge_info(xmlr);
 try{
     m_gauge.checkEqual(gauge_info);}
 catch(const std::exception& xp){ return false;}
 return true;
}

void SimpHandler::writeHeader(XMLHandler& xmlout)
{/*
 push(xmlout,"SimpHeader");
 m_gauge.output(xmlout);
 pop(xmlout);*/            // HERE
}

void SimpHandler::setMode(int mode)
{
 if ((mode<0)||(mode>2)){
    errorLaph("invalid mode in SimpHandler::setMode");}
 if (mode==m_mode) return;
 if (mode==0){
    m_mode=0;
    delete DHinsert; DHinsert=0;
    DHread=new DGetType(*this,m_file,"FileID");}
 else if (mode==1){
    m_mode=1;
    delete DHread; DHread=0;
    DHinsert=new DPutType(*this,m_file,"FileID",false);}
 else{
    m_mode=2; 
    delete DHread; DHread=0;
    DHinsert=new DPutType(*this,m_file,"FileID",true);}
}


void SimpHandler::putData(int val, int data)
{
 if (m_mode==0){
    errorLaph("attempt to insert data before setWriteMode called");}
 DHinsert->putData(RecordKey(val),DataType(data));
}

bool SimpHandler::queryData(int val)
{
 if (m_mode==0) return DHread->queryData(RecordKey(val));
 else return DHinsert->queryData(RecordKey(val));
}

string SimpHandler::getFileName() const
{
 if (m_mode==0) return DHread->getFileName();
 else return DHinsert->getFileName();
}

const SimpHandler::DataType& SimpHandler::getData(int val)
{
 if (m_mode!=0){
    errorLaph("attempt to get data before setReadMode called");}
 return DHread->getData(RecordKey(val));
}

void SimpHandler::removeData(int val)
{
 if (m_mode!=0){
    errorLaph("attempt to get data before setReadMode called");}
 DHread->removeData(RecordKey(val));
}

void SimpHandler::getKeys(XMLHandler& xmlout)
{/*
 if (m_mode==0){
    push(xmlout,"SimpHandlerKeys");
    set<RecordKey> keys(DHread->getKeys());
    for (set<RecordKey>::const_iterator 
            jt=keys.begin();jt!=keys.end();jt++)
       write(xmlout,"Key",jt->rkey);
    pop(xmlout);} */    // HERE
}

void SimpHandler::outputKeys(XMLHandler& xmlout)
{
  if (m_mode==0) DHread->outputKeys(xmlout);
}

/*
void write(BinaryWriter& xmlout, const SimpHandler::RecordKey& rkey)
{
 rkey.binaryWrite(xmlout);
}

void readPrimaryNode(BinaryReader& xmlin, SimpHandler::RecordKey& rkey)
{
 rkey.binaryReadPrimaryNode(xmlin);
}

void read(BinaryReader& xmlin, SimpHandler::RecordKey& rkey)
{
 rkey.binaryRead(xmlin);
}

void write(BinaryWriter& xmlout, const SimpHandler::DataType& data)
{
 write(xmlout,data.result);
}

void read(BinaryReader& xmlin, SimpHandler::DataType& data)
{
 read(xmlin,data.result);
}
*/
void write(IOHandlerFM& ioh, const SimpHandler::DataType& data)
{
 ioh.write(data.result);
}

void read(IOHandlerFM& ioh, SimpHandler::DataType& data)
{
 ioh.read(data.result);
}



// ***************************************************************************

void testDataIOHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestDataIOHandler")==0)
 return;

 XMLHandler xmlr(xml_in,"TestDataIOHandler");
 printLaph("\nIn TestDataIOHandler\n");

 GaugeConfigurationInfo gauge(xmlr);
 LaphNoiseInfo noise1(xmlr);

 if (xml_tag_count(xmlr,"DBTest")==1){
 FileListInfo files(xmlr);

printLaph("Start to assign SH");
 {SimpleHandler SH(gauge,files);  printLaph("Done constructor");
 SH.setWriteMode(false);
 SH.putData(3,2,5);   printLaph("insert (3,2)");
 SH.putData(3,3,-4);  printLaph("insert (3,3)");
 SH.putData(0,0,4);   printLaph("insert (0,0)");
 SH.putData(0,1,34);  printLaph("insert (0,1)");
 SH.putData(0,2,24);  printLaph("insert (0,2)");
 SH.putData(0,3,14);  printLaph("insert (0,3)");
 SH.putData(3,4,1);   printLaph("insert (3,4)");
 SH.putData(2,5,13);  printLaph("insert (2,5)");
// SH.putData(2,5,89);   printLaph("insert (2,5)");  // will fail due to no overwrite
 SH.putData(1,1,-65); printLaph("insert (1,1)"); }

 printLaph("first set of tests done");

printLaph("Start to assign SH 2nd time");
 {SimpleHandler SH(gauge,files); printLaph("Done constructor");
 SH.setWriteMode(true);
 SH.putData(5,12,321); printLaph("insert (5,12)");
 SH.putData(5,4,54);   printLaph("insert (5,4)");
 SH.putData(5,7,33);   printLaph("insert (5,7)");
 SH.putData(5,4,22);   printLaph("insert (5,4)");
 SH.putData(3,4,11);   printLaph("insert (3,4)");
 SH.putData(6,7,42);   printLaph("insert (6,7)");
 SH.putData(2,5,99);   printLaph("insert (2,5)");
 SH.putData(8,-3,-66); printLaph("insert (8,-3)");

 printLaph(make_strf("query 5,12 should be yes: %d",SH.queryData(5,12)));
 printLaph(make_strf("query 8,18 should be no: %d",SH.queryData(8,18)));
 printLaph(make_strf("query 3 file should be yes: %d",SH.queryFile(3)));
 printLaph(make_strf("query 11 file should be no: %d",SH.queryFile(11)));
 }

 printLaph("Starting reader\n");

 SimpleHandler SH(gauge,files);
 SH.setReadMode();

 const SimpleHandler::DataType& result1=SH.getData(3,3);
 printLaph(make_strf("result(3,3) = %d",result1.result));
/*
 XmlBufferWriter xmlout1;
 SH.getFileMap(xmlout1);
 QDPIO::cout << "SH file map:"<<endl<<xmlout1.str()<<endl;

// XmlBufferWriter xmlout2;
// SH.getKeys(xmlout2);
// QDPIO::cout << "SH keys:"<<endl<<xmlout2.str()<<endl;

 XmlBufferWriter xmlout3;
 SH.outputKeys(xmlout3);
 QDPIO::cout << "SH keys:"<<endl<<xmlout3.str()<<endl;

// SH.clearData(); cout << "cleared"<<endl;

 const SimpleHandler::DataType& result2=SH.getData(2,5);
 QDPIO::cout << "result(2,5) = "<<result2.result<<endl;

 const SimpleHandler::DataType& result3=SH.getData(0,0);
 QDPIO::cout << "result(0,0) = "<<result3.result<<endl;

 const SimpleHandler::DataType& result5=SH.getData(5,4);
 QDPIO::cout << "result(5,4) = "<<result5.result<<endl;

 const SimpleHandler::DataType& result6=SH.getData(5,12);
 QDPIO::cout << "result(5,12) = "<<result6.result<<endl;

 QDPIO::cout << "query 0 file should be yes: "<<SH.queryFile(0)<<endl;
 QDPIO::cout << "query 19 file should be no: "<<SH.queryFile(19)<<endl;
 QDPIO::cout << "query data should be yes: "<<SH.queryData(5,12)<<endl;
 QDPIO::cout << "query data should be no: "<<SH.queryData(5,11)<<endl;
 QDPIO::cout << "query data should be yes: "<<SH.queryData(5,4)<<endl;

 const SimpleHandler::DataType& result4=SH.getData(6,7);
 QDPIO::cout << "result(6,7) = "<<result4.result<<" should be 42"<<endl;

 FileListInfo flist1("DummySimpA",0,10);
 {SimpleHandler SH1(gauge,flist1);
 SH1.setWriteMode(true);
 SH1.putData(2,10,20);
 SH1.putData(3,11,33);  
 SH1.putData(4,12,48);  
 SH1.putData(2,13,26);  
 SH1.putData(3,14,42);  
 SH1.putData(4,15,60);  
 SH1.putData(6,16,96);  
 SH1.putData(6,17,102);  
 }
 
 FileListInfo flist2("DummySimpA",11,20);
 {SimpleHandler SH1(gauge,flist2);
 SH1.setWriteMode(true);
 SH1.putData(2,18,36);
 SH1.putData(3,19,57);  
 SH1.putData(4,20,80);  
 SH1.putData(2,21,42);  
 SH1.putData(3,22,66);  
 SH1.putData(4,23,96);  
 SH1.putData(6,24,96);  
 SH1.putData(6,25,150);  
 }

 FileListInfo flist3("DummySimpA",21,30);
 {SimpleHandler SH1(gauge,flist3);
 SH1.setWriteMode(true);
 SH1.putData(2,26,52);
 SH1.putData(3,27,81);  
 SH1.putData(4,28,112);  
 SH1.putData(2,29,58);  
 SH1.putData(3,30,90);  
 SH1.putData(4,31,124);  
 SH1.putData(6,32,192);  
 SH1.putData(6,33,198);  
 }

 SimpleHandler SHT(gauge,flist3);
 {const SimpleHandler::DataType& resultV=SHT.getData(6,32);
 QDPIO::cout << "result(6,32) = "<<resultV.result<<" should be 192"<<endl;}
 {const SimpleHandler::DataType& resultV=SHT.getData(3,27);
 QDPIO::cout << "result(3,27) = "<<resultV.result<<" should be 81"<<endl;}
*/
 }
/*
 if (xml_tag_count(xmlr,"DBMergeTest")==1){

    XMLHandler xmlrr(xmlr,"./descendant-or-self::HandlerList");
    FileListInfo handlerfiles(xmlrr);
    XMLHandler xmlrm(xmlr,"./descendant-or-self::MergeList");
    FileListInfo mergefiles(xmlrm);

    {SimpleHandler SH(gauge,handlerfiles);
    SH.setWriteMode();
    SH.merge(mergefiles);}

    XmlBufferWriter xmlout3;
    SimpleHandler SH(gauge,handlerfiles);
    SH.outputKeys(xmlout3);
    QDPIO::cout << "SH keys:"<<endl<<xmlout3.str()<<endl;

    {const SimpleHandler::DataType& resultV=SH.getData(6,16);
    QDPIO::cout << "result(6,16) = "<<resultV.result<<" should be 96"<<endl;}
    {const SimpleHandler::DataType& resultV=SH.getData(2,13);
    QDPIO::cout << "result(2,13) = "<<resultV.result<<" should be 26"<<endl;}
    {const SimpleHandler::DataType& resultV=SH.getData(3,11);
    QDPIO::cout << "result(3,11) = "<<resultV.result<<" should be 33"<<endl;}
    {const SimpleHandler::DataType& resultV=SH.getData(4,20);
    QDPIO::cout << "result(4,20) = "<<resultV.result<<" should be 80"<<endl;}
    {const SimpleHandler::DataType& resultV=SH.getData(3,27);
    QDPIO::cout << "result(3,27) = "<<resultV.result<<" should be 81"<<endl;}
 }

// *******************************************************************

* /

/ *
 if (xml_tag_count(xmlr,"DMTest")==1){
 FileListInfo files(xmlr);

// ***********************************************

 QDPIO::cout << "Testing DM IO Handler"<<endl<<endl;

 {SimpleDMHandler SH(gauge,files);
 SH.setWriteMode();
 SH.putData(3,2,5);   QDPIO::cout << "insert (3,2)"<<endl;
 SH.putData(3,3,-4);  QDPIO::cout << "insert (3,3)"<<endl;
 SH.putData(3,4,1);   QDPIO::cout << "insert (3,4)"<<endl;
 SH.putData(0,0,4);   QDPIO::cout << "insert (0,0)"<<endl;
 SH.putData(0,1,34);  QDPIO::cout << "insert (0,1)"<<endl;
 SH.putData(0,2,24);  QDPIO::cout << "insert (0,2)"<<endl;
 SH.putData(0,3,14);  QDPIO::cout << "insert (0,3)"<<endl;
 SH.putData(2,5,13);  QDPIO::cout << "insert (2,5)"<<endl;
 SH.putData(1,1,-65); QDPIO::cout << "insert (1,1)"<<endl; }

 {SimpleDMHandler SH(gauge,files);
 SH.setWriteMode(true);
 SH.putData(5,12,321); QDPIO::cout << "insert (5,12)"<<endl;
 SH.putData(5,4,54);   QDPIO::cout << "insert (5,4)"<<endl;
 SH.putData(5,7,33);   QDPIO::cout << "insert (5,7)"<<endl;
 SH.putData(5,4,22);   QDPIO::cout << "insert (5,4)"<<endl;
 SH.putData(8,-3,-66); QDPIO::cout << "insert (8,-3)"<<endl;

 QDPIO::cout << "query 8,-3 should be yes: "<<SH.queryDataInCurrentFile(-3)<<endl;
 QDPIO::cout << "query 8,8 should be no: "<<SH.queryDataInCurrentFile(8)<<endl;

 QDPIO::cout << "query 3 file should be yes: "<<SH.queryFile(3)<<endl;
 QDPIO::cout << "query 11 file should be no: "<<SH.queryFile(11)<<endl;
 }

 QDPIO::cout << "Starting reader"<<endl;

 SimpleDMHandler SHDM(gauge,files);
 SHDM.setReadMode();

 const SimpleDMHandler::DataType& DMresult1=SHDM.getData(3,3);
 QDPIO::cout << "DMresult(3,3) = "<<DMresult1.result<<endl;

 XmlBufferWriter DMxmlout1;
 SHDM.getFileMap(DMxmlout1);
 QDPIO::cout << "SHDM file map:"<<endl<<DMxmlout1.str()<<endl;

 XmlBufferWriter DMxmlout2;
 SHDM.getKeys(DMxmlout2);
 QDPIO::cout << "SHDM keys:"<<endl<<DMxmlout2.str()<<endl;

 XmlBufferWriter DMxmlout3;
 SHDM.outputKeys(DMxmlout3);
 QDPIO::cout << "SHDM keys:"<<endl<<DMxmlout3.str()<<endl;

// SHDM.clearData(); cout << "cleared"<<endl;

 const SimpleDMHandler::DataType& DMresult2=SHDM.getData(2,5);
 QDPIO::cout << "DMresult(2,5) = "<<DMresult2.result<<endl;

 const SimpleDMHandler::DataType& DMresult3=SHDM.getData(0,0);
 QDPIO::cout << "DMresult(0,0) = "<<DMresult3.result<<endl;

 const SimpleDMHandler::DataType& DMresult5=SHDM.getData(5,4);
 QDPIO::cout << "DMresult(5,4) = "<<DMresult5.result<<endl;

 SHDM.removeData(5,12);
 QDPIO::cout << "query 5,12 should be yes: "<<SHDM.queryData(5,12)<<endl;
 QDPIO::cout << "query 8,18 should be no: "<<SHDM.queryData(8,18)<<endl;

 QDPIO::cout << "query 0 file should be yes: "<<SHDM.queryFile(0)<<endl;
 QDPIO::cout << "query 19 file should be no: "<<SHDM.queryFile(19)<<endl;

// const SimpleDMHandler::DataType& DMresult4=SHDM.getData(6,7);
// QDPIO::cout << "DMresult(6,7) = "<<DMresult4.result<<endl;

 } */
/*
 FileListDBPutHandler<LaphNoiseInfo,QuarkLineEndInfo,multi2d<Complex> > 
        DBin(files,"DummyHeader",&IOHandlerDummyChecker,

 FileListDBGetHandler<LaphNoiseInfo,int,int,double> DB(files,"Header",
                                              &IOHandlerDummyChecker);

 DB.queryData(noise1);
 DB.queryData(noise1,3,4);

 double res;
 res=DB.getData(noise1,3,4);

 XmlBufferWriter DMxmlout;
 DB.getFileMap(DMxmlout);
*/

/*
 if (xml_tag_count(xmlr,"FileTest")==1){

 string filename;
 read(xmlr,"FileName",filename);

 {SimpHandler SH(gauge,filename,1);
 SH.putData(4,5);   QDPIO::cout << "insert (4)"<<endl;
 SH.putData(8,-4);  QDPIO::cout << "insert (8)"<<endl;
 SH.putData(7,4);   QDPIO::cout << "insert (7)"<<endl;
 SH.putData(5,34);  QDPIO::cout << "insert (5)"<<endl;
 SH.putData(0,24);  QDPIO::cout << "insert (0)"<<endl;
 SH.putData(1,14);  QDPIO::cout << "insert (1)"<<endl;
// SH.putData(8,-65); QDPIO::cout << "insert (8) again"<<endl; 
 }

 {SimpHandler SH(gauge,filename,2);
 QDPIO::cout << "query 5 should be yes: "<<SH.queryData(5)<<endl;
 SH.putData(5,321);  QDPIO::cout << "insert (5)"<<endl;
 SH.putData(4,54);   QDPIO::cout << "insert (4)"<<endl;
 SH.putData(7,33);   QDPIO::cout << "insert (7)"<<endl;
 SH.putData(6,22);   QDPIO::cout << "insert (6)"<<endl;

 QDPIO::cout << "query 5 should be yes: "<<SH.queryData(5)<<endl;
 QDPIO::cout << "query 18 should be no: "<<SH.queryData(18)<<endl;

 } 

 QDPIO::cout << "Does file "<<filename<<" exist? "<<fileExists(filename)<<endl;
 QDPIO::cout << "Is it empty file name?"<<emptyFileName(filename)<<endl;
 QDPIO::cout << "Here is empty file name: "<<emptyFileName("    ")<<endl;
 QDPIO::cout << "Here is empty file name: "<<emptyFileName(" NOM_   ")<<endl;
 

 QDPIO::cout << "Starting reader"<<endl;

 SimpHandler SH(gauge,filename);

 const SimpHandler::DataType& res1=SH.getData(4);
 QDPIO::cout << "result(4) = "<<res1.result<<endl;
 QDPIO::cout << "filename = "<<SH.getFileName()<<endl;

 QDPIO::cout << "query 5 should be yes: "<<SH.queryData(5)<<endl;

// XmlBufferWriter xmloutB;
// SH.getKeys(xmloutB);
// QDPIO::cout << "SH keys:"<<endl<<xmloutB.str()<<endl;

 XmlBufferWriter xmloutC;
 SH.outputKeys(xmloutC);
 QDPIO::cout << "SH keys:"<<endl<<xmloutC.str()<<endl;

 SH.clearData(); cout << "cleared"<<endl;

 const SimpHandler::DataType& res2=SH.getData(1);
 QDPIO::cout << "result(1) = "<<res2.result<<endl;

 const SimpHandler::DataType& res3=SH.getData(8);
 QDPIO::cout << "result(8) = "<<res3.result<<endl;

 const SimpHandler::DataType& res5=SH.getData(5);
 QDPIO::cout << "result(5) = "<<res5.result<<endl;

 SH.removeData(5);

// const SimpHandler::DataType& res4=SH.getData(6);
// QDPIO::cout << "result(6) = "<<res4.result<<endl; 

 }
 
 if (xml_tag_count(xmlr,"OutputTest")==1){

 string filename;
 read(xmlr,"FileName",filename);
 IOMap<SimpleHandler::RecordKey,SimpleHandler::DataType> iomap;
 iomap.openReadOnly(filename,"logfile",true);

 }



 QDPIO::cout << " ***********************"<<endl<<endl;  */
/* IOHandler ioh;
 ioh.openReadOnly("tester_smeared_gauge_file","Laph--SmearedGaugeField");
 unsigned long mapstart;
 ioh.read(mapstart);
 QDPIO::cout << "mapstart = "<<mapstart<<endl;
 char checksum;
 ioh.read(checksum);
 QDPIO::cout << "checksum = "<<checksum<<endl;
 string header;
 ioh.read(header);
 QDPIO::cout << "header string = "<<header<<endl;

 QDPIO::cout << " ***********************"<<endl<<endl; */

}


// ************************************************************************

















// **********************************************************************************************
}
#endif
