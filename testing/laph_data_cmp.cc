#include <vector>
#include <iostream>
#include "xml_handler.h"
#include "io_map.h"
#include "layout_info.h"
#include "latt_field.h"

using namespace std;
using namespace LaphEnv;

// *   <QudaLaphCompare>
// *       <LatticeLayoutInfo>                                          *
// *         <XYZTSizes>24 24 24 96</XYZTSizes>                         *
// *       </LatticeLayoutInfo>                                         *
// *          ....comparison info
// *   </QudaLaphCompare>

// *       <Type>SmearedGaugeField</Type>
// *       <FileName1>...</FileName1>
// *       <FileName2>...</FileName2>

// *       <Type>LaphEigenvectors</Type>
// *       <NumEigvecs>48</NumEigvecs>
// *       <FileStub1>...</FileStub1>
// *       <FileStub2>...</FileStub2>

#ifdef ARCH_SERIAL

namespace LaphEnv {

bool isPrimaryRank()
{
 return true;
}

void laph_abort()
{
 exit(1);
}

}

typedef vector<LattField> SmearedGaugeDataType;
typedef LattField LaphEigvecType;

template <typename K, typename D>
void do_vec_compare(const K* p1, const D* p2, int nreal, double eps)
{
 cout.precision(12);
 const K* pp1=p1;
 const D* pp2=p2;
 for (int k=0;k<nreal;++k,++pp1,++pp2){
    if (std::abs((*pp1)-(*pp2))>eps){
       cout << "MISMATCH on element number "<<k<<": "<<*pp1<<" "<<*pp2<<endl;}}
    //else{
    //   cout << "AGREE: "<<*pp1<<" "<<*pp2<<endl;}}
}

/*
template <typename K, typename D>
void checker(IOMap<K,D>& iom1, IOMap<K,D>& iom2)
{
 string header1,header2;
 header1=iom1.getHeader();
 header2=iom2.getHeader();
 cout << "Header 1: "<<header1<<endl;
 cout << "Header 2: "<<header2<<endl;
 vector<UIntKey> keys1,keys2;
 iom1.getKeys(keys1);
 iom2.getKeys(keys2);
 cout << "Keys in file 1:"<<endl;
 for (int k=0;k<int(keys1.size());++k){
    cout << "Key "<<k<<": "<<keys1[k].output()<<endl;}
 cout << "Keys in file 2:"<<endl;
 for (int k=0;k<int(keys2.size());++k){
    cout << "Key "<<k<<": "<<keys2[k].output()<<endl;}
 if ((keys1.size()!=1)||(keys2.size()!=1)){
    cout << "Invalid 4D smeared links file(s)"<<endl;
    return;}
 vector<LattField> sg1,sg2;
 iom1.get(UIntKey(0),sg1);
 iom2.get(UIntKey(0),sg2);
   // check that vectors have at least 3 elements
 if ((sg1.size()<3)||(sg2.size()<3)){
    cout << "Lattice fields read do not have at least 3 elements"<<endl;}
 for (int dir=0;dir<3;++dir){
    cout << "for dir = "<<dir<<" elems per site is "<<sg1[dir].elemsPerSite()<<endl;
    cout << "for dir = "<<dir<<" elems per site is "<<sg2[dir].elemsPerSite()<<endl;
    if  ((sg1[dir].elemsPerSite()!=FieldNcolor*FieldNcolor)
       ||(sg2[dir].elemsPerSite()!=FieldNcolor*FieldNcolor)){
        cout << "improper number of elements per site"<<endl;
        return;}}
 for (int dir=0;dir<3;++dir){
    cout << "for dir = "<<dir<<" bytes per word is "<<sg1[dir].bytesPerWord()<<endl;
    cout << "for dir = "<<dir<<" bytes per word is "<<sg2[dir].bytesPerWord()<<endl;}
 for (int dir=0;dir<3;++dir){
    const double *dp1, *dp2;
    const float *fp1, *fp2;
    bool dp1flag, dp2flag;
    if (sg1[dir].bytesPerWord()==sizeof(complex<double>)){
       dp1flag=true; dp1=reinterpret_cast<const double*>(sg1[dir].getDataConstPtr());}
    else{
       dp1flag=false; fp1=reinterpret_cast<const float*>(sg1[dir].getDataConstPtr());}
    if (sg2[dir].bytesPerWord()==sizeof(complex<double>)){
       dp2flag=true; dp2=reinterpret_cast<const double*>(sg2[dir].getDataConstPtr());}
    else{
       dp2flag=false; fp2=reinterpret_cast<const float*>(sg2[dir].getDataConstPtr());}
    int nreal=2*sg1[dir].elemsPerSite()*LayoutInfo::getLatticeNumSites();
    if (dp1flag&&dp2flag) do_compare(dp1,dp2,nreal,1e-12);
    else if ((!dp1flag)&&(!dp2flag)) do_compare(fp1,fp2,nreal,1e-6);
    else if ((dp1flag)&&(!dp2flag)) do_compare(dp1,fp2,nreal,1e-6);
    else do_compare(dp2,fp1,nreal,1e-6);
    cout <<" Comparisons for direction "<<dir<<" completed"<<endl;
    }

}
*/

template <typename D>
void do_data_compare(const D& data1, const D& data2)
{
 if (data1!=data2) cout <<"The data do not match"<<endl;
}

template <>
void do_data_compare(const LattField& data1, const LattField& data2)
{
 if (data1.elemsPerSite()!=data2.elemsPerSite()){
    cout << "Number of elements per site differs"<<endl;
    return;}
 const double *dp1, *dp2;
 const float *fp1, *fp2;
 bool dp1flag, dp2flag;
 if (data1.bytesPerWord()==sizeof(complex<double>)){
    dp1flag=true; dp1=reinterpret_cast<const double*>(data1.getDataConstPtr());}
 else{
    dp1flag=false; fp1=reinterpret_cast<const float*>(data1.getDataConstPtr());}
 if (data2.bytesPerWord()==sizeof(complex<double>)){
    dp2flag=true; dp2=reinterpret_cast<const double*>(data2.getDataConstPtr());}
 else{
    dp2flag=false; fp2=reinterpret_cast<const float*>(data2.getDataConstPtr());}
 int nreal=2*data1.elemsPerSite()*LayoutInfo::getLatticeNumSites();
 if (dp1flag&&dp2flag) do_vec_compare(dp1,dp2,nreal,1e-10);
 else if ((!dp1flag)&&(!dp2flag)) do_vec_compare(fp1,fp2,nreal,1e-6);
 else if ((dp1flag)&&(!dp2flag)) do_vec_compare(dp1,fp2,nreal,1e-6);
 else do_vec_compare(dp2,fp1,nreal,1e-6);
}

template <>
void do_data_compare(const vector<LattField>& data1, const vector<LattField>& data2)
{
 if (data1.size()!=data2.size()){
     cout << "vector sizes are different: "<<data1.size()<<" vs "<<data2.size()<<endl;
     cout << "...will compare up to minimum size"<<endl;}
 for (uint k=0;k<std::min(data1.size(),data2.size());++k){
     cout << "Comparing element "<<k<<endl;
     do_data_compare(data1[k],data2[k]);}
}



template <typename K, typename D>
void do_compare(IOMap<K,D>& iom1, IOMap<K,D>& iom2)
{
 string header1,header2;
 header1=iom1.getHeader();
 header2=iom2.getHeader();
 cout << "Header 1: "<<header1<<endl;
 cout << "Header 2: "<<header2<<endl;
 set<K> keys1,keys2;
 iom1.getKeys(keys1);
 iom2.getKeys(keys2);
 if (keys1!=keys2){
    cout << "Keys in the two files do not match"<<endl;
    return;}
// cout << "Keys in file 1:"<<endl;
 //for (int k=0;k<int(keys1.size());++k){
 //   cout << "Key "<<k<<": "<<keys1[k].output()<<endl;}
 //cout << "Keys in file 2:"<<endl;
 //for (int k=0;k<int(keys2.size());++k){
 //   cout << "Key "<<k<<": "<<keys2[k].output()<<endl;}
 //if ((keys1.size()!=1)||(keys2.size()!=1)){
 //   cout << "Number of keys in the files do no match"<<endl;
 //   return;}
 for (typename set<K>::const_iterator it=keys1.begin();it!=keys1.end();++it){
    cout <<" Doing comparison for key "<<it->output()<<endl;
    D data1,data2;
    iom1.get(*it,data1);
    iom2.get(*it,data2);
    do_data_compare(data1,data2);}
}
   // check that vectors have at least 3 elements
// if ((sg1.size()<3)||(sg2.size()<3)){
//    cout << "Lattice fields read do not have at least 3 elements"<<endl;}
// for (int dir=0;dir<3;++dir){
//    cout << "for dir = "<<dir<<" elems per site is "<<sg1[dir].elemsPerSite()<<endl;
//    cout << "for dir = "<<dir<<" elems per site is "<<sg2[dir].elemsPerSite()<<endl;
//    if  ((sg1[dir].elemsPerSite()!=FieldNcolor*FieldNcolor)
//       ||(sg2[dir].elemsPerSite()!=FieldNcolor*FieldNcolor)){
//        cout << "improper number of elements per site"<<endl;
//        return;}}
// for (int dir=0;dir<3;++dir){
//    cout << "for dir = "<<dir<<" bytes per word is "<<sg1[dir].bytesPerWord()<<endl;
//    cout << "for dir = "<<dir<<" bytes per word is "<<sg2[dir].bytesPerWord()<<endl;}
// for (int dir=0;dir<3;++dir){
//    const double *dp1, *dp2;
//    const float *fp1, *fp2;
/*    bool dp1flag, dp2flag;
    if (sg1[dir].bytesPerWord()==sizeof(complex<double>)){
       dp1flag=true; dp1=reinterpret_cast<const double*>(sg1[dir].getDataConstPtr());}
    else{
       dp1flag=false; fp1=reinterpret_cast<const float*>(sg1[dir].getDataConstPtr());}
    if (sg2[dir].bytesPerWord()==sizeof(complex<double>)){
       dp2flag=true; dp2=reinterpret_cast<const double*>(sg2[dir].getDataConstPtr());}
    else{
       dp2flag=false; fp2=reinterpret_cast<const float*>(sg2[dir].getDataConstPtr());}
    int nreal=2*sg1[dir].elemsPerSite()*LayoutInfo::getLatticeNumSites();
    if (dp1flag&&dp2flag) do_compare(dp1,dp2,nreal,1e-12);
    else if ((!dp1flag)&&(!dp2flag)) do_compare(fp1,fp2,nreal,1e-6);
    else if ((dp1flag)&&(!dp2flag)) do_compare(dp1,fp2,nreal,1e-6);
    else do_compare(dp2,fp1,nreal,1e-6);
    cout <<" Comparisons for direction "<<dir<<" completed"<<endl;
    }
*/




int main(int argc, const char** argv) 
{
 if (argc!=3){
    cout << "Required command line options:  -i input.xml"<<endl;
    return 0;}
 string inputxmlfile;
 int nargs=argc;
 for (int i=0; i<nargs; ++i){
    string argv_i( argv[i] );
      //  get the input xml file name
    if ((argv_i==std::string("-i"))&&((i+1)<nargs)){
       inputxmlfile=std::string( argv[i+1] );
       ++i;}}
 if (inputxmlfile.empty()){
    cout << "Invalid command line options"<<endl;
    return 1;}
 string cmptype;
 string filename1;
 string filename2;
 int Neigvecs;
 try{
    XMLHandler xml_in;
    xml_in.set_from_file(inputxmlfile);
    if (xml_in.fail()){
       errorLaph("  Unable to read/parse XML content in input file ... exiting...");}
    if (xml_in.get_node_name()!="QudaLaphCompare"){
       errorLaph("  Root tag of input XML must be named \"QudaLaphCompare\" ... exiting...");}
    XMLHandler xml_layoutinfo(xml_in,"LatticeLayoutInfo");
    vector<int> npartitions;
    LayoutInfo::init(xml_layoutinfo,npartitions,true);
    xmlread(xml_in,"Type",cmptype,"QudaLaphCompare");
    if (cmptype=="SmearedGaugeField"){
       xmlread(xml_in,"FileName1",filename1,"QudaLaphCompare");
       xmlread(xml_in,"FileName2",filename2,"QudaLaphCompare");}
    else if (cmptype=="LaphEigenvectors"){
       xmlread(xml_in,"FileStub1",filename1,"QudaLaphCompare");
       xmlread(xml_in,"FileStub2",filename2,"QudaLaphCompare");
       xmlread(xml_in,"NumEigvecs",Neigvecs,"QudaLaphCompare");}}
 catch(const std::exception& xp){
    cout << "Bad input "<<xp.what()<<endl;
    return 1;}

 if (cmptype=="SmearedGaugeField"){
    string gID("Laph--SmearedGaugeField4D");
    IOMap<UIntKey,SmearedGaugeDataType> iomq1;
    IOMap<UIntKey,SmearedGaugeDataType> iomq2;
    try{
       cout <<endl<< "Comparing LapH SmearedGaugeField files"<<endl;
       iomq1.openReadOnly(filename1,gID);
       iomq2.openReadOnly(filename2,gID);
       do_compare(iomq1,iomq2);
       cout <<"If no mismatches were reported, the comparison succeeded"<<endl;}
    catch(const std::exception& xp){
       cout << "Error:  "<<xp.what()<<endl;
       return 0;}}
 else if (cmptype=="LaphEigenvectors"){
    string gID("Laph--SmearedQuarkLevelFile");
    IOMap<UIntKey,LaphEigvecType> iomq1;
    IOMap<UIntKey,LaphEigvecType> iomq2;
    for (int v=0;v<Neigvecs;++v){
       try{
          cout <<endl<< "Comparing LapH Eigenvector "<<v<<" files"<<endl;
          iomq1.openReadOnly(filename1+"_level."+make_string(v),gID);
          iomq2.openReadOnly(filename2+"_level."+make_string(v),gID);
          do_compare(iomq1,iomq2);
          iomq1.close(); iomq2.close();}
       catch(const std::exception& xp){
          cout << "Error:  "<<xp.what()<<endl;
          return 0;}}
     cout <<"If no mismatches were reported, the comparison succeeded"<<endl;}

 return 0;
}

#endif
