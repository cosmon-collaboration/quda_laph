#ifdef TESTING
#include "task_tests.h"
#include "xml_handler.h"
#include "laph_stdio.h"
#include "latt_field.h"
#include "layout_info.h"
#include "utils.h"
#include "quda.h"
#include "color_spinor_field.h"
#include "quda_info.h"
#include <vector>

typedef std::complex<double> dcmplx;
typedef std::complex<float>  fcmplx;


using namespace std;
using namespace LaphEnv;
using namespace quda;

namespace QLTestEnv {


void get_spin_transform_matrix(Array<double>& transmat, const std::string& from_to)
{
 double vcp=1.0/sqrt(2.0);
 double vcm=-1.0/sqrt(2.0);
 if (from_to=="GR_to_DP"){
    transmat(0,0)=0.0; transmat(0,1)=vcm; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcp; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcm;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcp; transmat(3,3)=0.0;}
 else if (from_to=="DP_to_GR"){
    transmat(0,0)=0.0; transmat(0,1)=vcp; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcm; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcp;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcm; transmat(3,3)=0.0;}
 else if (from_to=="UK_to_GR"){
    transmat(0,0)=0.0; transmat(0,1)=vcm; transmat(0,2)=0.0; transmat(0,3)=vcm;
    transmat(1,0)=vcp; transmat(1,1)=0.0; transmat(1,2)=vcp; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcm; transmat(2,2)=0.0; transmat(2,3)=vcp;
    transmat(3,0)=vcp; transmat(3,1)=0.0; transmat(3,2)=vcm; transmat(3,3)=0.0;}
 else if (from_to=="GR_to_UK"){
    transmat(0,0)=0.0; transmat(0,1)=vcp; transmat(0,2)=0.0; transmat(0,3)=vcp;
    transmat(1,0)=vcm; transmat(1,1)=0.0; transmat(1,2)=vcm; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=vcp; transmat(2,2)=0.0; transmat(2,3)=vcm;
    transmat(3,0)=vcm; transmat(3,1)=0.0; transmat(3,2)=vcp; transmat(3,3)=0.0;}
 else if ((from_to=="CH_to_GR")||(from_to=="GR_to_CH")){
    transmat(0,0)=0.0; transmat(0,1)=0.0; transmat(0,2)=0.0; transmat(0,3)=-1.0;
    transmat(1,0)=0.0; transmat(1,1)=0.0; transmat(1,2)=1.0; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=1.0; transmat(2,2)=0.0; transmat(2,3)=0.0;
    transmat(3,0)=-1.0;transmat(3,1)=0.0; transmat(3,2)=0.0; transmat(3,3)=0.0;}
 else if (from_to=="UK_to_CH"){
    transmat(0,0)=vcm; transmat(0,1)=0.0; transmat(0,2)=vcp; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=vcm; transmat(1,2)=0.0; transmat(1,3)=vcp;
    transmat(2,0)=vcp; transmat(2,1)=0.0; transmat(2,2)=vcp; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=vcp; transmat(3,2)=0.0; transmat(3,3)=vcp;}
 else if (from_to=="CH_to_UK"){
    transmat(0,0)=vcm; transmat(0,1)=0.0; transmat(0,2)=vcp; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=vcm; transmat(1,2)=0.0; transmat(1,3)=vcp;
    transmat(2,0)=vcp; transmat(2,1)=0.0; transmat(2,2)=vcp; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=vcp; transmat(3,2)=0.0; transmat(3,3)=vcp;}
 else if ((from_to=="DP_to_UK")||(from_to=="UK_to_DP")){
    transmat(0,0)=-1.0;transmat(0,1)=0.0; transmat(0,2)=0.0; transmat(0,3)=0.0;
    transmat(1,0)=0.0; transmat(1,1)=-1.0;transmat(1,2)=0.0; transmat(1,3)=0.0;
    transmat(2,0)=0.0; transmat(2,1)=0.0; transmat(2,2)=1.0; transmat(2,3)=0.0;
    transmat(3,0)=0.0; transmat(3,1)=0.0; transmat(3,2)=0.0; transmat(3,3)=1.0;}
}

template <typename T>
void convertSpinBasis(complex<T>* zGR, const complex<T>* zDP,
                      const Array<double>& transmat)
{
 complex<T>* pL=zGR;
 const complex<T>* pR=zDP;
 for (int c=0;c<3;c++,++pL,++pR){
    complex<T>* ppL=pL;
    for (int s=0;s<4;++s,ppL+=3){
       complex<T> res=complex<T>(0.0,0.0);
       const complex<T>* ppR=pR;
       for (int ss=0;ss<4;++ss,ppR+=3){
          res+=T(transmat(s,ss))*(*ppR);}
       *ppL=res;}}
}

void convertSpinBasis(std::vector<char>& siteGR, const std::vector<char>& siteDP, bool dp,
                      const Array<double>& transmat)
{
 if (dp){
    complex<double>* zGR=reinterpret_cast<complex<double>*>(siteGR.data());
    const complex<double>* zDP=reinterpret_cast<const complex<double>*>(siteDP.data());
    convertSpinBasis<double>(zGR,zDP,transmat);}
 else{
    complex<float>* zGR=reinterpret_cast<complex<float>*>(siteGR.data());
    const complex<float>* zDP=reinterpret_cast<const complex<float>*>(siteDP.data());
    convertSpinBasis<float>(zGR,zDP,transmat);}
}


void convertSpinBasis(LattField& outferm, const LattField& inferm, const std::string& from_to)
{
 bool dp=(inferm.bytesPerWord()==sizeof(complex<double>)); 
 //int nelem_site=inferm.elemsPerSite();
 Array<double> transmat(4,4);
 get_spin_transform_matrix(transmat,from_to);
 std::vector<char> siteDataOut(inferm.bytesPerSite());
 std::vector<int> coord(LayoutInfo::Ndim);
 for (coord[3]=0;coord[3]<LayoutInfo::getLattSizes()[3];++coord[3])
 for (coord[2]=0;coord[2]<LayoutInfo::getLattSizes()[2];++coord[2])
 for (coord[1]=0;coord[1]<LayoutInfo::getLattSizes()[1];++coord[1])
 for (coord[0]=0;coord[0]<LayoutInfo::getLattSizes()[0];++coord[0]){
    std::vector<char> siteDataIn(inferm.getSiteData(coord));
    convertSpinBasis(siteDataOut,siteDataIn,dp,transmat);
    outferm.putSiteData(coord,siteDataOut);}
}

// ****************************************************


void testSpinConversions(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestSpinConversions")==0)
 return;

 printLaph("Creating starting color-spin fields");
 LattField StartFermDP(FieldSiteType::ColorSpinVector);
 LattField StartFermGR(FieldSiteType::ColorSpinVector);
 LattField StartFermUK(FieldSiteType::ColorSpinVector);
 LattField StartFermCH(FieldSiteType::ColorSpinVector);

 LatticeAssigner Q(807374182,0.0,9.0,true);
 Q.assign_field(StartFermDP,"ColorSpinVectorDP");
 Q.reSeed(621739425);
 Q.assign_field(StartFermCH,"ColorSpinVectorCH");
 Q.reSeed(317481003);
 Q.assign_field(StartFermGR,"ColorSpinVectorGR");
 Q.reSeed(461359627);
 Q.assign_field(StartFermUK,"ColorSpinVectorUK");

 printLaph("Allocating space for host converted fields");
 LattField ConvertedFermDP(FieldSiteType::ColorSpinVector);
 LattField ConvertedFermGR(FieldSiteType::ColorSpinVector);
 LattField ConvertedFermUK(FieldSiteType::ColorSpinVector);
 LattField ConvertedFermCH(FieldSiteType::ColorSpinVector);

 printLaph("Creating Quda fields");
 LattField QudaFermDP(FieldSiteType::ColorSpinVector);
 LattField QudaFermGR(FieldSiteType::ColorSpinVector);
 LattField QudaFermUK(FieldSiteType::ColorSpinVector);
 LattField QudaFermCH(FieldSiteType::ColorSpinVector);

 QudaInvertParam inv_param_DP=newQudaInvertParam();
 inv_param_DP.cpu_prec = QudaInfo::get_cpu_prec();
 inv_param_DP.cuda_prec = QudaInfo::get_cuda_prec();
 inv_param_DP.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_DP.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_DP.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS;
 inv_param_DP.dirac_order = QUDA_DIRAC_ORDER;

 QudaInvertParam inv_param_GR=newQudaInvertParam();
 inv_param_GR.cpu_prec = QudaInfo::get_cpu_prec();
 inv_param_GR.cuda_prec = QudaInfo::get_cuda_prec();
 inv_param_GR.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_GR.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_GR.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
 inv_param_GR.dirac_order = QUDA_DIRAC_ORDER;

 QudaInvertParam inv_param_UK=newQudaInvertParam();
 inv_param_UK.cpu_prec = QudaInfo::get_cpu_prec();
 inv_param_UK.cuda_prec = QudaInfo::get_cuda_prec();
 inv_param_UK.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_UK.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_UK.gamma_basis = QUDA_UKQCD_GAMMA_BASIS;
 inv_param_UK.dirac_order = QUDA_DIRAC_ORDER;

 QudaInvertParam inv_param_CH=newQudaInvertParam();
 inv_param_CH.cpu_prec = QudaInfo::get_cpu_prec();
 inv_param_CH.cuda_prec = QudaInfo::get_cuda_prec();
 inv_param_CH.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_CH.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 inv_param_CH.gamma_basis = QUDA_CHIRAL_GAMMA_BASIS;
 inv_param_CH.dirac_order = QUDA_DIRAC_ORDER;

 lat_dim_t x = {LayoutInfo::getRankLattSizes()[0], LayoutInfo::getRankLattSizes()[1], 
                LayoutInfo::getRankLattSizes()[2], LayoutInfo::getRankLattSizes()[3]};

 printLaph("Forming the starting fields for quda");
 void* sfermDPptr=reinterpret_cast<void*>(StartFermDP.getDataPtr());
 ColorSpinorParam sfermDP_param(sfermDPptr, inv_param_DP, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField sqfermDP(sfermDP_param);

 void* sfermGRptr=reinterpret_cast<void*>(StartFermGR.getDataPtr());
 ColorSpinorParam sfermGR_param(sfermGRptr, inv_param_GR, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField sqfermGR(sfermGR_param);

 void* sfermUKptr=reinterpret_cast<void*>(StartFermUK.getDataPtr());
 ColorSpinorParam sfermUK_param(sfermUKptr, inv_param_UK, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField sqfermUK(sfermUK_param);

 void* sfermCHptr=reinterpret_cast<void*>(StartFermCH.getDataPtr());
 ColorSpinorParam sfermCH_param(sfermCHptr, inv_param_CH, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField sqfermCH(sfermCH_param);

 printLaph("Allocating the converted fields for quda");

 void* fermDPptr=reinterpret_cast<void*>(QudaFermDP.getDataPtr());
 ColorSpinorParam fermDP_param(fermDPptr, inv_param_DP, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField qfermDP(fermDP_param);

 void* fermGRptr=reinterpret_cast<void*>(QudaFermGR.getDataPtr());
 ColorSpinorParam fermGR_param(fermGRptr, inv_param_GR, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField qfermGR(fermGR_param);

 void* fermUKptr=reinterpret_cast<void*>(QudaFermUK.getDataPtr());
 ColorSpinorParam fermUK_param(fermUKptr, inv_param_UK, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField qfermUK(fermUK_param);

 void* fermCHptr=reinterpret_cast<void*>(QudaFermCH.getDataPtr());
 ColorSpinorParam fermCH_param(fermCHptr, inv_param_CH, x, false, QUDA_CPU_FIELD_LOCATION);
 ColorSpinorField qfermCH(fermCH_param);
 
    // First, check straight copies

 printLaph("Testing Dirac-Pauli straight copy");
 qfermDP.copy(sqfermDP);
 compare_fields(StartFermDP,QudaFermDP);
 printLaph("Testing UKQCD straight copy");
 qfermUK.copy(sqfermUK);
 compare_fields(StartFermUK,QudaFermUK);
 printLaph("Testing DeGrand-Rossi straight copy");
 qfermGR.copy(sqfermGR);
 compare_fields(StartFermGR,QudaFermGR);
 printLaph("Testing Chiral straight copy");
 qfermCH.copy(sqfermCH);
 compare_fields(StartFermCH,QudaFermCH);

    // now do the conversions

 string from_to("DP_to_GR");
 string from("Dirac-Pauli");
 string to("DeGrand-Rossi");

 printLaph("\n\n  TESTS BETWEEN DIRAC-PAULI AND DEGRAND-ROSSI\n\n");
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermGR,StartFermDP,from_to);
 qfermGR.copy(sqfermDP);
 compare_fields(QudaFermGR,ConvertedFermGR);

 from_to="GR_to_DP";
 from="DeGrand-Rossi"; to="Dirac-Pauli";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermDP,StartFermGR,from_to);
 qfermDP.copy(sqfermGR);
 compare_fields(QudaFermDP,ConvertedFermDP);

 qfermGR.copy(sqfermDP);
 qfermDP.copy(qfermGR);
 compare_fields(QudaFermDP,StartFermDP);
 qfermDP.copy(sqfermGR);
 qfermGR.copy(qfermDP);
 compare_fields(QudaFermGR,StartFermGR);

 printLaph("\n\n  TESTS BETWEEN UKQCD AND DEGRAND-ROSSI\n\n");
 from_to="UK_to_GR";
 from="UKQCD"; to="DeGrand-Rossi";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermGR,StartFermUK,from_to);
 qfermGR.copy(sqfermUK);
 compare_fields(QudaFermGR,ConvertedFermGR);

 from_to="GR_to_UK";
 from="DeGrand-Rossi"; to="UKQCD";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermUK,StartFermGR,from_to);
 qfermUK.copy(sqfermGR);
 compare_fields(QudaFermUK,ConvertedFermUK);

 qfermGR.copy(sqfermUK);
 qfermUK.copy(qfermGR);
 compare_fields(QudaFermUK,StartFermUK);
 qfermUK.copy(sqfermGR);
 qfermGR.copy(qfermUK);
 compare_fields(QudaFermGR,StartFermGR);


 printLaph("\n\n  TESTS BETWEEN CHIRAL AND DEGRAND-ROSSI\n\n");
 from_to="CH_to_GR";
 from="Chiral"; to="DeGrand-Rossi";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermGR,StartFermCH,from_to);
 qfermGR.copy(sqfermCH);
 compare_fields(QudaFermGR,ConvertedFermGR);

 from_to="GR_to_CH";
 from="DeGrand-Rossi"; to="Chiral";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermCH,StartFermGR,from_to);
 qfermCH.copy(sqfermGR);
 compare_fields(QudaFermCH,ConvertedFermCH);

 qfermGR.copy(sqfermCH);
 qfermCH.copy(qfermGR);
 compare_fields(QudaFermCH,StartFermCH);
 qfermCH.copy(sqfermGR);
 qfermGR.copy(qfermCH);
 compare_fields(QudaFermGR,StartFermGR);


 printLaph("\n\n  TESTS BETWEEN UKQCD AND DIRAC-PAULI\n\n");
 from_to="UK_to_DP";
 from="UKQCD"; to="Dirac-Pauli";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermDP,StartFermUK,from_to);
 qfermDP.copy(sqfermUK);
 compare_fields(QudaFermDP,ConvertedFermDP);

 from_to="DP_to_UK";
 from="Dirac-Pauli"; to="UKQCD";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermUK,StartFermDP,from_to);
 qfermUK.copy(sqfermDP);
 compare_fields(QudaFermUK,ConvertedFermUK);

 qfermUK.copy(sqfermDP);
 qfermDP.copy(qfermUK);
 compare_fields(QudaFermDP,StartFermDP);
 qfermDP.copy(sqfermUK);
 qfermUK.copy(qfermDP);
 compare_fields(QudaFermUK,StartFermUK);

 qfermGR.copy(sqfermDP);
 qfermUK.copy(qfermGR);
 convertSpinBasis(ConvertedFermUK,StartFermDP,"DP_to_UK");
 compare_fields(QudaFermUK,ConvertedFermUK);
 qfermGR.copy(sqfermUK);
 qfermDP.copy(qfermGR);
 convertSpinBasis(ConvertedFermDP,StartFermUK,"UK_to_DP");
 compare_fields(QudaFermDP,ConvertedFermDP);

 printLaph("\n\n  TESTS BETWEEN UKQCD AND CHIRAL\n\n");
 from_to="UK_to_CH";
 from="UKQCD"; to="Chiral";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermCH,StartFermUK,from_to);
 qfermCH.copy(sqfermUK);
 compare_fields(QudaFermCH,ConvertedFermCH);

 from_to="CH_to_UK";
 from="Chiral"; to="UKQCD";
 printLaph(make_str("Testing conversion from ",from," to ",to));
 convertSpinBasis(ConvertedFermUK,StartFermCH,from_to);
 qfermUK.copy(sqfermCH);
 compare_fields(QudaFermUK,ConvertedFermUK);

 qfermCH.copy(sqfermUK);
 qfermUK.copy(qfermCH);
 compare_fields(QudaFermUK,StartFermUK);
 qfermUK.copy(sqfermCH);
 qfermCH.copy(qfermUK);
 compare_fields(QudaFermCH,StartFermCH);

 LattField FermTmp(FieldSiteType::ColorSpinVector);
 printLaph("Convert from UK to GR, then GR to CH");
 convertSpinBasis(ConvertedFermGR,StartFermUK,"UK_to_GR");
 convertSpinBasis(FermTmp,ConvertedFermGR,"GR_to_CH");
 compare_fields(FermTmp,ConvertedFermCH);

 printLaph("Convert from CH to GR, then GR to UK");
 convertSpinBasis(ConvertedFermGR,StartFermCH,"CH_to_GR");
 convertSpinBasis(FermTmp,ConvertedFermGR,"GR_to_UK");
 compare_fields(FermTmp,ConvertedFermUK);

}

// ***********************************************
}
#endif
