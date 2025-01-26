#include "quark_action_info.h"
#include "quda_info.h"

using namespace std;

namespace LaphEnv {

// For all fermion actions, we must have
//
//    ivalues[0]= 0 or 1   flavor
//                   0 = light u,d  1 = strange
//    ivalues[1]= 0 or 1   temporal boundary conditions
//                   0 =antiperiodic   1 = periodic
//
//    rvalues[0] =  factor to multiply solution by (if want
//                   normalization different from quda)

// ************************************************************

QuarkActionInfo::QuarkActionInfo(const XMLHandler& xml_in)
{
 XMLHandler xmlr(xml_in, "QuarkActionInfo");
 string name;
 xmlread(xmlr,"Name",name,"QuarkActionInfo");
 if (name=="WILSON_CLOVER"){
    set_info_wilson_clover(xmlr);}
 else{
    xmlreadfail(xmlr,"QuarkActionInfo","Unsupported name in QuarkActionInfo");}
}


     // copy constructor

QuarkActionInfo::QuarkActionInfo(const QuarkActionInfo& rhs)
                : svalues(rhs.svalues), ivalues(rhs.ivalues), 
                  rvalues(rhs.rvalues) {}

	
QuarkActionInfo& QuarkActionInfo::operator=(const QuarkActionInfo& rhs)
{
 svalues = rhs.svalues;
 ivalues = rhs.ivalues;
 rvalues = rhs.rvalues;
 return *this;
}


void QuarkActionInfo::checkEqual(const QuarkActionInfo& rhs) const
{
 if (!((*this)==rhs)){
    throw(std::invalid_argument("QuarkActionInfo contents do not match"));}
}


bool QuarkActionInfo::operator==(const QuarkActionInfo& rhs) const
{
 for (uint k=0;k<svalues.size();++k){
    if (svalues[k]!=rhs.svalues[k]){
       return false;}}
 for (uint k=0;k<ivalues.size();++k){
    if (ivalues[k]!=rhs.ivalues[k]){
       return false;}}
 for (uint k=0;k<rvalues.size();++k){
    if (std::abs(rvalues[k]-rhs.rvalues[k])>1e-12){
       return false;}}
 return true;
}


string QuarkActionInfo::output(int indent) const
{
 XMLHandler xmlh;
 output(xmlh);
 return xmlh.output(indent);
}


void QuarkActionInfo::output(XMLHandler& xmlout) const
{
 if (svalues[0]=="WILSON_CLOVER"){
    output_wilson_clover(xmlout);}
}


void QuarkActionInfo::setQudaInvertParam(QudaInvertParam& invParam) const
{
 if (svalues[0]=="WILSON_CLOVER"){
    setQudaInvertParam_wilson_clover(invParam);}
}


// ****************************************************************

//    rvalues[0] =  factor to multiply solution by (if want
//                   normalization different from quda)
//    rvalues[1]= mass
//    rvalues[2]= kappa
//    rvalues[3]= anisotropy
//    rvalues[4]= clover coefficient space-space
//    rvalues[5]= clover coefficient space-time
//    rvalues[6]= tadpole coefficient

//    ivalues[0]= 0 or 1   flavor
//                   0 = light u,d  1 = strange
//    ivalues[1]= 0 or 1   temporal boundary conditions
//                   0 =antiperiodic   1 = periodic

void QuarkActionInfo::set_info_wilson_clover(XMLHandler& xmlr)
{
 svalues.resize(1);
 ivalues.resize(2);
 rvalues.resize(7);
 svalues[0]="WILSON_CLOVER";
 bool calcmass=true;
 bool calckappa=true;
 if (xml_tag_count(xmlr,"Mass")==1){
    xmlread(xmlr,"Mass",rvalues[1],"QuarkActionInfo");
    calcmass=false;}
 if (xml_tag_count(xmlr,"Kappa")==1){
    xmlread(xmlr,"Kappa",rvalues[2],"QuarkActionInfo");
    calckappa=false;}
 if ((calcmass)&&(calckappa)){
    xmlreadfail(xmlr,"QuarkActionInfo","At least one of <Mass> or <Kappa> must be presented");}
 if (calcmass){
    rvalues[1]=0.5/rvalues[2]-4.0;}
 if (calckappa){
    rvalues[2]=0.5/(rvalues[1]+4.0);}
// rvalues[0]=2.0*rvalues[2];
 rvalues[0]=1.0;
 rvalues[3]=1.0;
 xmlreadif(xmlr,"Anisotropy",rvalues[3],"QuarkActionInfo");
 if (std::abs(rvalues[3]-1.0)<1e-12){
    xmlread(xmlr,"clovCoeff",rvalues[4],"QuarkActionInfo");
    rvalues[5]=rvalues[4];}
 else{
    xmlread(xmlr,"clovCoeffSS",rvalues[4],"QuarkActionInfo");
    xmlread(xmlr,"clovCoeffST",rvalues[5],"QuarkActionInfo");}
 rvalues[6]=1.0;
 xmlreadif(xmlr,"Tadpole",rvalues[6],"QuarkActionInfo");
 string flavor;
 xmlread(xmlr,"Flavor",flavor,"QuarkActionInfo");
 if ((flavor=="light")||(flavor=="ud")||(flavor=="u")||(flavor=="d")||(flavor=="l")){
    ivalues[0]=0;}
 else if ((flavor=="s")||(flavor=="strange")){
    ivalues[0]=1;}
 else{
    xmlreadfail(xmlr,"QuarkActionInfo","Invalid flavor");}
 ivalues[1]=0;
 string timebc;
 if (xmlreadif(xmlr,"TimeBC",timebc,"QuarkActionInfo")){
    if (timebc=="periodic"){
       ivalues[1]=1;}
    else if (timebc=="antiperiodic"){
       ivalues[1]=0;}
    else{
       xmlreadfail(xmlr,"QuarkActionInfo","Invalid temporal boundary condition");}}
}


void QuarkActionInfo::output_wilson_clover(XMLHandler& xmlout) const
{
 xmlout.set_root("QuarkActionInfo");
 string flavor=(ivalues[0]==0)?"ud":"s";
 xmlout.put_child("Name","WILSON_CLOVER");
 xmlout.put_child("Flavor",flavor);
 xmlout.put_child("Mass",make_string(rvalues[1]));
 xmlout.put_child("Kappa",make_string(rvalues[2]));
 xmlout.put_child("Anisotropy",make_string(rvalues[3]));
 if (std::abs(rvalues[3]-1.0)<1e-12){
    xmlout.put_child("clovCoeff",make_string(rvalues[4]));}
 else{
    xmlout.put_child("clovCoeffSS",make_string(rvalues[4]));
    xmlout.put_child("clovCoeffST",make_string(rvalues[5]));}
 xmlout.put_child("Tadpole",make_string(rvalues[6]));
 string timebc=(ivalues[1]==0)?"antiperiodic":"periodic";
 xmlout.put_child("TimeBC",timebc);
}


    // Use Dirac-Pauli spin basis

void QuarkActionInfo::setQudaInvertParam_wilson_clover(QudaInvertParam& invParam) const
{
 invParam.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS;
// invParam.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
 invParam.dirac_order = QUDA_DIRAC_ORDER;
 invParam.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
 //invParam.dslash_type = QUDA_WILSON_DSLASH;
 //invParam.mass_normalization = QUDA_KAPPA_NORMALIZATION;
 invParam.mass_normalization = QUDA_MASS_NORMALIZATION;
 invParam.kappa = rvalues[2];
 invParam.mass = rvalues[1];
 invParam.clover_cpu_prec = QudaInfo::get_cpu_prec();
 invParam.clover_cuda_prec = QudaInfo::get_cuda_prec();
 invParam.clover_cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.clover_cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.clover_cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.clover_cuda_prec_eigensolver = QudaInfo::get_cuda_prec_sloppy();
 invParam.clover_order = QUDA_PACKED_CLOVER_ORDER;
// invParam.clover_order = QUDA_FLOAT8_CLOVER_ORDER;
 invParam.clover_coeff = rvalues[4]*rvalues[2];
 invParam.clover_rho = 0.0;
// invParam.clover_csw = clover_csw;
// invParam.compute_clover_trlog = compute_clover_trlog ? 1 : 0;
 invParam.Ls = 1;    
 invParam.compute_clover = 1;
 invParam.compute_clover_inverse = 1;
 invParam.compute_clover_trlog = 1;
 invParam.clover_location = QUDA_CUDA_FIELD_LOCATION;
 invParam.struct_size = sizeof(invParam);
}

  // ************************************************************
}
