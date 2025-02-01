#include "quark_action_info.h"
#include "quda_info.h"

namespace LaphEnv {

// For all fermion actions, we must have
//
//    ivalues[0]= 0 or 1 or 2 or 3 flavor
//                   0 = light u,d  1 = strange
//                   2 = charm      3 = bottom
//    ivalues[1]= 0 or 1   temporal boundary conditions
//                   0 =antiperiodic   1 = periodic
//
//    rvalues[0] =  factor to multiply solution by (if want
//                   normalization different from quda)

QuarkActionInfo::QuarkActionInfo(const XMLHandler &xml_in) {
  XMLHandler xmlr(xml_in, "QuarkActionInfo");
  std::string name;
  xmlread(xmlr, "Name", name, "QuarkActionInfo");
#ifdef LAPH_DOMAIN_WALL
  if (name == "DOMAIN_WALL") {    
    set_info_domain_wall(xmlr);
#else
  if (name == "WILSON_CLOVER") {
    set_info_wilson_clover(xmlr);
#endif
  } else {
    xmlreadfail(xmlr, "QuarkActionInfo", "Unsupported name in QuarkActionInfo");
  }
}

// copy constructor

QuarkActionInfo::QuarkActionInfo(const QuarkActionInfo &rhs)
    : svalues(rhs.svalues), ivalues(rhs.ivalues), rvalues(rhs.rvalues) {}

QuarkActionInfo &QuarkActionInfo::operator=(const QuarkActionInfo &rhs) {
  svalues = rhs.svalues;
  ivalues = rhs.ivalues;
  rvalues = rhs.rvalues;
  return *this;
}

void QuarkActionInfo::checkEqual(const QuarkActionInfo &rhs) const {
  if (!((*this) == rhs)) {
    throw(std::invalid_argument("QuarkActionInfo contents do not match"));
  }
}

bool QuarkActionInfo::operator==(const QuarkActionInfo &rhs) const {
  for (uint k = 0; k < svalues.size(); ++k) {
    if (svalues[k] != rhs.svalues[k]) {
      return false;
    }
  }
  for (uint k = 0; k < ivalues.size(); ++k) {
    if (ivalues[k] != rhs.ivalues[k]) {
      return false;
    }
  }
  for (uint k = 0; k < rvalues.size(); ++k) {
    if (std::abs(rvalues[k] - rhs.rvalues[k]) > 1e-12) {
      return false;
    }
  }
  return true;
}

std::string QuarkActionInfo::output(const int indent) const {
  XMLHandler xmlh;
  output(xmlh);
  return xmlh.output(indent);
}

void QuarkActionInfo::output(XMLHandler &xmlout) const {
#ifdef LAPH_DOMAIN_WALL
  if (svalues[0] == "DOMAIN_WALL") {
    output_domain_wall(xmlout);
  }
#else
  if (svalues[0] == "WILSON_CLOVER") {
    output_wilson_clover(xmlout);
  }
#endif
}

void QuarkActionInfo::setQudaInvertParam(QudaInvertParam &invParam) const {
#ifdef LAPH_DOMAIN_WALL
  if (svalues[0] == "DOMAIN_WALL") {
    setQudaInvertParam_domain_wall(invParam);
#else
  if (svalues[0] == "WILSON_CLOVER") {
    setQudaInvertParam_wilson_clover(invParam);
#endif
  } else {
    exit(1) ;
  }
}

enum { FLAVOR_IDX , TIMEBC_IDX } ;

// get the flavor and temporal bc
void QuarkActionInfo::readflavortbc( XMLHandler &xmlr )
{
  // expanded for charm and bottom because why the hell not?
  std::string flavor ;
  xmlread(xmlr, "Flavor", flavor, "QuarkActionInfo");
  if ((flavor == "light") || (flavor == "ud") || (flavor == "u") ||
      (flavor == "d") || (flavor == "l")) {
    ivalues[ FLAVOR_IDX ] = 0;
  } else if ((flavor == "s") || (flavor == "strange")) {
    ivalues[ FLAVOR_IDX ] = 1;
  } else if ((flavor == "c") || (flavor == "charm")) {
    ivalues[ FLAVOR_IDX ] = 2;
  } else if ((flavor == "b") || (flavor == "bottom")) {
    ivalues[ FLAVOR_IDX ] = 3;
  } else {
    xmlreadfail(xmlr, "QuarkActionInfo", "Invalid flavor");
  }
  ivalues[ TIMEBC_IDX ] = 0;
  std::string timebc;
  if (xmlreadif(xmlr, "TimeBC", timebc, "QuarkActionInfo")) {
    if (timebc == "periodic") {
      ivalues[TIMEBC_IDX] = 1;
    } else if (timebc == "antiperiodic") {
      ivalues[TIMEBC_IDX] = 0;
    } else {
      xmlreadfail(xmlr, "QuarkActionInfo",
                  "Invalid temporal boundary condition");
    }
  }
}

#ifdef LAPH_DOMAIN_WALL
 enum { NORM , DWF_MASS , DWF_LS , DWF_M5 , DWF_B5 , DWF_C5 } ;

 void QuarkActionInfo::set_info_domain_wall(XMLHandler &xmlr) {
  svalues.resize(1);
  ivalues.resize(2);
  rvalues.resize(6);  
  svalues[0] = "DOMAIN_WALL";
  if (xml_tag_count(xmlr, "Mass") == 1) {
    xmlread(xmlr, "Mass", rvalues[ DWF_MASS ], "QuarkActionInfo");
  }
  if (xml_tag_count(xmlr, "Ls") == 1) {
    xmlread(xmlr, "Ls", rvalues[ DWF_LS ], "QuarkActionInfo");
  }  
  if (xml_tag_count(xmlr, "m5") == 1) {
    xmlread(xmlr, "m5", rvalues[ DWF_M5 ], "QuarkActionInfo");
  }
  if (xml_tag_count(xmlr, "b5") == 1) {
    xmlread(xmlr, "b5", rvalues[ DWF_B5 ], "QuarkActionInfo");
  }
  if (xml_tag_count(xmlr, "c5") == 1) {
    xmlread(xmlr, "c5", rvalues[ DWF_C5 ], "QuarkActionInfo");
  }
  readflavortbc( xmlr ) ;
 }

void QuarkActionInfo::output_domain_wall(XMLHandler &xmlout) const {
  xmlout.set_root("QuarkActionInfo");
  std::string flavor ;
  switch( ivalues[0] ) {
  case 0 : flavor = "ud" ; break ;
  case 1 : flavor = "s" ; break ;
  case 2 : flavor = "c" ; break ;
  case 3 : flavor = "b" ; break ;
  }
  xmlout.put_child("Name", "DOMAIN_WALL");
  xmlout.put_child("Flavor", flavor);
  xmlout.put_child("Mass", make_string(rvalues[DWF_MASS]));
  xmlout.put_child("Ls", make_string(rvalues[DWF_LS]));
  xmlout.put_child("m5", make_string(rvalues[DWF_M5]));
  xmlout.put_child("b5", make_string(rvalues[DWF_B5]));
  xmlout.put_child("c5", make_string(rvalues[DWF_C5]));
  std::string timebc = (ivalues[TIMEBC_IDX] == 0) ? "antiperiodic" : "periodic";
  xmlout.put_child("TimeBC", timebc);
}

// Use Dirac-Pauli spin basis
void QuarkActionInfo::setQudaInvertParam_domain_wall(
     QudaInvertParam &invParam) const {
  invParam.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS;
  invParam.dirac_order = QUDA_DIRAC_ORDER;
  invParam.dslash_type = QUDA_MOBIUS_DWF_DSLASH;
  invParam.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  invParam.mass  =  rvalues[DWF_MASS];
  invParam.Ls    =  rvalues[DWF_LS];
  invParam.m5    = -rvalues[DWF_M5];
  invParam.kappa = 1./(2*(4.-rvalues[DWF_M5])+1) ;
  for( int k = 0 ; k < invParam.Ls ; k++ ) {
    invParam.b_5[k] = rvalues[DWF_B5] ;
    invParam.c_5[k] = rvalues[DWF_C5] ;
  }
  invParam.clover_cpu_prec = QudaInfo::get_cpu_prec();
  invParam.clover_cuda_prec = QudaInfo::get_cuda_prec();
  invParam.clover_cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_eigensolver = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_location = QUDA_CUDA_FIELD_LOCATION;
  invParam.struct_size = sizeof(invParam);
}

#else
 
void QuarkActionInfo::set_info_wilson_clover(XMLHandler &xmlr) {
  svalues.resize(1);
  ivalues.resize(2);
  rvalues.resize(7);
  svalues[0] = "WILSON_CLOVER";
  bool calcmass = true;
  bool calckappa = true;
  if (xml_tag_count(xmlr, "Mass") == 1) {
    xmlread(xmlr, "Mass", rvalues[1], "QuarkActionInfo");
    calcmass = false;
  }
  if (xml_tag_count(xmlr, "Kappa") == 1) {
    xmlread(xmlr, "Kappa", rvalues[2], "QuarkActionInfo");
    calckappa = false;
  }
  if ((calcmass) && (calckappa)) {
    xmlreadfail(xmlr, "QuarkActionInfo",
                "At least one of <Mass> or <Kappa> must be presented");
  }
  if (calcmass) {
    rvalues[1] = 0.5 / rvalues[2] - 4.0;
  }
  if (calckappa) {
    rvalues[2] = 0.5 / (rvalues[1] + 4.0);
  }
  // rvalues[0]=2.0*rvalues[2];
  rvalues[0] = 1.0;
  rvalues[3] = 1.0;
  xmlreadif(xmlr, "Anisotropy", rvalues[3], "QuarkActionInfo");
  if (std::abs(rvalues[3] - 1.0) < 1e-12) {
    xmlread(xmlr, "clovCoeff", rvalues[4], "QuarkActionInfo");
    rvalues[5] = rvalues[4];
  } else {
    xmlread(xmlr, "clovCoeffSS", rvalues[4], "QuarkActionInfo");
    xmlread(xmlr, "clovCoeffST", rvalues[5], "QuarkActionInfo");
  }
  rvalues[6] = 1.0;
  xmlreadif(xmlr, "Tadpole", rvalues[6], "QuarkActionInfo");

  readflavortbc( xmlr ) ;
}

void QuarkActionInfo::output_wilson_clover(XMLHandler &xmlout) const {
  xmlout.set_root("QuarkActionInfo");
  std::string flavor ;
  switch( ivalues[0] ) {
  case 0 : flavor = "ud" ; break ;
  case 1 : flavor = "s" ; break ;
  case 2 : flavor = "c" ; break ;
  case 3 : flavor = "b" ; break ;
  }
  xmlout.put_child("Name", "WILSON_CLOVER");
  xmlout.put_child("Flavor", flavor);
  xmlout.put_child("Mass", make_string(rvalues[1]));
  xmlout.put_child("Kappa", make_string(rvalues[2]));
  xmlout.put_child("Anisotropy", make_string(rvalues[3]));
  if (std::abs(rvalues[3] - 1.0) < 1e-12) {
    xmlout.put_child("clovCoeff", make_string(rvalues[4]));
  } else {
    xmlout.put_child("clovCoeffSS", make_string(rvalues[4]));
    xmlout.put_child("clovCoeffST", make_string(rvalues[5]));
  }
  xmlout.put_child("Tadpole", make_string(rvalues[6]));
  std::string timebc = (ivalues[1] == 0) ? "antiperiodic" : "periodic";
  xmlout.put_child("TimeBC", timebc);
}

// Use Dirac-Pauli spin basis

void QuarkActionInfo::setQudaInvertParam_wilson_clover(
    QudaInvertParam &invParam) const {
  invParam.gamma_basis = QUDA_DIRAC_PAULI_GAMMA_BASIS;
  invParam.dirac_order = QUDA_DIRAC_ORDER;
  invParam.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
  invParam.mass_normalization = QUDA_MASS_NORMALIZATION;
  invParam.kappa = rvalues[2];
  invParam.mass = rvalues[1];
  invParam.clover_cpu_prec = QudaInfo::get_cpu_prec();
  invParam.clover_cuda_prec = QudaInfo::get_cuda_prec();
  invParam.clover_cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_refinement_sloppy =
      QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_cuda_prec_eigensolver = QudaInfo::get_cuda_prec_sloppy();
  invParam.clover_order = QUDA_PACKED_CLOVER_ORDER;
  invParam.clover_coeff = rvalues[4] * rvalues[2];
  invParam.clover_rho = 0.0;
  invParam.Ls = 1;
  invParam.compute_clover = 1;
  invParam.compute_clover_inverse = 1;
  invParam.compute_clover_trlog = 1;
  invParam.clover_location = QUDA_CUDA_FIELD_LOCATION;
  invParam.struct_size = sizeof(invParam);
}
#endif
} // namespace LaphEnv
