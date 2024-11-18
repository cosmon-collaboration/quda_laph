#include "laph_eigen_info.h"
#include "layout_info.h"
#include "quda_info.h"
#include <cstring>

using namespace std;

namespace LaphEnv {

// XmlReader constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(const XMLHandler &xmlin) {
  XMLHandler xml_in(xmlin);
  xml_tag_assert(xml_in, "LaphEigenSolverInfo");
  XMLHandler xmlr(xml_in, "LaphEigenSolverInfo");
  extract_info_from_reader(xmlr);
}

void LaphEigenSolverInfo::extract_info_from_reader(XMLHandler &xmlr) {
  xmlread(xmlr, "ResidualTolerance", tolerance, "LaphEigenSolverInfo");
  xmlread(xmlr, "MaxIterations", maxIterations, "LaphEigenSolverInfo");
  xmlread(xmlr, "KrylovDimension", dimKrylov, "LaphEigenSolverInfo");

  chebyshevOrder = 1;
  if (xml_tag_count(xmlr, "ChebyshevOrder") == 1)
    xmlread(xmlr, "ChebyshevOrder", chebyshevOrder, "LaphEigenSolverInfo");
  if (chebyshevOrder < 1)
    chebyshevOrder = 1;

  maxEigenvalue = 15.0;
  if (xml_tag_count(xmlr, "MaxEigenvalue") == 1) {
    xmlread(xmlr, "MaxEigenvalue", maxEigenvalue, "LaphEigenSolverInfo");
    if (maxEigenvalue < 12.0) {
      xmlreadfail(
          xmlr, "LaphEigenSolverInfo",
          "invalid MaxEigenvalue in LaphEigenSolverInfo: should exceed 12");
    }
  } else {
    if (chebyshevOrder > 1) {
      xmlreadfail(
          xmlr, "LaphEigenSolverInfo",
          "when Chebyshev acceleration used, must supply a maximum eigenvalue");
    }
  }

  cutoffEigenvalue = 0.5;
  if (xml_tag_count(xmlr, "CutoffEigenvalue") == 1) {
    xmlread(xmlr, "CutoffEigenvalue", cutoffEigenvalue, "LaphEigenSolverInfo");
    if ((cutoffEigenvalue < 0.03) ||
        (cutoffEigenvalue > (maxEigenvalue - 0.9))) {
      xmlreadfail(xmlr, "LaphEigenSolverInfo",
                  "invalid CutoffEigenvalue in LaphEigenSolverInfo");
    }
  } else {
    if (chebyshevOrder > 1) {
      xmlreadfail(
          xmlr, "LaphEigenSolverInfo",
          "when Chebyshev acceleration used, must supply a cutoff eigenvalue");
    }
  }

  startVector = "equal_components";
  if (xml_tag_count(xmlr, "StartingVectorType") == 1) {
    string svread;
    xmlread(xmlr, "StartingVectorType", svread, "LaphEigenSolverInfo");
    svread = tidyString(svread);
    if (svread == "random")
      startVector = "random";
    else if (svread == "equal_components")
      startVector = "equal_components";
    else {
      xmlreadfail(xmlr, "LaphEigenSolverInfo",
                  "Invalid <StartingVectorType> tag in LaphEigenSolverInfo");
    }
  }

  if (xml_tag_count(xmlr, "CheckSolution") > 0) {
    check_solution = true;
  } else {
    check_solution = false;
  }

  if (xml_tag_count(xmlr, "OutputVerbosity") == 1) {
    xmlread(xmlr, "OutputVerbosity", outputVerbosity, "LaphEigenSolverInfo");
    if (outputVerbosity < 0)
      outputVerbosity = 0;
    if (outputVerbosity > 2)
      outputVerbosity = 2;
  } else
    outputVerbosity = 0;

  if ((maxIterations < 4) || (dimKrylov < 6) || (tolerance <= 0)) {
    xmlreadfail(xmlr, "LaphEigenSolverInfo",
                "invalid input parameters in LaphEigenSolverInfo");
  }
}

// copy constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(const LaphEigenSolverInfo &in)
    : maxIterations(in.maxIterations), dimKrylov(in.dimKrylov),
      tolerance(in.tolerance), chebyshevOrder(in.chebyshevOrder),
      maxEigenvalue(in.maxEigenvalue), cutoffEigenvalue(in.cutoffEigenvalue),
      startVector(in.startVector), check_solution(in.check_solution),
      outputVerbosity(in.outputVerbosity) {}

LaphEigenSolverInfo &
LaphEigenSolverInfo::operator=(const LaphEigenSolverInfo &in) {
  maxIterations = in.maxIterations;
  startVector = in.startVector;
  dimKrylov = in.dimKrylov;
  tolerance = in.tolerance;
  chebyshevOrder = in.chebyshevOrder;
  maxEigenvalue = in.maxEigenvalue;
  cutoffEigenvalue = in.cutoffEigenvalue;
  check_solution = in.check_solution;
  outputVerbosity = in.outputVerbosity;
  return *this;
}

string LaphEigenSolverInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

void LaphEigenSolverInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("LaphEigenSolverInfo");
  xmlout.put_child("ResidualTolerance", make_string(tolerance));
  xmlout.put_child("MaxIterations", make_string(maxIterations));
  xmlout.put_child("KrylovDimension", make_string(dimKrylov));
  if (chebyshevOrder > 1) {
    xmlout.put_child("ChebyshevOrder", make_string(chebyshevOrder));
    xmlout.put_child("MaxEigenvalue", make_string(maxEigenvalue));
    xmlout.put_child("CutoffEigenvalue", make_string(cutoffEigenvalue));
  }
  xmlout.put_child("StartingVectorType", startVector);
  if (check_solution) {
    xmlout.put_child("CheckSolution", "true");
  } else {
    xmlout.put_child("CheckSolution", "false");
  }
  xmlout.put_child("OutputVerbosity", make_string(outputVerbosity));
}

void LaphEigenSolverInfo::setQudaParam(QudaInvertParam &eig_inv_param,
                                       QudaEigParam &eig_param,
                                       const QuarkSmearingInfo &qsmear) const {
  // build eig_inv_param first
  eig_inv_param = newQudaInvertParam();
  eig_inv_param.dslash_type = QUDA_LAPLACE_DSLASH;
  eig_inv_param.inv_type = QUDA_CG_INVERTER;
  eig_inv_param.tol = tolerance;
  eig_inv_param.maxiter = maxIterations;
  eig_inv_param.reliable_delta = 0.1;
  eig_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  eig_inv_param.mass = 0.0; // for massless Laplacian
  eig_inv_param.kappa = 1.0 / (6.0 + eig_inv_param.mass);

  eig_inv_param.dagger = QUDA_DAG_NO;
  eig_inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
  eig_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  eig_inv_param.dirac_order = QUDA_DIRAC_ORDER;

  eig_inv_param.laplace3D = LayoutInfo::Ndim - 1; // spatial laplacian
  eig_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  eig_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
  eig_inv_param.cpu_prec = QudaInfo::get_cpu_prec();
  eig_inv_param.cuda_prec = QudaInfo::get_cuda_prec();
  if (getOutputVerbosity() == 0) {
    eig_inv_param.verbosity = QUDA_SILENT;
  } else if (getOutputVerbosity() == 1) {
    eig_inv_param.verbosity = QUDA_SUMMARIZE;
  } else {
    eig_inv_param.verbosity = QUDA_VERBOSE;
  }

  eig_inv_param.struct_size = sizeof(eig_inv_param);

  // build eig_param next
  eig_param = newQudaEigParam();
  eig_param.invert_param = &eig_inv_param;
  eig_param.eig_type = QUDA_EIG_TR_LANCZOS_3D;
  eig_param.spectrum = QUDA_SPECTRUM_SR_EIG; // -Laplacian eigenvalues all
                                             // positive (SR = smallest, real)
  eig_param.n_conv = qsmear.getNumberOfLaplacianEigenvectors();
  eig_param.max_restarts = maxIterations;
  eig_param.use_smeared_gauge = QUDA_BOOLEAN_TRUE;
  eig_param.max_ortho_attempts = 4; // not used at the moment
  eig_param.tol = tolerance;
  eig_param.n_kr = dimKrylov;
  eig_param.n_ev = eig_param.n_conv;
  eig_param.use_poly_acc =
      (chebyshevOrder > 1) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  eig_param.poly_deg = chebyshevOrder;
  eig_param.a_min = cutoffEigenvalue / 6.0; // due to quda operator definition
  eig_param.a_max = maxEigenvalue / 6.0;    // due to quda operator definition
  eig_param.n_ev_deflate = 0;

  strcpy(eig_param.vec_infile, "");
  strcpy(eig_param.vec_outfile, "");

  eig_param.block_size = 1 ;
  eig_param.require_convergence = QUDA_BOOLEAN_TRUE;

  eig_param.use_norm_op = QUDA_BOOLEAN_FALSE;
  eig_param.use_dagger = QUDA_BOOLEAN_FALSE;
  eig_param.use_pc = QUDA_BOOLEAN_FALSE;
  eig_param.compute_gamma5 = QUDA_BOOLEAN_FALSE;
  eig_param.compute_svd = QUDA_BOOLEAN_FALSE;
  eig_param.use_eigen_qr = QUDA_BOOLEAN_FALSE;
  eig_param.ortho_dim = LayoutInfo::Ndim - 1;
  eig_param.ortho_dim_size_local =
      LayoutInfo::getRankLattExtents()[LayoutInfo::Ndim - 1];

  eig_param.preserve_deflation =
      QUDA_BOOLEAN_FALSE; // do not keep eigenvectors in gpu device memory
  // eig_param.preserve_deflation_space =
  eig_param.preserve_evals = QUDA_BOOLEAN_FALSE;

  eig_param.struct_size = sizeof(eig_param);
}
} // namespace LaphEnv
