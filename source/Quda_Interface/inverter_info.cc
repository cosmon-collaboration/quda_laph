#include "inverter_info.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include "quda_info.h"
#include "util_quda.h"
#include "utils.h"
#include <map>

using namespace std;

namespace LaphEnv {

InverterInfo::InverterInfo(const XMLHandler &xml_in) {
  XMLHandler xmlr(xml_in, "InverterInfo");
  string name;
  xmlread(xmlr, "Name", name, "InverterInfo");
  if (name == "CGNR") {
    set_info_cgnr(xmlr);
  } else if (name == "BICGSTAB") {
    set_info_bicgstab(xmlr);
  } else if (name == "GCR") {
    set_info_gcr(xmlr);
  } else if (name == "GCR_MULTIGRID") {
    set_info_gcr_multigrid(xmlr);
  } else {
    xmlreadfail(xmlr, "InverterInfo", "Unsupported name in InverterInfo");
  }
}

// copy constructor

InverterInfo::InverterInfo(const InverterInfo &rhs)
    : svalues(rhs.svalues), ivalues(rhs.ivalues), rvalues(rhs.rvalues) {}

InverterInfo &InverterInfo::operator=(const InverterInfo &rhs) {
  svalues = rhs.svalues;
  ivalues = rhs.ivalues;
  rvalues = rhs.rvalues;
  return *this;
}

InverterInfo::~InverterInfo() {}

void InverterInfo::checkEqual(const InverterInfo &rhs) const {
  if (!((*this) == rhs)) {
    throw(std::invalid_argument("InverterInfo contents do not match"));
  }
}

bool InverterInfo::operator==(const InverterInfo &rhs) const {
  for (int k = 0; k < int(svalues.size()); ++k) {
    if (svalues[k] != rhs.svalues[k]) {
      return false;
    }
  }
  for (int k = 0; k < int(ivalues.size()); ++k) {
    if (ivalues[k] != rhs.ivalues[k]) {
      return false;
    }
  }
  for (int k = 0; k < int(rvalues.size()); ++k) {
    if (std::abs(rvalues[k] - rhs.rvalues[k]) > 1e-12) {
      return false;
    }
  }
  return true;
}

string InverterInfo::output(int indent) const {
  XMLHandler xmlh;
  output(xmlh);
  return xmlh.output(indent);
}

void InverterInfo::output(XMLHandler &xmlout) const {
  if (svalues[0] == "CGNR") {
    output_cgnr(xmlout);
  } else if (svalues[0] == "BICGSTAB") {
    output_bicgstab(xmlout);
  } else if (svalues[0] == "GCR") {
    output_gcr(xmlout);
  } else if (svalues[0] == "GCR_MULTIGRID") {
    output_gcr_multigrid(xmlout);
  }
}

void InverterInfo::setQudaInvertParam(
    QudaInvertParam &invParam, const QuarkActionInfo &qactioninfo) const {
  invParam = newQudaInvertParam();
  qactioninfo.setQudaInvertParam(invParam);
  if (svalues[0] == "CGNR") {
    setQudaInvertParam_cgnr(invParam);
  } else if (svalues[0] == "BICGSTAB") {
    setQudaInvertParam_bicgstab(invParam);
  } else if (svalues[0] == "GCR") {
    setQudaInvertParam_gcr(invParam);
  } else if (svalues[0] == "GCR_MULTIGRID") {
    setQudaInvertParam_gcr_multigrid(invParam, qactioninfo);
  }
}

// these are always the same
void InverterInfo::commonQudaInvertParam(QudaInvertParam &invParam) const {
  invParam.cpu_prec = QudaInfo::get_cpu_prec();
  invParam.cuda_prec = QudaInfo::get_cuda_prec();
  invParam.solution_type = QUDA_MAT_SOLUTION;
  invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
  invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
  int ivalindex = 0, rvalindex = 0;
  invParam.tol = xmlputQLReal("Tolerance", rvalues, rvalindex);
  invParam.reliable_delta = 0.01;
  invParam.maxiter = xmlputQLInt("MaxIterations", ivalues, ivalindex);
  invParam.pipeline = 0;
  invParam.dagger = QUDA_DAG_NO;
  invParam.verbosity = getVerbosity();
  invParam.compute_true_res = true;
  invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
  invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
}

//
//     Conjugate-Gradient on Normal Equations:
//       (all tags except <Name> optional; default values shown)
//
//    <InvertInfo>
//      <Name>CGNR</Name>
//      <Tolerance>1.0e-10</Tolerance>
//      <MaxIterations>10000</MaxIterations>
//    </InvertInfo>
//
//          rvalues[0]=tolerance in residual
//          ivalues[0]=maximum iterations

void InverterInfo::set_info_cgnr(XMLHandler &xmlr) {
  svalues.resize(1);
  rvalues.resize(1);
  ivalues.resize(1);
  svalues[0] = "CGNR";
  int rvalindex = 0, ivalindex = 0;
  xmlsetQLReal(xmlr, "Tolerance", rvalues, rvalindex, true, 1e-10);
  xmlsetQLInt(xmlr, "MaxIterations", ivalues, ivalindex, true, 10000);
}

void InverterInfo::output_cgnr(XMLHandler &xmlout) const {
  xmlout.set_root("InverterInfo");
  xmlout.put_child("Name", "CGNR");
  int ivalindex = 0, rvalindex = 0;
  xmlout.put_child(xmloutputQLReal("Tolerance", rvalues, rvalindex));
  xmlout.put_child(xmloutputQLInt("MaxIterations", ivalues, ivalindex));
}

void InverterInfo::setQudaInvertParam_cgnr(QudaInvertParam &invParam) const {
  commonQudaInvertParam(invParam);
  invParam.inv_type = QUDA_CGNR_INVERTER;
  invParam.struct_size = sizeof(invParam);
}

//
//     Biconjugate Gradient Stabilized:
//
//       (all tags except <Name> optional; default values shown)
//
//    <InvertInfo>
//      <Name>BICGSTAB</Name>
//      <Tolerance>1.0e-10</Tolerance>
//      <MaxIterations>6000</MaxIterations>
//    </InvertInfo>
//
//          rvalues[0]=tolerance in residual
//          ivalues[0]=maximum iterations
//

void InverterInfo::set_info_bicgstab(XMLHandler &xmlr) {
  svalues.resize(1);
  rvalues.resize(1);
  ivalues.resize(1);
  svalues[0] = "BICGSTAB";
  int rvalindex = 0, ivalindex = 0;
  xmlsetQLReal(xmlr, "Tolerance", rvalues, rvalindex, true, 1e-10);
  xmlsetQLInt(xmlr, "MaxIterations", ivalues, ivalindex, true, 6000);
}

void InverterInfo::output_bicgstab(XMLHandler &xmlout) const {
  xmlout.set_root("InverterInfo");
  xmlout.put_child("Name", "BICGSTAB");
  int ivalindex = 0, rvalindex = 0;
  xmlout.put_child(xmloutputQLReal("Tolerance", rvalues, rvalindex));
  xmlout.put_child(xmloutputQLInt("MaxIterations", ivalues, ivalindex));
}

void InverterInfo::setQudaInvertParam_bicgstab(
    QudaInvertParam &invParam) const {
  commonQudaInvertParam(invParam);
  invParam.inv_type = QUDA_BICGSTAB_INVERTER;
  invParam.struct_size = sizeof(invParam);
}

//
//     Generalized Conjugate Residual (Restarted OrthoDir)
//
//     With restarting, the Krylov subspace can be limited and does
//     not need to keep growing.  This version can be used with
//     variable preconditioning (such as Multigrid).
//
//       (all tags except <Name> optional; default values shown)
//
//    <InvertInfo>
//      <Name>GCR</Name>
//      <Tolerance>1.0e-10</Tolerance>
//      <MaxIterations>5000</MaxIterations>
//      <NKrylov>16</NKrylov>
//    </InvertInfo>
//
//          rvalues[0]=tolerance in residual
//          ivalues[0]=maximum iterations
//          ivalues[1]=NKrylov

void InverterInfo::set_info_gcr(XMLHandler &xmlr) {
  svalues.resize(1);
  rvalues.resize(1);
  ivalues.resize(2);
  svalues[0] = "GCR";
  int rvalindex = 0, ivalindex = 0;
  xmlsetQLReal(xmlr, "Tolerance", rvalues, rvalindex, true, 1e-10);
  xmlsetQLInt(xmlr, "MaxIterations", ivalues, ivalindex, true, 5000);
  xmlsetQLInt(xmlr, "NKrylov", ivalues, ivalindex, true, 16);
}

void InverterInfo::output_gcr(XMLHandler &xmlout) const {
  xmlout.set_root("InverterInfo");
  xmlout.put_child("Name", "GCR");
  int ivalindex = 0, rvalindex = 0;
  xmlout.put_child(xmloutputQLReal("Tolerance", rvalues, rvalindex));
  xmlout.put_child(xmloutputQLInt("MaxIterations", ivalues, ivalindex));
  xmlout.put_child(xmloutputQLInt("NKrylov", ivalues, ivalindex));
}

void InverterInfo::setQudaInvertParam_gcr(QudaInvertParam &invParam) const {
  commonQudaInvertParam(invParam);
  invParam.inv_type = QUDA_GCR_INVERTER;
  invParam.struct_size = sizeof(invParam);
}

//  GCR Inverter with Adaptive Multigrid Preconditioning
//
//    XML input should have the following form:
//       (all tags except initial <Name> optional; default values shown,
//        except as discussed below)
//
//    <InvertInfo>
//      <Name>GCR_MULTIGRID</Name>
//      <Tolerance>1.0e-11</Tolerance>
//      <MaxIterations>200</MaxIterations>
//      <NKrylov>24</NKrylov>
//      <MGPreconditioner>
//         <NumLevels>2</NumLevels>    2,3,4 (2 or 3 usually best)
//         <Level0>
//            <XYZTBlockExtents>4 4 4 4</XYZTBlockExtents>
//            <NullSpaceDim>S</NullSpaceDim>  (S=small 24, L=large 32)
//            <NullSpaceSolveSteps>500</NullSpaceSolveSteps>
//            <CoarseSolverTolerance>0.25</CoarseSolverTolerance>
//            <CoarseSolverMaxIterations>50</CoarseSolverMaxIterations>
//            <NumPreSmooth>0</NumPreSmooth>
//            <NumPostSmooth>4</NumPostSmooth>
//         </Level0>
//         <Level1> ...similar to level 0 </Level1>
//         <Level2>  (the coarsest level)
//            <CoarseSolverTolerance>0.25</CoarseSolverTolerance>
//            <CoarseSolverMaxIterations>50</CoarseSolverMaxIterations>
//            <TRLMDeflation>...</TRLMDeflation>  (default: absent)
//         </Level2>
//         <LoadNullVectors>false</LoadNullVectors>     compute or reload
//         <GenerateAllLevels>true</GenerateAllLevels>  all or level 0 only
//         <PreOrthoNullVectors>true</PreOrthoNullVectors>
//         <PostOrthoNullVectors>true</PostOrthoNullVectors>
//         <SetupMinimizeMemory>false<SetupMinimizeMemory>
//         <RunVerify>false</RunVerify>   (use true for initial runs)
//      </MGPreconditioner>
//    </InvertInfo>
//
//   Input XML for the deflation on coarsest level (if present)
//
//    <TRLMDeflation>
//       <NumVectors> 128 </NumVectors>
//       <Tolerance> 1e-6 </Tolerance>
//       <MaxIterations> 100 </MaxIterations>
//       <ChebyshevOrder> 8 </ChebyshevOrder>
//       <CutoffEigenvalue> 0.1 </CutoffEigenvalue>
//    </TRLMDeflation>
//
//   On an L^3 x T lattice:
//      - default value for <NumLevels> is 2
//      - default <XYZTBlockExtents> values determined based on L,T
//

void InverterInfo::set_info_gcr_multigrid(XMLHandler &xmlr) {
  svalues.resize(1);
  const int mg_max_levels = 4;
  int deflate = 1;
  rvalues.resize(2 + (mg_max_levels - 1) + 2 * deflate);
  ivalues.resize(11 + 9 * (mg_max_levels - 1) + 3 * deflate);
  svalues[0] = "GCR_MULTIGRID";
  int rvalindex = 0, ivalindex = 0;
  xmlsetQLReal(xmlr, "Tolerance", rvalues, rvalindex, true, 1e-11);
  xmlsetQLInt(xmlr, "MaxIterations", ivalues, ivalindex, true, 200);
  xmlsetQLInt(xmlr, "NKrylov", ivalues, ivalindex, true, 24);
  XMLHandler xmlmg(getXML_nofail(xmlr, "MGPreconditioner"));
  vector<int> localextents = LayoutInfo::getRankLattExtents();
  // set default value of nlevels
  int nlevels = 2;
  // read nlevels if given and reset nlevels
  xmlsetQLInt(xmlmg, "NumLevels", ivalues, ivalindex, true, nlevels);
  nlevels = ivalues[ivalindex - 1];
  if ((nlevels < 2) || (nlevels > mg_max_levels)) {
    throw(std::invalid_argument(
        "Unsupported number of levels in MultiGrid inverter"));
  }
  rvalues.resize(2 + (nlevels - 1) + 2 * deflate);
  ivalues.resize(11 + 9 * (nlevels - 1) + 3 * deflate);

  // loop over levels, except for coarsest last level
  for (int level = 0; level < (nlevels - 1); ++level) {
    const string tag = "Level" + make_string(level);
    XMLHandler xmllevel(getXML_nofail(xmlmg, tag));
    vector<int> block_extents(LayoutInfo::Ndim);
    do_blocking(localextents, block_extents); // set default blocking extents
    xmlsetQLIntVector(xmllevel, "XYZTBlockExtents", ivalues, ivalindex,
                      LayoutInfo::Ndim, true, block_extents);
    for (int d = 0; d < LayoutInfo::Ndim; ++d) {
      block_extents[d] = ivalues[ivalindex - LayoutInfo::Ndim + d];
    }
    // check if blocking values are valid
    string errmsg;
    if (!check_blocking(block_extents, localextents, errmsg)) {
      string errmessage =
          string("Invalid level ") + make_string(level) + " block extents";
      for (int d = 0; d < LayoutInfo::Ndim; ++d) {
        errmessage += " " + make_string(block_extents[d]);
      }
      errmessage += " in MultiGrid inverter: " + errmsg;
      throw(std::invalid_argument(errmessage));
    }
    // read remaining parameters
    xmlsetQLEnum(xmllevel, "NullSpaceDim", {"S", "L"}, ivalues, ivalindex, true,
                 0);
    xmlsetQLInt(xmllevel, "NullSpaceSolveSteps", ivalues, ivalindex, true, 500);
    xmlsetQLReal(xmllevel, "CoarseSolverTolerance", rvalues, rvalindex, true,
                 0.25);
    xmlsetQLInt(xmllevel, "CoarseSolverMaxIterations", ivalues, ivalindex, true,
                50);
    xmlsetQLInt(xmllevel, "NumPreSmooth", ivalues, ivalindex, true, 0);
    xmlsetQLInt(xmllevel, "NumPostSmooth", ivalues, ivalindex, true, 4);
    if (xml_tag_count(xmllevel, "TRLMDeflation") > 0) {
      throw(std::invalid_argument("Deflation only allowed on coarsest level"));
    }
  }

  // now probe coarsest last level
  const int level = nlevels - 1;
  string tag = "Level" + make_string(level);
  XMLHandler xmllevel(getXML_nofail(xmlmg, tag));
  xmlsetQLReal(xmllevel, "CoarseSolverTolerance", rvalues, rvalindex, true,
               0.25);
  xmlsetQLInt(xmllevel, "CoarseSolverMaxIterations", ivalues, ivalindex, true,
              50);
  if (xml_tag_count(xmllevel, "TRLMDeflation") == 1) {
    XMLHandler xmldeflate(xmllevel, "TRLMDeflation");
    ivalues[ivalindex] = 1;
    ++ivalindex;
    xmlsetQLReal(xmldeflate, "Tolerance", rvalues, rvalindex, true, 1e-6);
    xmlsetQLInt(xmldeflate, "MaxIterations", ivalues, ivalindex, true, 100);
    xmlsetQLInt(xmldeflate, "NumVectors", ivalues, ivalindex, true, 128);
    xmlsetQLInt(xmldeflate, "ChebyshevOrder", ivalues, ivalindex, true, 8);
    xmlsetQLReal(xmldeflate, "CutoffEigenvalue", rvalues, rvalindex, true, 0.1);
  } else {
    deflate = 0;
    ivalues[ivalindex] = 0;
    ++ivalindex;
    rvalues.resize(2 + (nlevels - 1) + 2 * deflate);
    ivalues.resize(11 + 9 * (nlevels - 1) + 3 * deflate);
  }

  xmlsetQLBool(xmlmg, "LoadNullVectors", ivalues, ivalindex, true, false);
  xmlsetQLBool(xmlmg, "GenerateAllLevels", ivalues, ivalindex, true, true);
  xmlsetQLBool(xmlmg, "PreOrthoNullVectors", ivalues, ivalindex, true, true);
  xmlsetQLBool(xmlmg, "PostOrthoNullVectors", ivalues, ivalindex, true, true);
  xmlsetQLBool(xmlmg, "SetupMinimizeMemory", ivalues, ivalindex, true, false);
  xmlsetQLBool(xmlmg, "RunVerify", ivalues, ivalindex, true, true);
}

void InverterInfo::output_gcr_multigrid(XMLHandler &xmlout) const {
  int ivalindex = 0, rvalindex = 0;
  xmlout.set_root("InverterInfo");
  xmlout.put_child("Name", "GCR_MULTIGRID");
  xmlout.put_child(xmloutputQLReal("Tolerance", rvalues, rvalindex));
  xmlout.put_child(xmloutputQLInt("MaxIterations", ivalues, ivalindex));
  xmlout.put_child(xmloutputQLInt("NKrylov", ivalues, ivalindex));
  XMLHandler xmlmg("MGPreconditioner");
  xmlmg.put_child(xmloutputQLInt("NumLevels", ivalues, ivalindex));
  int nlevels = ivalues[ivalindex - 1];
  for (int level = 0; level < (nlevels - 1); ++level) {
    string tag = "Level";
    tag += make_string(level);
    XMLHandler xmllevel(tag);
    xmllevel.put_child(xmloutputQLIntVector("XYZTBlockExtents", ivalues,
                                            ivalindex, LayoutInfo::Ndim));
    xmllevel.put_child(
        xmloutputQLEnum("NullSpaceDim", {"S", "L"}, ivalues, ivalindex));
    xmllevel.put_child(
        xmloutputQLInt("NullSpaceSolveSteps", ivalues, ivalindex));
    xmllevel.put_child(
        xmloutputQLReal("CoarseSolverTolerance", rvalues, rvalindex));
    xmllevel.put_child(
        xmloutputQLInt("CoarseSolverMaxIterations", ivalues, ivalindex));
    xmllevel.put_child(xmloutputQLInt("NumPreSmooth", ivalues, ivalindex));
    xmllevel.put_child(xmloutputQLInt("NumPostSmooth", ivalues, ivalindex));
    xmlmg.put_child(xmllevel);
  }
  const int level = nlevels - 1;
  string tag = "Level";
  tag += make_string(level);
  XMLHandler xmllevel(tag);
  xmllevel.put_child(
      xmloutputQLReal("CoarseSolverTolerance", rvalues, rvalindex));
  xmllevel.put_child(
      xmloutputQLInt("CoarseSolverMaxIterations", ivalues, ivalindex));
  int deflate = ivalues[ivalindex];
  ivalindex++;
  if (deflate) {
    XMLHandler xmldeflate("TRLMDeflation");
    xmldeflate.put_child(xmloutputQLReal("Tolerance", rvalues, rvalindex));
    xmldeflate.put_child(xmloutputQLInt("MaxIterations", ivalues, ivalindex));
    xmldeflate.put_child(xmloutputQLInt("NumVectors", ivalues, ivalindex));
    xmldeflate.put_child(xmloutputQLInt("ChebyshevOrder", ivalues, ivalindex));
    xmldeflate.put_child(
        xmloutputQLReal("CutoffEigenvalue", rvalues, rvalindex));
    xmllevel.put_child(xmldeflate);
  }
  xmlmg.put_child(xmllevel);
  xmlmg.put_child(xmloutputQLBool("LoadNullVectors", ivalues, ivalindex));
  xmlmg.put_child(xmloutputQLBool("GenerateAllLevels", ivalues, ivalindex));
  xmlmg.put_child(xmloutputQLBool("PreOrthoNullVectors", ivalues, ivalindex));
  xmlmg.put_child(xmloutputQLBool("PostOrthoNullVectors", ivalues, ivalindex));
  xmlmg.put_child(xmloutputQLBool("SetupMinimizeMemory", ivalues, ivalindex));
  xmlmg.put_child(xmloutputQLBool("RunVerify", ivalues, ivalindex));
  xmlout.put_child(xmlmg);
}

void InverterInfo::setQudaInvertParam_gcr_multigrid(
    QudaInvertParam &invParam, const QuarkActionInfo &qactioninfo) const {
  int ivalindex = 0, rvalindex = 0;
  const double outer_tolerance = rvalues[rvalindex];
  ++rvalindex;
  const int outer_maxiter = ivalues[ivalindex];
  ++ivalindex;
  const int outer_Nkrylov = ivalues[ivalindex];
  ++ivalindex;
  const int nlevels = ivalues[ivalindex];
  ++ivalindex;

  vector<vector<int>> blockextents(nlevels, vector<int>(LayoutInfo::Ndim));
  vector<int> nullspacedim(nlevels, 24);
  vector<int> setup_solver_numsteps(nlevels, 100);
  vector<double> coarse_solver_tol(nlevels, 0.25);
  vector<int> coarse_solver_maxiter(nlevels, 100);
  vector<int> presmooth_num(nlevels, 0);
  vector<int> postsmooth_num(nlevels, 4);
  for (int level = 0; level < (nlevels - 1); ++level) {
    for (int d = 0; d < LayoutInfo::Ndim; ++d) {
      blockextents[level][d] = ivalues[ivalindex];
      ++ivalindex;
    }
    nullspacedim[level] = (ivalues[ivalindex] == 0) ? 24 : 32;
    ++ivalindex;
    setup_solver_numsteps[level] = ivalues[ivalindex];
    ++ivalindex;
    coarse_solver_tol[level] = rvalues[rvalindex];
    ++rvalindex;
    coarse_solver_maxiter[level] = ivalues[ivalindex];
    ++ivalindex;
    presmooth_num[level] = ivalues[ivalindex];
    ++ivalindex;
    postsmooth_num[level] = ivalues[ivalindex];
    ++ivalindex;
  }

  const int level = nlevels - 1;
  coarse_solver_tol[level] = rvalues[rvalindex];
  ++rvalindex;
  coarse_solver_maxiter[level] = ivalues[ivalindex];
  ++ivalindex;
  const int coarsest_deflate = ivalues[ivalindex];
  ivalindex++;
  double deflate_tol = 0.0, deflate_cutoff = 0.0;
  int deflate_maxits = 0, deflate_nvectors = 0, deflate_chebyshev = 0;
  if (coarsest_deflate) {
    deflate_tol = rvalues[rvalindex];
    ++rvalindex;
    deflate_maxits = ivalues[ivalindex];
    ++ivalindex;
    deflate_nvectors = ivalues[ivalindex];
    ++ivalindex;
    deflate_chebyshev = ivalues[ivalindex];
    ++ivalindex;
    deflate_cutoff = rvalues[rvalindex];
    ++rvalindex;
  }

  const bool load_null_vectors = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;
  const bool generate_all = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;
  const bool pre_ortho_setup = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;
  const bool post_orth_setup = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;
  const bool minimize_mem = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;
  const bool setup_verify = (ivalues[ivalindex]) ? true : false;
  ++ivalindex;

  double setup_tolerance = 1e-10; // so setup will continue to max iteration
  double smoother_tol = 1e-6;     // so will go to max iteration

  // set the QudaInvertParam
  invParam.cpu_prec = QudaInfo::get_cpu_prec();
  invParam.cuda_prec = QudaInfo::get_cuda_prec();
  invParam.solution_type = QUDA_MAT_SOLUTION;
  invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
  invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
  invParam.inv_type = QUDA_GCR_INVERTER; // the outer solver
  invParam.tol = outer_tolerance;
  invParam.reliable_delta = 0.01;
  invParam.maxiter = outer_maxiter;
  invParam.pipeline = 0;
  invParam.dagger = QUDA_DAG_NO;
  invParam.verbosity = getVerbosity();
  invParam.verbosity_precondition = getVerbosity();
  invParam.compute_true_res = true;
  invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
  invParam.cuda_prec_precondition = QUDA_HALF_PRECISION;
  invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec_sloppy();
  invParam.gcrNkrylov = outer_Nkrylov;
  invParam.input_location = QUDA_CPU_FIELD_LOCATION;
  invParam.output_location = QUDA_CPU_FIELD_LOCATION;
  invParam.inv_type_precondition =
      QUDA_MG_INVERTER; // the crucial preconditioner
  invParam.preconditioner = nullptr;
  invParam.deflation_op = nullptr;
  invParam.eig_param = nullptr;
  invParam.deflate = QUDA_BOOLEAN_FALSE;
  invParam.struct_size = sizeof(invParam);

  // create new structure to hold the auxiliary data, then
  // assign information to the QudaMultigridParam structure
  // needed for the preconditioner set up
  QudaMGInfoPtr.reset(new QudaMGInfo);
  QudaInvertParam &mg_inv_param(
      QudaMGInfoPtr->mg_inv_param);                      // to shorten notation
  QudaMultigridParam &mg_param(QudaMGInfoPtr->mg_param); // to shorten notation
  mg_param.invert_param = &mg_inv_param;

  qactioninfo.setQudaInvertParam(mg_inv_param);
  mg_inv_param.cpu_prec = QudaInfo::get_cpu_prec();
  mg_inv_param.cuda_prec = QudaInfo::get_cuda_prec();
  mg_inv_param.solution_type = QUDA_MAT_SOLUTION;
  mg_inv_param.solve_type = QUDA_DIRECT_SOLVE;
  mg_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  mg_inv_param.inv_type = QUDA_GCR_INVERTER;
  mg_inv_param.gamma_basis =
      QUDA_DEGRAND_ROSSI_GAMMA_BASIS; // Must use this basis here
  mg_inv_param.tol = invParam.tol;
  mg_inv_param.reliable_delta = 0.01;
  mg_inv_param.maxiter = invParam.maxiter;
  mg_inv_param.gcrNkrylov = invParam.gcrNkrylov;
  mg_inv_param.pipeline = 0;
  mg_inv_param.dagger = QUDA_DAG_NO;
  mg_inv_param.verbosity = getVerbosity();
  mg_inv_param.compute_true_res = true;
  mg_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  mg_inv_param.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
  mg_inv_param.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
  mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
  mg_inv_param.cuda_prec_eigensolver = QudaInfo::get_cuda_prec_sloppy();
  mg_inv_param.struct_size = sizeof(mg_inv_param);

  mg_param.n_level = nlevels;
  mg_param.setup_type = QUDA_NULL_VECTOR_SETUP;
  mg_param.compute_null_vector = (load_null_vectors)
                                     ? QUDA_COMPUTE_NULL_VECTOR_NO
                                     : QUDA_COMPUTE_NULL_VECTOR_YES;
  mg_param.generate_all_levels =
      (generate_all) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  mg_param.pre_orthonormalize =
      (pre_ortho_setup) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  mg_param.post_orthonormalize =
      (post_orth_setup) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  mg_param.setup_minimize_memory =
      (minimize_mem) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
  mg_param.run_verify = (setup_verify) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;

  for (int level = 0; level < nlevels; ++level) {

    mg_param.setup_inv_type[level] = QUDA_CGNR_INVERTER;
    mg_param.num_setup_iter[level] = 1;
    mg_param.setup_tol[level] = setup_tolerance;
    mg_param.setup_maxiter[level] = setup_solver_numsteps[level];
    mg_param.setup_maxiter_refresh[level] = 0;
    mg_param.n_vec_batch[level] = 1;

    mg_param.coarse_solver[level] =
        (level < (nlevels - 1)) ? QUDA_GCR_INVERTER : QUDA_CA_GCR_INVERTER;
    mg_param.coarse_grid_solution_type[level] = QUDA_MATPC_SOLUTION;
    mg_param.coarse_solver_tol[level] = coarse_solver_tol[level];
    mg_param.coarse_solver_maxiter[level] = coarse_solver_maxiter[level];

    mg_param.coarse_solver_ca_basis[level] = QUDA_CHEBYSHEV_BASIS;
    mg_param.coarse_solver_ca_basis_size[level] = 4;
    mg_param.coarse_solver_ca_lambda_min[level] = 0.0;
    mg_param.coarse_solver_ca_lambda_max[level] =
        -1.0; // determine max by power method

    for (int dir = 0; dir < LayoutInfo::Ndim; ++dir) {
      mg_param.geo_block_size[level][dir] = blockextents[level][dir];
    }

    mg_param.spin_block_size[level] = (level == 0) ? 2 : 1;
    mg_param.n_vec[level] = nullspacedim[level];
    mg_param.smoother[level] = QUDA_CA_GCR_INVERTER;
    mg_param.smoother_solve_type[level] = QUDA_DIRECT_PC_SOLVE;
    mg_param.smoother_tol[level] = smoother_tol;
    mg_param.smoother_solver_ca_basis[level] = QUDA_POWER_BASIS;
    mg_param.smoother_solver_ca_lambda_min[level] = 0.0;
    mg_param.smoother_solver_ca_lambda_max[level] = -1.0;
    mg_param.nu_pre[level] = presmooth_num[level];
    mg_param.nu_post[level] = postsmooth_num[level];
    mg_param.omega[level] = 1.0;
    mg_param.cycle_type[level] = QUDA_MG_CYCLE_RECURSIVE;
    mg_param.location[level] = QUDA_CUDA_FIELD_LOCATION;
    mg_param.precision_null[level] = QUDA_HALF_PRECISION;
    mg_param.verbosity[level] = QUDA_SUMMARIZE;
    mg_param.global_reduction[level] = QUDA_BOOLEAN_YES;
    mg_param.use_eig_solver[level] = QUDA_BOOLEAN_FALSE;
    mg_param.setup_use_mma[level] = QUDA_BOOLEAN_FALSE;
    mg_param.dslash_use_mma[level] = QUDA_BOOLEAN_FALSE;
    strcpy(mg_param.vec_infile[level], "");
    strcpy(mg_param.vec_outfile[level], "");
    mg_param.n_block_ortho[level] = 2;
    mg_param.block_ortho_two_pass[level] = QUDA_BOOLEAN_TRUE;
    mg_param.transfer_type[level] = QUDA_TRANSFER_AGGREGATE;

    if ((level < (nlevels - 1)) || (!coarsest_deflate)) {
      mg_param.eig_param[level] = nullptr;
    } else {
      initDeflationParam_gcr_multigrid();
      mg_param.use_eig_solver[level] = QUDA_BOOLEAN_TRUE;
      mg_param.vec_load[level] = QUDA_BOOLEAN_FALSE;
      mg_param.vec_store[level] = QUDA_BOOLEAN_FALSE;
      QudaEigParam &coarse_deflate_param(QudaMGInfoPtr->coarse_deflate_param);
      coarse_deflate_param.use_poly_acc =
          (deflate_chebyshev > 1) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
      coarse_deflate_param.poly_deg = deflate_chebyshev;
      coarse_deflate_param.a_min = deflate_cutoff;
      coarse_deflate_param.n_ev = deflate_nvectors + 1;
      coarse_deflate_param.n_kr = (3 * deflate_nvectors) / 2;
      coarse_deflate_param.n_ev_deflate = deflate_nvectors;
      coarse_deflate_param.n_conv = deflate_nvectors + 1;
      coarse_deflate_param.tol = deflate_tol;
      coarse_deflate_param.max_restarts = deflate_maxits;
      coarse_deflate_param.struct_size = sizeof(coarse_deflate_param);
      mg_param.eig_param[level] = &coarse_deflate_param;
    }
  }
  mg_param.struct_size = sizeof(mg_param);
}

// Initializes the QudaEigParam and its needed QudaInvertParam for the
// coarsest level deflation eigensolver.  This sets all parameters
// that don't need to be set by the end user.

void InverterInfo::initDeflationParam_gcr_multigrid() const {
  QudaEigParam &coarse_deflate_param(
      QudaMGInfoPtr->coarse_deflate_param); // to shorten notation
  coarse_deflate_param.invert_param = 0;    // not used
  coarse_deflate_param.a_max = -1;          // determine by power method
  coarse_deflate_param.eig_type = QUDA_EIG_TR_LANCZOS;
  coarse_deflate_param.use_smeared_gauge = false;
  coarse_deflate_param.spectrum = QUDA_SPECTRUM_SR_EIG;
  coarse_deflate_param.use_norm_op = QUDA_BOOLEAN_TRUE;
  coarse_deflate_param.require_convergence = QUDA_BOOLEAN_FALSE;
  strcpy(coarse_deflate_param.vec_infile, "");
  strcpy(coarse_deflate_param.vec_outfile, "");
}

// Factorize the integer "ivalue" into its prime factors,
// returned in increasing order in "factors".  "ivalue"
// cannot exceed 1024.

void InverterInfo::factorize(const uint ivalue, list<uint> &factors) const {
  if (ivalue > 1024) {
    throw(std::invalid_argument("Cannot factorize for value > 1024"));
  }
  const vector<uint> primes({2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37});
  factors.clear();
  uint n = ivalue, k = 0, f = primes[k];
  while (f * f <= n) {
    if (n % f == 0) {
      factors.push_back(f);
      n /= f;
    } else {
      f = primes[++k];
    }
  }
  if (n != 1) {
    factors.push_back(n);
  }
}

// Given current local fine or pre-coarsened lattice sizes in "extents",
// return a suggested block sizes in "block extents".  The block extent in each
// direction must be a divisor of the lattice extent in that direction,
// and the product of block extents cannot exceed the quda limit of 1024.
// Typically, one wants the largest blocking that is possible.

void InverterInfo::do_blocking(const vector<int> &extents,
                               vector<int> &block_extents) const {
  const uint Nd = LayoutInfo::Ndim;
  const uint max_block_val = 864;
  block_extents.resize(Nd);
  for (uint d = 0; d < Nd; ++d)
    block_extents[d] = 1;
  vector<list<uint>> factors(Nd);
  for (uint d = 0; d < Nd; ++d) {
    factorize(extents[d], factors[d]);
  }
  uint d = 0, block_vol = 1, empty = 0;
  for (uint d = 0; d < Nd; ++d) {
    if (factors[d].empty())
      ++empty;
  }
  while (empty < Nd) {
    if (!factors[d].empty()) {
      uint f = factors[d].front();
      if ((f * block_extents[d] <= 8) && (block_vol * f <= max_block_val)) {
        block_extents[d] *= f;
        block_vol *= f;
        factors[d].pop_front();
        if (factors[d].empty())
          ++empty;
      } else {
        factors[d].clear();
        ++empty;
      }
    }
    ++d;
    if (d == Nd)
      d = 0;
  }
}

// checks that the requested blocking on the current level is valid;
// updates "extents" to be ready for checking next level

bool InverterInfo::check_blocking(std::vector<int> &block_extents,
                                  std::vector<int> &extents,
                                  std::string &error_message) const {
  uint prod = 1;
  for (int d = 0; d < LayoutInfo::Ndim; ++d) {
    int nextvalue = extents[d] / block_extents[d];
    if ((nextvalue * block_extents[d]) != (extents[d])) {
      error_message = make_string(extents[d]) + " not divisible by " +
                      make_string(block_extents[d]);
      return false;
    }
    prod *= block_extents[d];
    extents[d] = nextvalue;
  }
  if ((prod < mg_blockprodmin) || (prod > mg_blockprodmax)) {
    error_message = "Block volume " + make_string(prod) + " outside range " +
                    make_string(mg_blockprodmin) + "-" +
                    make_string(mg_blockprodmax);
    return false;
  }
  return true;
}
} // namespace LaphEnv
