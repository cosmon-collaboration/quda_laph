#include "inverter_info.h"
#include "quda_info.h"
#include "util_quda.h"
#include "utils.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include <map>

using namespace std;


namespace LaphEnv {




/*  Balint questions:

--  go over each of the inverters...what info needed for each?  when used?

--  MR inverter?

--  inv_param instead mg_param....  how different from  top-level inv_param?

--  steps in each cycle (smoothing,...)

*/








InverterInfo::InverterInfo(const XMLHandler& xml_in)
{ 
 XMLHandler xmlr(xml_in, "InverterInfo");
 string name;
 xmlread(xmlr,"Name",name,"InverterInfo");
 if (name=="CGNR"){
    set_info_cgnr(xmlr);}
 else if (name=="BICGSTAB"){
    set_info_bicgstab(xmlr);}
 else if (name=="GCR"){
    set_info_gcr(xmlr);}
 else if (name=="GCR_MULTIGRID"){
    set_info_gcr_multigrid(xmlr);}
 else{
    xmlreadfail(xmlr,"InverterInfo","Unsupported name in InverterInfo");}
}


     // copy constructor

InverterInfo::InverterInfo(const InverterInfo& rhs)
                : svalues(rhs.svalues), ivalues(rhs.ivalues), 
                  rvalues(rhs.rvalues) {}

InverterInfo& InverterInfo::operator=(const InverterInfo& rhs)
{
 svalues = rhs.svalues;
 ivalues = rhs.ivalues;
 rvalues = rhs.rvalues;
 return *this;
}

InverterInfo::~InverterInfo() {}

void InverterInfo::checkEqual(const InverterInfo& rhs) const
{
 if (!((*this)==rhs)){
    throw(std::invalid_argument("InverterInfo contents do not match"));}
}


bool InverterInfo::operator==(const InverterInfo& rhs) const
{
 for (int k=0;k<int(svalues.size());++k){
    if (svalues[k]!=rhs.svalues[k]){
       return false;}}
 for (int k=0;k<int(ivalues.size());++k){
    if (ivalues[k]!=rhs.ivalues[k]){
       return false;}}
 for (int k=0;k<int(rvalues.size());++k){
    if (std::abs(rvalues[k]-rhs.rvalues[k])>1e-12){
       return false;}}
 return true;
}


string InverterInfo::output(int indent) const
{
 XMLHandler xmlh;
 output(xmlh);
 return xmlh.output(indent);
}


void InverterInfo::output(XMLHandler& xmlout) const
{
 if (svalues[0]=="CGNR"){
    output_cgnr(xmlout);}
 else if (svalues[0]=="BICGSTAB"){
    output_bicgstab(xmlout);}
 else if (svalues[0]=="GCR"){
    output_gcr(xmlout);}
 else if (svalues[0]=="GCR_MULTIGRID"){
    output_gcr_multigrid(xmlout);}
}


void InverterInfo::setQudaInvertParam(QudaInvertParam& invParam,
                                      const QuarkActionInfo& qactioninfo) const
{
 invParam=newQudaInvertParam();
 qactioninfo.setQudaInvertParam(invParam);
 if (svalues[0]=="CGNR"){
    setQudaInvertParam_cgnr(invParam);}
 else if (svalues[0]=="BICGSTAB"){
    setQudaInvertParam_bicgstab(invParam);}
 else if (svalues[0]=="GCR"){
    setQudaInvertParam_gcr(invParam);}
 else if (svalues[0]=="GCR_MULTIGRID"){
    setQudaInvertParam_gcr_multigrid(invParam,qactioninfo);}
}


// **********************************************************************
// *                                                                    *
// *    Conjugate-Gradient on Normal Equations:  (all tags except <Name> optional;          *
// *                          default values shown)                     *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>CGNR</Name>                                              *
// *     <Tolerance>1.0e-10</Tolerance>                                 *
// *     <MaxIterations>10000</MaxIterations>                           *
// *   </InvertInfo>                                                    *
// *                                                                    *
// *         rvalues[0]=tolerance in residual                           *
// *         ivalues[0]=maximum iterations                              *
// *                                                                    *
// **********************************************************************

void InverterInfo::set_info_cgnr(XMLHandler& xmlr)
{
 svalues.resize(1);
 rvalues.resize(1);
 ivalues.resize(1);
 svalues[0]="CGNR";
 int rvalindex=0;
 int ivalindex=0;
 xmlsetQLReal(xmlr,"Tolerance",rvalues,rvalindex,true,1e-10);
 xmlsetQLInt(xmlr,"MaxIterations",ivalues,ivalindex,true,10000);
}


void InverterInfo::output_cgnr(XMLHandler& xmlout) const
{
 xmlout.set_root("InverterInfo");
 xmlout.put_child("Name","CGNR");
 int ivalindex=0;
 int rvalindex=0;
 xmlout.put_child(xmloutputQLReal("Tolerance",rvalues,rvalindex));
 xmlout.put_child(xmloutputQLInt("MaxIterations",ivalues,ivalindex));
}
 
/*
void InverterInfo::setQudaInvertParam_cgnr(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_CGNR_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 0.1;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}
*/

    //  TESTING CA_CGNR

void InverterInfo::setQudaInvertParam_cgnr(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_CA_CGNR_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 0.1;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.gcrNkrylov = 4;
 invParam.ca_basis = QUDA_CHEBYSHEV_BASIS;   // QUDA_POWER_BASIS
 invParam.ca_lambda_min = 0.0;
 invParam.ca_lambda_max = -1.0;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}

// *************************************************************************
// *                                                                       *
// *    Biconjugate Gradient Stabilized:  (all tags except <Name> optional;*
// *                                       default values shown)           *
// *                                                                       *
// *   <InvertInfo>                                                        *
// *     <Name>BICGSTAB</Name>                                             *
// *     <Tolerance>1.0e-10</Tolerance>                                    *
// *     <MaxIterations>6000</MaxIterations>                               *
// *   </InvertInfo>                                                       *
// *                                                                       *
// *         rvalues[0]=tolerance in residual                              *
// *         ivalues[0]=maximum iterations                                 *
// *                                                                       *
// *************************************************************************

void InverterInfo::set_info_bicgstab(XMLHandler& xmlr)
{
 svalues.resize(1);
 rvalues.resize(1);
 ivalues.resize(1);
 svalues[0]="BICGSTAB";
 int rvalindex=0;
 int ivalindex=0;
 xmlsetQLReal(xmlr,"Tolerance",rvalues,rvalindex,true,1e-10);
 xmlsetQLInt(xmlr,"MaxIterations",ivalues,ivalindex,true,6000);
}


void InverterInfo::output_bicgstab(XMLHandler& xmlout) const
{
 xmlout.set_root("InverterInfo");
 xmlout.put_child("Name","BICGSTAB");
 int ivalindex=0;
 int rvalindex=0;
 xmlout.put_child(xmloutputQLReal("Tolerance",rvalues,rvalindex));
 xmlout.put_child(xmloutputQLInt("MaxIterations",ivalues,ivalindex));
}
 

void InverterInfo::setQudaInvertParam_bicgstab(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_BICGSTAB_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 0.1;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}


// **********************************************************************
// *                                                                    *
// *    Generalized Conjugate Residual (GCR)  (all tags except <Name>   *
// *                                           optional; default values *
// *                                           shown)                   *
// *   <InvertInfo>                                                     *
// *     <Name>GCR</Name>                                               *
// *     <Tolerance>1.0e-10</Tolerance>                                 *
// *     <MaxIterations>5000</MaxIterations>                            *
// *     <NKrylov>16</NKrylov>                                          *
// *   </InvertInfo>                                                    *
// *                                                                    *
// *         rvalues[0]=tolerance in residual                           *
// *         ivalues[0]=maximum iterations                              *
// *         ivalues[1]=NKrylov                                         *
// *                                                                    *
// **********************************************************************


/*
CA_GCR  (Communication Avoiding Generalized Conjugate Residual)
      -- needs basis in addition to NKrylov

MR  (Minimal Residual method)   Hermitian only

GMRES
*/

void InverterInfo::set_info_gcr(XMLHandler& xmlr)
{
 svalues.resize(1);
 rvalues.resize(1);
 ivalues.resize(2);
 svalues[0]="GCR";
 int rvalindex=0;
 int ivalindex=0;
 xmlsetQLReal(xmlr,"Tolerance",rvalues,rvalindex,true,1e-10);
 xmlsetQLInt(xmlr,"MaxIterations",ivalues,ivalindex,true,5000);
 xmlsetQLInt(xmlr,"NKrylov",ivalues,ivalindex,true,16);
}


void InverterInfo::output_gcr(XMLHandler& xmlout) const
{
 xmlout.set_root("InverterInfo");
 xmlout.put_child("Name","GCR");
 int ivalindex=0;
 int rvalindex=0;
 xmlout.put_child(xmloutputQLReal("Tolerance",rvalues,rvalindex));
 xmlout.put_child(xmloutputQLInt("MaxIterations",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLInt("NKrylov",ivalues,ivalindex));
}
 
/*
void InverterInfo::setQudaInvertParam_gcr(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_GCR_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 0.1;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.gcrNkrylov = xmlputQLInt("NKrylov",ivalues,ivalindex);
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}
*/

/*   //  TESTING CA_GCR
void InverterInfo::setQudaInvertParam_gcr(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_CA_GCR_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 0.1;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.gcrNkrylov = xmlputQLInt("NKrylov",ivalues,ivalindex);
 invParam.ca_basis = QUDA_CHEBYSHEV_BASIS;
 invParam.ca_lambda_min = 0.0;
 invParam.ca_lambda_max = -1.0;
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}
*/
   //  TESTING MR
void InverterInfo::setQudaInvertParam_gcr(QudaInvertParam& invParam) const
{
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_MR_INVERTER;
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 1e-5;    //  mixed precision parameter (how often
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
 invParam.Nsteps = 2; //xmlputQLInt("MaxIterations",ivalues,ivalindex);  // compute high precision residual
// invParam.gcrNkrylov = xmlputQLInt("NKrylov",ivalues,ivalindex);
// invParam.ca_basis = QUDA_CHEBYSHEV_BASIS;
// invParam.ca_lambda_min = 0.0;
// invParam.ca_lambda_max = -1.0;
// invParam.use_init_guess = QUDA_USE_INIT_GUESS_YES;
// invParam.global_reduction = QUDA_BOOLEAN_TRUE;
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
// invParam.residual_type = QUDA_L2_RELATIVE_RESIDUAL;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.struct_size = sizeof(invParam);
}

// **************************************************************************
// *                                                                        *
// *    GCR Inverter with Adaptive Multigrid Preconditioning                *
// *                                                                        *
// *      Limitations: - only QUDA_MG_CYCLE_RECURSIVE cycle type allowed    *
// *                   - deflation only allowed on coarsest grid            *
// *      (all tags except <Name> optional; default values shown)           *  
// *                                                                        *
// *   <InvertInfo>                                                         *
// *     <Name>GCR_MULTIGRID</Name>                                         *
// *     <Tolerance>1.0e-11</Tolerance>                                     *
// *     <MaxIterations>1000</MaxIterations>                                *
// *     <NKrylov>24</NKrylov>                                              *
// *     <MGPreconditioner>                                                 *
// *        <CycleType>K</CycleType>  (or V,W,F)                            *
// *        <NumLevels>3</NumLevels>    1,2,3,4 (3 usually best)            *
// *        <Level0>                                                        *
// *           <Blocking>4 4 4 4<Blocking>                                  *
// *           <CoarseSolver> ... </CoarseSolver>                           *
// *           <NullSpaceDim>S</NullSpaceDim>  (S=small 24, L=large 32)     *
// *           <OverRelaxOmega>0.8</OverRelaxOmega>     --- part of MR inverter???                    *
// *           <SmootherSolver> ... </SmootherSolver>                       *
// *        </Level0>                                                       *
// *        <Level1> ...similar to level 0 </Level1>                        *
// *        <Level2>  (the coarsest level)                                  *
// *           <CoarseSolver> ... </CoarseSolver>                           *
// *           <TRLMDeflation> ... </TRLMDeflation>                         *
// *        </Level2>                                                       *
// *        <LoadNullVectors>false</LoadNullVectors>     compute or reload  *
// *        <GenerateAllLevels>true</GenerateAllLevels>  all or level 0 only*
// *        <RunVerify>true</RunVerify>                                     *
// *        <PreOrthonormalize>true</PreOrthonormalize>                     *
// *        <PostOrthonormalize>true</PostOrthonormalize>                   *
// *        <SetupMinimizeMemory>true<SetupMinimizeMemory>                  *
// *     </MGPreconditioner>                                                *
// *   </InvertInfo>                                                        *
// *                                                                        *
// *                                                                        *
// *  Allowed XML for the solvers used above
// *
// *   <...Solver>         (Conjugate gradient)                              *
// *     <Name>CG</Name>  (or CGNE or CGNR)                                              *
// *     <Tolerance>0.1</Tolerance>                                  *
// *     <MaxIterations>200</MaxIterations>                           *
// *   </...Solver>                                                    *
// *
// *   <...Solver>          (Biconjugate stabilized)                                      *
// *     <Name>BICGSTAB</Name>                                          *
// *     <Tolerance>0.1</Tolerance>                                  *
// *     <MaxIterations>200</MaxIterations>                           *
// *   </...Solver>                                                    *
// *
// *   <...Solver>          (Generalized conjugate residual)                                           *
// *     <Name>GCR</Name>                                               *
// *     <Tolerance>0.1</Tolerance>                                  *
// *     <MaxIterations>200</MaxIterations>                           *
// *     <NKrylov>16</NKrylov>                                          *
// *   </...Solver>                                                    *
// *
// *   <...Solver>          (Communications avoiding GCR)                                           *
// *     <Name>CA_GCR</Name>                                               *
// *     <Tolerance>0.1</Tolerance>                                  *
// *     <MaxIterations>200</MaxIterations>                           *
// *     <NKrylov>16</NKrylov>                                          *
// *     <CommAvoidBasis>chebyshev</CommAvoidBasis> (or power)
// *   </...Solver>                                                    *
// *
// *   <...Solver>          (Communications avoiding CGNR)                                           *
// *     <Name>CA_CGNR</Name>                      *
// *     <Tolerance>0.1</Tolerance>                                  *
// *     <MaxIterations>200</MaxIterations>                           *
// *     <CommAvoidBasis>chebyshev</CommAvoidBasis> (or power)
// *     <CommAvoidBasisDim>4</CommAvoidBasisDim>
// *   </...Solver>                                                    *
// *


//CGNR

//CGNE

// *  Allowed XML for the Schwarz preconditioning
// *
// *   <Schwarz>
// *      <Type>none</Type>      (default)
// *   </Schwarz>
// *
// *   <Schwarz>
// *      <Type>additive</Type>
// *      <NumLevels>1</NumLevels>    ????? need pre post??
// *   </Schwarz>
// *
// *   <Schwarz>
// *      <Type>multiplicative</Type>
// *      <NumPreLevels>1</NumPreLevels>          ????????
// *      <NumPostLevels>1</NumPostLevels>
// *   </Schwarz>



// *   <TRLMDeflation>                                         *
// *      <Tolerance> 1e-8 </Tolerance>              *
// *      <MaxIterations> 200 </MaxIterations>                       *
// *      <KrylovDimension> 90 </KrylovDimension>                    *
// *      <ChebyshevOrder> 8 </ChebyshevOrder>                         *
// *      <MaxEigenvalue> 15.0 </MaxEigenvalue>       *
// *      <CutoffEigenvalue> 3.4 </CutoffEigenvalue>        *
// *   </TRLMDeflation>                                         *







// *           <BlockingLevel1>4 4 4 4<BlockingLevel1> (none for last level)


// *           <Solver>   QudaInverterType    coarse_solver[QUDA_MAX_MG_LEVEL] 	The solver that wraps around the coarse grid correction and smoother 	Set for all levels except 0. Suggest using QUDA_GCR_INVERTER on all intermediate grids and QUDA_CA_GCR_INVERTER on the bottom. 	--mg-coarse-solver <level gcr/etc.>
// *           <SolverTolerance>0.25</CoarseSolverTolerance>       double coarse_solver_tol[QUDA_MAX_MG_LEVEL]        Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting each level to 0.25 	--mg-coarse-solver-tol <level gcr/etc.>
// *           <SolverMaxIterations>12</CoarseSolverMaxIterations>    int coarse_solver_maxiter[QUDA_MAX_MG_LEVEL]   Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting in the range 8-100 	--mg-coarse-solver-maxiter <level n>
// *           <NullSpaceDim>S L L</NullSpaceDim>   (S=small 24, L=large 32)      n_vec[QUDA_MAX_MG_LEVEL];


// * 
// *        </CoarsenProlongateCycles>
// *        <Smoother>

// *        </Smoother>
// *        <CoarsestGridDeflation>


// *        </CoarsestGridDeflation>
// *        <RunVerify>true</RunVerify>                   (optional: true)   verification checks after set up
// *     </MGPreconditioner>
// *   </InvertInfo>                                                     *






// **********************************************************************
// *                                                                    *
// *    GCR Inverter with Adaptive Multigrid Preconditioning            *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>GCR_MULTIGRID</Name>                                     *
// *     <Tolerance>1.0e-12</Tolerance>                                 *
// *     <MaxIterations>10000</MaxIterations>                           *
// *     <NKrylov>30</NKrylov>
// *     <MGPreconditioner>
// *        <NumLevels>3</NumLevels>    1,2,3,4 (3 usually best) 
// *        <CoarsenProlongateCycles>
// *           <BlockingLevel0>4 4 4 4<BlockingLevel0> 
// *           <BlockingLevel1>4 4 4 4<BlockingLevel1> (none for last level)


// *           <Solver>   QudaInverterType    coarse_solver[QUDA_MAX_MG_LEVEL] 	The solver that wraps around the coarse grid correction and smoother 	Set for all levels except 0. Suggest using QUDA_GCR_INVERTER on all intermediate grids and QUDA_CA_GCR_INVERTER on the bottom. 	--mg-coarse-solver <level gcr/etc.>
// *           <SolverTolerance>0.25</CoarseSolverTolerance>       double coarse_solver_tol[QUDA_MAX_MG_LEVEL]        Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting each level to 0.25 	--mg-coarse-solver-tol <level gcr/etc.>
// *           <SolverMaxIterations>12</CoarseSolverMaxIterations>    int coarse_solver_maxiter[QUDA_MAX_MG_LEVEL]   Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting in the range 8-100 	--mg-coarse-solver-maxiter <level n>
// *           <NullSpaceDim>S L L</NullSpaceDim>   (S=small 24, L=large 32)      n_vec[QUDA_MAX_MG_LEVEL];


// * 
// *        </CoarsenProlongateCycles>
// *        <Smoother>

// *        </Smoother>
// *        <CoarsestGridDeflation>


// *        </CoarsestGridDeflation>
// *        <RunVerify>true</RunVerify>                   (optional: true)   verification checks after set up
// *     </MGPreconditioner>
// *   </InvertInfo>                                                     *



// *     <PreCondGCRNKrylov>20</PreCondGCRNKrylov> (optional: 20)
// *     <NumLevels>3</NumLevels>    1,2,3,4 (3 usually best) 
// *     <ComputeNullVectors>true</ComputeNullVectors> (optional: true)   compute or reload
// *     <GenerateAllLevels>true</GenerateAllLevels>   (optional: true)   all or level 0 only
// *     <RunVerify>true</RunVerify>                   (optional: true)   verification checks after set up
// *     <PreOrthonormalize>true</PreOrthonormalize>   (optional: true)   set up phase
// *     <PostOrthonormalize>true</PostOrthonormalize> (optional: true)   set up phase
// *     <SetupMinimizeMemory>true<SetupMinimizeMemory>(optional: ture)   set up phase










// **********************************************************************
// *                                                                    *
// *    GCR Inverter with Adaptive Multigrid Preconditioning            *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>GCR_MULTIGRID</Name>                                     *
// *     <Tolerance>1.0e-12</Tolerance>                                 *
// *     <MaxIterations>10000</MaxIterations>                           *
// *     <OuterGCRNKrylov>30</OuterGCRNKrylov>   (optional: 20)
// *     <PreCondGCRNKrylov>20</PreCondGCRNKrylov> (optional: 20)
// *     <NumLevels>3</NumLevels>    1,2,3,4 (3 usually best) 
// *     <ComputeNullVectors>true</ComputeNullVectors> (optional: true)   compute or reload
// *     <GenerateAllLevels>true</GenerateAllLevels>   (optional: true)   all or level 0 only
// *     <RunVerify>true</RunVerify>                   (optional: true)   verification checks after set up
// *     <PreOrthonormalize>true</PreOrthonormalize>   (optional: true)   set up phase
// *     <PostOrthonormalize>true</PostOrthonormalize> (optional: true)   set up phase
// *     <SetupMinimizeMemory>true<SetupMinimizeMemory>(optional: ture)   set up phase




    

    
//    QudaBoolean run_low_mode_check;   / ** Whether to run null Vs eigen vector overlap checks once set up is complete * /

   
//    QudaBoolean run_oblique_proj_check;   / ** Whether to run null vector oblique checks once set up is complete * /


    
//    QudaBoolean coarse_guess;   / ** Whether to use and initial guess during coarse grid deflation * /

    
//    QudaBoolean preserve_deflation;   / ** Whether to preserve the deflation space during MG update * /


   
//    QudaBoolean allow_truncation;    / ** Whether or not to let MG coarsening drop improvements, for ex dropping long links in small aggregation dimensions * /

    


// *     <Level>
// *       <IdNum>0</IdNum>
// *          ....
// *     </Level>
// *     <Level>
// *       <IdNum>1</IdNum>
// *          ....
// *     </Level>
// *     <Level>
// *       <IdNum>2</IdNum>
// *          ....
// *     </Level>


// *   </InvertInfo>                                                    *
// *                                                                    *
// *   For each level, input XML is
// *
// *     <Level>
// *       <IdNum>0</IdNum>
// *
// *
// *       <CoarseSolverType>   QudaInverterType    coarse_solver[QUDA_MAX_MG_LEVEL] 	The solver that wraps around the coarse grid correction and smoother 	Set for all levels except 0. Suggest using QUDA_GCR_INVERTER on all intermediate grids and QUDA_CA_GCR_INVERTER on the bottom. 	--mg-coarse-solver <level gcr/etc.>

// *       <CoarseSolutionType>  QudaSolutionType coarse_grid_solution_type[QUDA_MAX_MG_LEVEL] 	The type of residual to send to the next coarse grid, and thus the type of solution to receive back from this coarse grid 	Use QUDA_MATPC_SOLUTION for all levels; if solving the full unpreconditioned system on level 0 then set coarse_grid_solution_type[0]=QUDA_MAT_SOLUTION 	Implicitly depending on --mg-solve-type and --solve-type

// *       <CoarseSolverTolerance>0.25</CoarseSolverTolerance>       double coarse_solver_tol[QUDA_MAX_MG_LEVEL]        Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting each level to 0.25 	--mg-coarse-solver-tol <level gcr/etc.>

// *       <CoarseSolverMaxIterations>12</CoarseSolverMaxIterations>    int coarse_solver_maxiter[QUDA_MAX_MG_LEVEL]   Tolerance for the solver that wraps around the coarse grid correction and smoother 	Suggest setting in the range 8-100 	--mg-coarse-solver-maxiter <level n>

// *       <Blocking>[QUDA_MAX_MG_LEVEL][QUDA_MAX_DIM] 	             int geo_block_size[QUDA_MAX_MG_LEVEL][QUDA_MAX_DIM];Geometric block sizes to use on each level 	Generally favor more aggressive coarsening for first level, rule of thumb is 4^4 	--mg-block-size <level x y z t>

// *       <SpinBlockSize>    int spin_block_size[QUDA_MAX_MG_LEVEL] 	Spin block sizes to use on each level 	2 for level 0, and 1 thereafter 	N/A

// *       <NumNullVectors>24</NumNullVectors>   (24 or 32; 24 is default)      n_vec[QUDA_MAX_MG_LEVEL];

// *       <SmootherType>   QudaInverterType smoother[QUDA_MAX_MG_LEVEL] 	Smoother to use on each level 	Set to QUDA_CA_GCR_INVERTER for each level 	--mg-smoother <level mr/ca-gcr/etc.>

// *       <SmootherTolerance>   double smoother_tol[QUDA_MAX_MG_LEVEL] 	Tolerance to use for the smoother / solver on each level 	Suggest setting each level to 0.25 	--mg-smoother-tol <level resid_tol>

// *       <NumPreSmooth> int nu_pre[QUDA_MAX_MG_LEVEL] 	Number of pre-smoother applications on each level 	Suggest setting to 0 	--mg-nu-pre  <level 1-20>

// *       <NumPostSmooth> int nu_post[QUDA_MAX_MG_LEVEL] 	Number of post-smoother applications on each level 	Suggest setting to 8 	--mg-nu-post <level 1-20>

// *       <OverRelaxOmega> double omega[QUDA_MAX_MG_LEVEL] 	Over/under relaxation factor for the smoother at each level 	Set to 0.8-1.0 	--mg-omega

// *       <SmootherSolveType>    QudaSolveType smoother_solve_type[QUDA_MAX_MG_LEVEL] 	The type of smoother solve to do on each grid (e/o preconditioning or not) 	Suggest setting to QUDA_DIRECT_PC_SOLVE for all levels 	--mg-solve-type <level solve>

// *       <CycleType>   QudaMultigridCycleType cycle_type[QUDA_MAX_MG_LEVEL] 	The type of multigrid cycle to perform at each level 	Always set to QUDA_MG_CYCLE_RECURSIVE (this sets the MG cycles to be a K-cycle which is generally superior to a V-cycle for non-Hermitian systems) 	

// *       <SmootherSchwarzCycle>    int smoother_schwarz_cycle[QUDA_MAX_MG_LEVEL] 	Number of Schwarz cycles to apply 	Experimental, set to 1 for each level 	--mg-schwarz-cycle <level cycle>

/*   
    QudaCABasis coarse_solver_ca_basis[QUDA_MAX_MG_LEVEL];   / ** Basis to use for CA coarse solvers * /

   
    int coarse_solver_ca_basis_size[QUDA_MAX_MG_LEVEL];   / ** Basis size for CA coarse solvers * /

    
    double coarse_solver_ca_lambda_min[QUDA_MAX_MG_LEVEL];  / ** Minimum eigenvalue for Chebyshev CA basis * /

   
    double coarse_solver_ca_lambda_max[QUDA_MAX_MG_LEVEL];    / ** Maximum eigenvalue for Chebyshev CA basis * /


    
    QudaCABasis smoother_solver_ca_basis[QUDA_MAX_MG_LEVEL];  / ** Basis to use for CA smoother solvers * /

    
    double smoother_solver_ca_lambda_min[QUDA_MAX_MG_LEVEL];   / ** Minimum eigenvalue for Chebyshev CA smoother basis * /

    
    double smoother_solver_ca_lambda_max[QUDA_MAX_MG_LEVEL];   / ** Maximum eigenvalue for Chebyshev CA smoother basis * /

    
 <SmootherHaloPrecision>    QudaPrecision smoother_halo_precision[QUDA_MAX_MG_LEVEL];        / ** Precision to use for halo communication in the smoother * /

QudaSchwarzType smoother_schwarz_type[QUDA_MAX_MG_LEVEL] 	Whether to use additive or multiplicative Schwarz preconditioning in the smoother 	Experimental, set to QUDA_INVALID_SCHWARZ for each level unless you know what you're doing 	--mg-schwarz-type <level false/add/mul>

    
    QudaSolutionType coarse_grid_solution_type[QUDA_MAX_MG_LEVEL];   / ** The type of residual to send to the next coarse grid, and thus the type of solution to receive back from this coarse grid * /Use QUDA_MATPC_SOLUTION for all levels; if solving the full unpreconditioned system on level 0 then set coarse_grid_solution_type[0]=QUDA_MAT_SOLUTION

QudaBoolean global_reduction[QUDA_MAX_MG_LEVEL] 	Whether to use global reductions or not for the smoother / solver at each level 	Experimental, keep set to QUDA_BOOLEAN_YES for all levels unless using a Schwarz smoother 	

    
QudaFieldLocation location[QUDA_MAX_MG_LEVEL] 	Location where each level should be done 	Set to QUDA_CUDA_FIELD_LOCATION for all levels 	N/A

   
    QudaFieldLocation setup_location[QUDA_MAX_MG_LEVEL];   / ** Location where the coarse-operator construction will be computedn * /

    
    QudaBoolean use_eig_solver[QUDA_MAX_MG_LEVEL];   / ** Whether to use eigenvectors for the nullspace or, if the coarsest instance deflate* /

QudaPrecision precision_null[QUDA_MAX_MG_LEVEL] 	Precision to store the null-space vectors and preconditioned coarse-link variables 	Use QUDA_HALF_PRECISION for optimal performance 	--prec-null <double/single/half>

   
    int n_block_ortho[QUDA_MAX_MG_LEVEL];   / ** Number of times to repeat Gram-Schmidt in block orthogonalization * /

   
    QudaBoolean block_ortho_two_pass[QUDA_MAX_MG_LEVEL];   / ** Whether to do passes at block orthogonalize in fixed point for improved accuracy * /

QudaVerbosity verbosity[QUDA_MAX_MG_LEVEL] 	Verbosity on each level of the multigrid 		--mg-verbosity <level verb>

    
    QudaBoolean setup_use_mma[QUDA_MAX_MG_LEVEL];   / ** Setup MMA usage on each level of the multigrid * /

  
    QudaBoolean dslash_use_mma[QUDA_MAX_MG_LEVEL];    / ** Dslash MMA usage on each level of the multigrid * /

-----------------------------
INVERSE ITERATIONS

    QudaInverterType setup_inv_type[QUDA_MAX_MG_LEVEL];  / ** Inverter to use in the setup phase * /QUDA_BICGSTAB_INVERTER or QUDA_CGNR_INVERTER generally preferred
   
    int num_setup_iter[QUDA_MAX_MG_LEVEL];   / ** Number of setup iterations * /
    
    double setup_tol[QUDA_MAX_MG_LEVEL];  / ** Tolerance to use in the setup phase * /1e-6 for CGNR/CGNE, 1e-5 otherwise.

   
    int setup_maxiter[QUDA_MAX_MG_LEVEL];   / ** Maximum number of iterations for each setup solver * /500-1000 should work for most systems

  
    int setup_maxiter_refresh[QUDA_MAX_MG_LEVEL];    / ** Maximum number of iterations for refreshing the null-space vectors * /100 as a default

   
    QudaCABasis setup_ca_basis[QUDA_MAX_MG_LEVEL];    / ** Basis to use for CA solver setup * /QUDA_CHEBYSHEV_BASIS for stability

    
    int setup_ca_basis_size[QUDA_MAX_MG_LEVEL];  / ** Basis size for CA solver setup * /   4, higher requires Chebyshev basis

   
    double setup_ca_lambda_min[QUDA_MAX_MG_LEVEL];   / ** Minimum eigenvalue for Chebyshev CA basis * /  0.0 to be safe

    
    double setup_ca_lambda_max[QUDA_MAX_MG_LEVEL];  / ** Maximum eigenvalue for Chebyshev CA basis * /   -1, triggers power iterations

---------------------    
    
    
    QudaBoolean vec_load[QUDA_MAX_MG_LEVEL];  / ** Whether to load the null-space vectors to disk (requires QIO) * /  NO

   
    char vec_infile[QUDA_MAX_MG_LEVEL][256];   / ** Filename prefix where to load the null-space vectors * /

    
    QudaBoolean vec_store[QUDA_MAX_MG_LEVEL];  / ** Whether to store the null-space vectors to disk (requires QIO) * /  NO

    
    char vec_outfile[QUDA_MAX_MG_LEVEL][256];   / ** Filename prefix for where to save the null-space vectors * /

   
    QudaBoolean mg_vec_partfile[QUDA_MAX_MG_LEVEL];   / ** Whether to store the null-space vectors in singlefile or partfile format * /

double mu_factor[QUDA_MAX_MG_LEVEL] 	Multiplicative factor for the mu parameter 	Only applicable for twisted-mass and twisted-clover fermions 	--mg-mu-factor <level factor>

    
    QudaTransferType transfer_type[QUDA_MAX_MG_LEVEL];  / ** Boolean for aggregation type, implies staggered or not * /


    QudaEigParam *eig_param[QUDA_MAX_MG_LEVEL];


*/




// *         rvalues[0]=tolerance in residual                           *
// *         ivalues[0]=maximum iterations                              *
// *         ivalues[1]=NKrylov
// *                                                                    *
/*

    QudaSetupType setup_type = QUDA_NULL_VECTOR_SETUP;  (ONLY CHOICE: web site wrong) / ** Null-space type to use in the setup phase * / QUDA_SETUP_NULL_VECTOR_INVERSE_ITERATIONS, (no support for others?) QUDA_SETUP_NULL_VECTOR_CHEBYSHEV_FILTER, QUDA_SETUP_NULL_VECTOR_EIGENVECTORS, QUDA_SETUP_NULL_VECTOR_TEST_VECTORS, QUDA_SETUP_NULL_VECTOR_RESTRICT_FINE, or QUDA_SETUP_NULL_VECTOR_FREE_FIELD
    QudaBoolean staggered_kd_dagger_approximation;   / ** Whether or not to use the dagger approximation for the KD preconditioned operator * /
    QudaBoolean thin_update_only;   / ** Whether to do a full (false) or thin (true) update in the context of updateMultigridQuda * /

*/




















/*

          //       <invType>QUDA_MULTIGRID_CLOVER_INVERTER</invType>
          //       <CloverParams>
          //         <!--Mass>-0.356557937</Mass-->
          //         <Kappa>0.137232867</Kappa>
          //         <clovCoeff>1.8248654</clovCoeff>
          //         <AnisoParam>
          //           <anisoP>false</anisoP>
          //           <t_dir>3</t_dir>
          //           <xi_0>1</xi_0>
          //           <nu>1</nu>
          //         </AnisoParam>
          //       </CloverParams>
          //       <RsdTarget>1e-11</RsdTarget>
                 <Delta>1.0e-3</Delta>
                 <Pipeline>0</Pipeline>
          //       <MaxIter>1000</MaxIter>
                 <RsdToleranceFactor>100.0</RsdToleranceFactor>
                 <AntiPeriodicT>true</AntiPeriodicT>
                 <SolverType>GCR</SolverType>
                 <Verbose>false</Verbose>
                 <AsymmetricLinop>true</AsymmetricLinop>
                 <CudaReconstruct>RECONS_12</CudaReconstruct>
                 <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>
                 <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>
                 <AxialGaugeFix>false</AxialGaugeFix>
                <AutotuneDslash>true</AutotuneDslash>
                 <MULTIGRIDParams>
                   <Verbosity>true</Verbosity>
                   <Precision>HALF</Precision>
                   <Reconstruct>RECONS_12</Reconstruct>
                   <Blocking>
                     <elem>3 3 3 2</elem>
                     <elem>2 2 2 2</elem>
                     <elem>2 2 2 2</elem>
                   </Blocking>
                   <CoarseSolverType>
                     <elem>GCR</elem>
                     <elem>GCR</elem>
                     <elem>CA_GCR</elem>
                   </CoarseSolverType>
                   <CoarseResidual>0.1 0.1 0.1 0.1</CoarseResidual>
                   <MaxCoarseIterations>12 12 12 8</MaxCoarseIterations>
                   <RelaxationOmegaMG>1.0 1.0 1.0 1.0</RelaxationOmegaMG>
                   <SmootherType>
                     <elem>CA_GCR</elem>
                     <elem>CA_GCR</elem>
                     <elem>CA_GCR</elem>
                     <elem>CA_GCR</elem>
                   </SmootherType>
                   <SmootherTol>0.25 0.25 0.25 0.25</SmootherTol>
                   <SmootherHaloPrecision>
                     <elem>DEFAULT</elem>
                     <elem>HALF</elem>
                     <elem>HALF</elem>
                     <elem>DEFAULT</elem>
                   </SmootherHaloPrecision>
                   <NullVectors>24 32 32</NullVectors>
                   <Pre-SmootherApplications>0 0 0</Pre-SmootherApplications>
                   <Post-SmootherApplications>8 8 8</Post-SmootherApplications>
                   <SubspaceSolver>
                     <elem>CG</elem>
                     <elem>CG</elem>
                     <elem>CG</elem>
                   </SubspaceSolver>
                   <RsdTargetSubspaceCreate>5e-06 5e-06 5e-06</RsdTargetSubspaceCreate>
                   <MaxIterSubspaceCreate>1000 1000 1000</MaxIterSubspaceCreate>
                   <MaxIterSubspaceRefresh>500 500 500</MaxIterSubspaceRefresh>
                   <OuterGCRNKrylov>30</OuterGCRNKrylov>
                   <PrecondGCRNKrylov>20</PrecondGCRNKrylov>
                   <GenerateNullspace>true</GenerateNullspace>
                   <CheckMultigridSetup>false</CheckMultigridSetup>
                   <GenerateAllLevels>true</GenerateAllLevels>
                   <CycleType>MG_RECURSIVE</CycleType>
                   <SchwarzType>ADDITIVE_SCHWARZ</SchwarzType>
                   <RelaxationOmegaOuter>1.0</RelaxationOmegaOuter>
                   <SetupOnGPU>1 1 1</SetupOnGPU>
                 </MULTIGRIDParams>
                 <SubspaceID>subspace_208</SubspaceID>
                 <ThresholdCount>500</ThresholdCount>
              </InvertParam>
 */






// **********************************************************************

void InverterInfo::set_info_gcr_multigrid(XMLHandler& xmlr)
{ /*
 svalues.resize(1);
 rvalues.resize(10);
 ivalues.resize(50);
 svalues[0]="GCR_MULTIGRID";
 int rvalindex=0;
 int ivalindex=0;
 xmlsetQLReal(xmlr,"Tolerance",rvalues,rvalindex,true,1e-10);
 xmlsetQLInt(xmlr,"MaxIterations",ivalues,ivalindex,true,200);
 xmlsetQLInt(xmlr,"NKrylov",ivalues,ivalindex,true,24);
 const int mg_max_levels=4;
 xmlsetQLInt(xmlr,"NumLevels",ivalues,ivalindex);
 uint nlevels=ivalues[ivalindex-1];
 if ((ivalues[1]<1)||(ivalues[1]>mg_max_levels)){
    throw(std::invalid_argument("Unsupported number of levels in MultiGrid inverter"));}
 list<XMLHandler> xmllevels=xmlr.find("Level");
 if (xmllevels.size()!=nlevels){
    throw(std::invalid_argument("Number of <Level> tags does not match <NumLevels> in MultiGrid inverter"));}
 std::map<int,list<XMLHandler>::iterator> lmap;
 for (list<XMLHandler>::iterator it=xmllevels.begin();it!=xmllevels.end();++it){
    xmlsetQLInt(*it,"IdNum",ivalues,ivalindex);
    int idnum=ivalues[ivalindex-1];
    if ((idnum<0)||(idnum>=int(nlevels))){
       throw(std::invalid_argument("Invalid level <IdNum> in MultiGrid inverter"));}
    lmap.insert(make_pair(idnum,it));
    --ivalindex;}
 if (lmap.size()!=nlevels){
    throw(std::invalid_argument("<Level> tags do not all have appropriate <IdNum> values in MultiGrid inverter"));}
 ivalues.resize(2+16*nlevels);
 rvalues.resize(1+3*nlevels);
 for (uint k=0;k<nlevels;++k){
    list<XMLHandler>::iterator it=lmap[k];
    xmlsetQLInt(*it,"IdNum",ivalues,ivalindex);
    xmlsetQLEnum(*it,"CoarseSolverType",{"CG","BICGSTAB","GCR","CA_GCR","CA_CG"},
                 ivalues,ivalindex,true,(k<(nlevels-1))?2:3);
    xmlsetQLEnum(*it,"CoarseSolutionType",{"MAT","MAT_PC"},ivalues,ivalindex);
    xmlsetQLReal(*it,"CoarseSolverTolerance",rvalues,rvalindex,true,0.25);
    xmlsetQLInt(*it,"CoarseSolverMaxIterations",ivalues,ivalindex,true,40);
    xmlsetQLIntVector(*it,"Blocking",ivalues,ivalindex,LayoutInfo::Ndim);
    xmlsetQLInt(*it,"SpinBlockSize",ivalues,ivalindex,true,(k==0)?2:1);
    xmlsetQLInt(*it,"NumNullVectors",ivalues,ivalindex,true,24); // other choice 32
    xmlsetQLEnum(*it,"SmootherType",{"CA_GCR","CA_CG"},ivalues,ivalindex,true,0);
    xmlsetQLEnum(*it,"SmootherSolveType",{"DIRECT_PC","NORMOP_PC"},ivalues,ivalindex,true,0);
    xmlsetQLReal(*it,"SmootherTolerance",rvalues,rvalindex,true,0.25);
    xmlsetQLInt(*it,"NumPreSmooth",ivalues,ivalindex,true,0);
    xmlsetQLInt(*it,"NumPostSmooth",ivalues,ivalindex,true,8);
    xmlsetQLReal(*it,"OverRelaxOmega",rvalues,rvalindex,true,0.8);
    xmlsetQLEnum(*it,"CycleType",{"K","V","F","W"},ivalues,ivalindex,true,0);
    xmlsetQLInt(*it,"SmootherSchwarzCycle",ivalues,ivalindex,true,1);   
    }

 xmlsetQLBool(xmlr,"ComputeNullVectors",ivalues,ivalindex,true,true);
 xmlsetQLBool(xmlr,"GenerateAllLevels",ivalues,ivalindex,true,true);
 xmlsetQLBool(xmlr,"RunVerify",ivalues,ivalindex,true,true);
// printLaph(make_str("rvalindex = ",rvalindex));
// printLaph(make_str("ivalindex = ",ivalindex)); */
}

// *        <Level0>                                                        *
// *           <Blocking>4 4 4 4<Blocking>                                  *
// *           <CoarseSolver> ... </CoarseSolver>                           *
// *           <NullSpaceDim>S</NullSpaceDim>  (S=small 24, L=large 32)     *
// *           <OverRelaxOmega>0.8</OverRelaxOmega>                         *
// *           <SmootherSolver> ... </SmootherSolver>                       *
// *        </Level0>                                                       *
// *        <Level1> ...similar to level 0 </Level1>                        *
// *        <Level2>  (the coarsest level)                                  *
// *           <CoarseSolver> ... </CoarseSolver>                           *
// *           <TRLMDeflation> ... </TRLMDeflation>                         *
// *        </Level2>                                                       *

/*
void InverterInfo::set_info_gcr_multigrid_level(XMLHandler& xmlr, int level,
                                                int& ivalindex, int& rvalindex)
{
 XMLHandler xmllevel(getXML_nofail(xmlr,string("Level")+make_string(level)));
 xmlsetQLIntVector(xmllevel,"Blocking",ivalues,ivalindex,LayoutInfo::Ndim,true,{4,4,4,4});
 set_info_gcr_multigrid_solver(xmlr,"Coarse",2,ivalindex,rvalindex);

 set_info_gcr_multigrid_solver(xmlr,"Smoother",2,ivalindex,rvalindex);

    // "type" below should be "Coarse" or "Smoother"

void InverterInfo::set_info_gcr_multigrid_solver(XMLHandler& xmlr, const string& type,
                                                 int default_name_choice,
                                                 int& ivalindex, int& rvalindex)
{
 XMLHandler xmlsolve(getXML_nofail(xmlr,type+"Solver"));
 xmlsetQLEnum(xmlsolve,"Name",{"CG","BICGSTAB","GCR","CA_GCR","CA_CG"},
              ivalues,ivalindex,true,default_name_choice);
 int namechoice=ivalues[ivalindex-1];
 xmlsetQLReal(xmlsolve,"Tolerance",rvalues,rvalindex,true,0.25);
 xmlsetQLInt(xmlsolve,"MaxIterations",ivalues,ivalindex,true,40);
 if ((namechoice==2)||(namechoice==3)){
    xmlsetQLInt(xmlsolve,"NKrylov",ivalues,ivalindex,true,16);}
 if ((namechoice==3)||(namechoice==4)){
    xmlsetQLInt(xmlsolve,"ChebyshevOrder",ivalues,ivalindex,true,4);}
}

*/



void InverterInfo::output_gcr_multigrid(XMLHandler& xmlout) const
{ /*
 int ivalindex=0;
 int rvalindex=0;
 xmlout.set_root("InverterInfo");
 xmlout.put_child("Name","GCR_MULTIGRID");
 xmlout.put_child(xmloutputQLReal("Tolerance",rvalues,rvalindex));
 xmlout.put_child(xmloutputQLInt("MaxIterations",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLInt("OuterGCRNKrylov",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLInt("PreCondGCRNKrylov",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLInt("NumLevels",ivalues,ivalindex));
 for (int idnum=0;idnum<ivalues[1];++idnum){
    XMLHandler xmllevel("Level");
    xmllevel.put_child(xmloutputQLInt("IdNum",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLEnum("CoarseSolverType",{"CG","BICGSTAB","GCR","CA_GCR","CA_CG"},ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLEnum("CoarseSolutionType",{"MAT","MAT_PC"},ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLReal("CoarseSolverTolerance",rvalues,rvalindex));
    xmllevel.put_child(xmloutputQLInt("CoarseSolverMaxIterations",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLIntVector("Blocking",ivalues,ivalindex,LayoutInfo::Ndim));
    xmllevel.put_child(xmloutputQLInt("SpinBlockSize",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLInt("NumNullVectors",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLEnum("SmootherType",{"CA_GCR","CA_CG"},ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLEnum("SmootherSolveType",{"DIRECT_PC","NORMOP_PC"},ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLReal("SmootherTolerance",rvalues,rvalindex));
    xmllevel.put_child(xmloutputQLInt("NumPreSmooth",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLInt("NumPostSmooth",ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLReal("OverRelaxOmega",rvalues,rvalindex));
    xmllevel.put_child(xmloutputQLEnum("CycleType",{"R","V","F","W"},ivalues,ivalindex));
    xmllevel.put_child(xmloutputQLInt("SmootherSchwarzCycle",ivalues,ivalindex));   


    xmlout.put_child(xmllevel);
    }
 xmlout.put_child(xmloutputQLBool("ComputeNullVectors",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLBool("GenerateAllLevels",ivalues,ivalindex));
 xmlout.put_child(xmloutputQLBool("RunVerify",ivalues,ivalindex)); */
}
 
/*
void InverterInfo::output_gcr_multigrid_solver(XMLHandler& xmlout, const string& type,
                                               int& ivalindex, int& rvalindex)
{
 XMLHandler xmlsolve(type+"Solver");
 xmlsolve.put_child(xmloutputQLEnum("Name",{"CG","BICGSTAB","GCR","CA_GCR"},ivalues,ivalindex));
 int namechoice=ivalues[ivalindex-1];
 xmlsolve.put_child(xmloutputQLReal("Tolerance",rvalues,rvalindex));
 xmlsolve.put_child(xmloutputQLInt("MaxIterations",ivalues,ivalindex));
 if (namechoice>1){
    xmlsolve.put_child(xmloutputQLInt("NKrylov",ivalues,ivalindex));}
 if (namechoice>2){
    xmlsolve.put_child(xmloutputQLInt("ChebyshevOrder",ivalues,ivalindex));}
 xmlout.put_child(xmlsolve);
}
*/

void InverterInfo::setQudaInvertParam_gcr_multigrid(QudaInvertParam& invParam,
                                         const QuarkActionInfo& qactioninfo) const
{ /*
    // set the QudaInvertParam
 invParam.cpu_prec = QudaInfo::get_cpu_prec();
 invParam.cuda_prec = QudaInfo::get_cuda_prec();
 invParam.solution_type = QUDA_MAT_SOLUTION;
 invParam.solve_type = QUDA_DIRECT_PC_SOLVE;
 invParam.matpc_type = QUDA_MATPC_EVEN_EVEN;
 invParam.tune = QUDA_TUNE_YES;
 invParam.inv_type = QUDA_GCR_INVERTER;    // the outer solver
 int ivalindex=0;
 int rvalindex=0;
 invParam.tol = xmlputQLReal("Tolerance",rvalues,rvalindex);
 invParam.reliable_delta = 1e-5;  
 invParam.maxiter = xmlputQLInt("MaxIterations",ivalues,ivalindex);
 invParam.pipeline = 0;
 invParam.dagger = QUDA_DAG_NO;
 invParam.verbosity = getVerbosity();
 invParam.compute_true_res=true;
 invParam.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 invParam.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 invParam.cuda_prec_precondition = QUDA_HALF_PRECISION;
 invParam.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 invParam.gcrNkrylov = xmlputQLInt("OuterGCRNKrylov",ivalues,ivalindex); //   CHECK???

 invParam.inv_type_precondition=QUDA_MG_INVERTER;   // the crucial preconditioner

    / *
     * The following parameters are related to the solver
     * preconditioner, if enabled.
     * /

    / **
     * The inner Krylov solver used in the preconditioner.  Set to
     * QUDA_INVALID_INVERTER to disable the preconditioner entirely.
     * /
    QudaInverterType inv_type_precondition;

    / ** Preconditioner instance, e.g., multigrid * /
    void *preconditioner;

    / ** Deflation instance * /
    void *deflation_op;

    / ** defines deflation * /
    void *eig_param;

    / ** If true, deflate the initial guess * /
    QudaBoolean deflate;

typedef enum QudaBoolean_s {
  QUDA_BOOLEAN_FALSE = 0,
  QUDA_BOOLEAN_TRUE = 1,
  QUDA_BOOLEAN_INVALID = QUDA_INVALID_ENUM
} QudaBoolean;



    / ** Dirac Dslash used in preconditioner * /
    QudaDslashType dslash_type_precondition;


    / ** Verbosity of the inner Krylov solver * /
    QudaVerbosity verbosity_precondition;

    / ** Tolerance in the inner solver * /
    double tol_precondition;

    / ** Maximum number of iterations allowed in the inner solver * /
    int maxiter_precondition;

    / ** Relaxation parameter used in GCR-DD (default = 1.0) * /
    double omega;

    / ** Basis for CA algorithms * /
    QudaCABasis ca_basis;

typedef enum QudaCABasis_s {
  QUDA_POWER_BASIS,
  QUDA_CHEBYSHEV_BASIS,
  QUDA_INVALID_BASIS = QUDA_INVALID_ENUM
} QudaCABasis;


    / ** Minimum eigenvalue for Chebyshev CA basis * /
    double ca_lambda_min;

    / ** Maximum eigenvalue for Chebyshev CA basis * /
    double ca_lambda_max;

    / ** Basis for CA algorithms in a preconditioned solver * /
    QudaCABasis ca_basis_precondition;

    / ** Minimum eigenvalue for Chebyshev CA basis in a preconditioner solver * /
    double ca_lambda_min_precondition;

    / ** Maximum eigenvalue for Chebyshev CA basis in a preconditioner solver * /
    double ca_lambda_max_precondition;

    / ** Number of preconditioner cycles to perform per iteration * /
    int precondition_cycle;

    / ** Whether to use additive or multiplicative Schwarz preconditioning * /
    QudaSchwarzType schwarz_type;

typedef enum QudaSchwarzType_s {
  QUDA_ADDITIVE_SCHWARZ = 0,
  QUDA_MULTIPLICATIVE_SCHWARZ = 1,
  QUDA_INVALID_SCHWARZ = QUDA_INVALID_ENUM
} QudaSchwarzType;



    / ** The type of accelerator type to use for preconditioner * /
    QudaAcceleratorType accelerator_type_precondition;

typedef enum QudaAcceleratorType_s {
  QUDA_MADWF_ACCELERATOR = 0, // Use the MADWF accelerator
  QUDA_INVALID_ACCELERATOR = QUDA_INVALID_ENUM
} QudaAcceleratorType;



    / **
     * The following parameters are the ones used to perform the adaptive MADWF in MSPCG
     * See section 3.3 of [arXiv:2104.05615]
     * /

    / ** The diagonal constant to suppress the low modes when performing 5D transfer * /
    double madwf_diagonal_suppressor;

    / ** The target MADWF Ls to be used in the accelerator * /
    int madwf_ls;

    / ** The minimum number of iterations after which to generate the null vectors for MADWF * /
    int madwf_null_miniter;

    / ** The maximum tolerance after which to generate the null vectors for MADWF * /
    double madwf_null_tol;

    / ** The maximum number of iterations for the training iterations * /
    int madwf_train_maxiter;

    / ** Whether to load the MADWF parameters from the file system * /
    QudaBoolean madwf_param_load;

    / ** Whether to save the MADWF parameters to the file system * /
    QudaBoolean madwf_param_save;

    / ** Path to load from the file system * /
    char madwf_param_infile[256];

    / ** Path to save to the file system * /
    char madwf_param_outfile[256];

    / **
     * Whether to use the L2 relative residual, Fermilab heavy-quark
     * residual, or both to determine convergence.  To require that both
     * stopping conditions are satisfied, use a bitwise OR as follows:
     *
     * p.residual_type = (QudaResidualType) (QUDA_L2_RELATIVE_RESIDUAL
     *                                     | QUDA_HEAVY_QUARK_RESIDUAL);
     * /
    QudaResidualType residual_type;

typedef enum QudaResidualType_s {
  QUDA_L2_RELATIVE_RESIDUAL = 1, // L2 relative residual [default]
  QUDA_L2_ABSOLUTE_RESIDUAL = 2, // L2 absolute residual
  QUDA_HEAVY_QUARK_RESIDUAL = 4, // Fermilab heavy quark residual
  QUDA_INVALID_RESIDUAL = QUDA_INVALID_ENUM
} QudaResidualType;


    / **Parameters for deflated solvers* /
    / ** The precision of the Ritz vectors * /
    QudaPrecision cuda_prec_ritz;
    / ** How many vectors to compute after one solve
     *  for eigCG recommended values 8 or 16
    * /
    int n_ev;
    / ** EeigCG  : Search space dimension
     *  gmresdr : Krylov subspace dimension
     * /
    int max_search_dim;
    / ** For systems with many RHS: current RHS index * /
    int rhs_idx;
    / ** Specifies deflation space volume: total number of eigenvectors is n_ev*deflation_grid * /
    int deflation_grid;
    / ** eigCG: selection criterion for the reduced eigenvector set * /
    double eigenval_tol;
    / ** mixed precision eigCG tuning parameter:  minimum search vector space restarts * /
    int eigcg_max_restarts;
    / ** initCG tuning parameter:  maximum restarts * /
    int max_restart_num;
    / ** initCG tuning parameter:  tolerance for cg refinement corrections in the deflation stage * /
    double inc_tol;

    / ** Whether to make the solution vector(s) after the solve * /
    int make_resident_solution;

    / ** Whether to use the resident solution vector(s) * /
    int use_resident_solution;

    / ** Whether to use the solution vector to augment the chronological basis * /
    int chrono_make_resident;

    / ** Whether the solution should replace the last entry in the chronology * /
    int chrono_replace_last;

    / ** Whether to use the resident chronological basis * /
    int chrono_use_resident;

    / ** The maximum length of the chronological history to store * /
    int chrono_max_dim;

    / ** The index to indicate which chrono history we are augmenting * /
    int chrono_index;

    / ** Precision to store the chronological basis in * /
    QudaPrecision chrono_precision;

    / ** Which external library to use in the linear solvers (Eigen) * /
    QudaExtLibType extlib_type;

typedef enum QudaExtLibType_s {
  QUDA_CUSOLVE_EXTLIB,
  QUDA_EIGEN_EXTLIB,
  QUDA_EXTLIB_INVALID = QUDA_INVALID_ENUM
} QudaExtLibType;


    / ** Whether to use the platform native or generic BLAS / LAPACK * /
    QudaBoolean native_blas_lapack;

    / ** Whether to use fused kernels for mobius * /
    QudaBoolean use_mobius_fused_kernel;

    / **
     * Parameters for distance preconditioning algorithm proposed in arXiv:1006.4028,
     * which is useful to solve a precise heavy quark propagator.
     * alpha0 and t0 follow Eq.(9) in the article.
     * /

    / ** The alpha0 parameter for distance preconditioning, related to the pseudoscalar meson mass * /
    double distance_pc_alpha0;
    / ** The t0 parameter for distance preconditioning, the timeslice where the source is located * /
    int distance_pc_t0;









 invParam.struct_size = sizeof(invParam);

    // create new structure to hold the auxiliary data, then
    // assign information to the QudaMultigridParam structure
    // needed for the preconditioner set up
 QudaMGInfoPtr.reset(new QudaMGInfo);
 QudaInvertParam& mg_inv_param(QudaMGInfoPtr->mg_inv_param);  // to shorten notation
 QudaMultigridParam& mg_param(QudaMGInfoPtr->mg_param);       // to shorten notation
 mg_param.invert_param=&mg_inv_param;

 qactioninfo.setQudaInvertParam(mg_inv_param);
 mg_inv_param.cpu_prec = QudaInfo::get_cpu_prec();
 mg_inv_param.cuda_prec = QudaInfo::get_cuda_prec();
 mg_inv_param.solution_type = QUDA_MAT_SOLUTION;
 mg_inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
 mg_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
 mg_inv_param.tune = QUDA_TUNE_YES;
 mg_inv_param.inv_type = QUDA_GCR_INVERTER;    // the inner solver
 mg_inv_param.tol = invParam.tol;                  // CHECK??
 mg_inv_param.reliable_delta = 1e-5;  
 mg_inv_param.maxiter = invParam.maxiter;     // CHECK??
 mg_inv_param.gcrNkrylov =  xmlputQLInt("PreCondGCRNKrylov",ivalues,ivalindex); //   CHECK???
 mg_inv_param.pipeline = 0;
 mg_inv_param.dagger = QUDA_DAG_NO;
 mg_inv_param.verbosity = getVerbosity();
 mg_inv_param.compute_true_res=true;
 mg_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
 mg_inv_param.cuda_prec_sloppy = QudaInfo::get_cuda_prec_sloppy();
 mg_inv_param.cuda_prec_refinement_sloppy = QudaInfo::get_cuda_prec_sloppy();
 mg_inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
 mg_inv_param.cuda_prec_eigensolver = QudaInfo::get_cuda_prec();
 mg_inv_param.struct_size = sizeof(mg_inv_param);

 mg_param.n_level=xmlputQLInt("NumLevels",ivalues,ivalindex);
 mg_param.setup_type = QUDA_NULL_VECTOR_SETUP; 

 for (int idnum=0;idnum<ivalues[1];++idnum){
    if (xmlputQLInt("IdNum",ivalues,ivalindex)!=idnum){
       errorLaph("Problem occurred in set multigrid parameters");}
    mg_param.coarse_solver[idnum]=static_cast<QudaInverterType>(xmlputQLEnum("CoarseSolverType",
                 {QUDA_CG_INVERTER,QUDA_BICGSTAB_INVERTER,QUDA_GCR_INVERTER,
                  QUDA_CA_GCR_INVERTER,QUDA_CA_CG_INVERTER},ivalues,ivalindex));
    mg_param.coarse_grid_solution_type[idnum]=static_cast<QudaSolutionType>(xmlputQLEnum("CoarseSolutionType",
                 {QUDA_MAT_SOLUTION,QUDA_MATPC_SOLUTION},ivalues,ivalindex));
    mg_param.coarse_solver_tol[idnum]=xmlputQLReal("CoarseSolverTolerance",rvalues,rvalindex);
    mg_param.coarse_solver_maxiter[idnum]=xmlputQLInt("CoarseSolverMaxIterations",ivalues,ivalindex);
    vector<int> blocksizes=xmlputQLIntVector("Blocking",ivalues,ivalindex,LayoutInfo::Ndim);
    for (int dir=0;dir<LayoutInfo::Ndim;++dir){
       mg_param.geo_block_size[idnum][dir]=blocksizes[dir];}
    mg_param.spin_block_size[idnum]=xmlputQLInt("SpinBlockSize",ivalues,ivalindex);
    mg_param.n_vec[idnum]=xmlputQLInt("NumNullVectors",ivalues,ivalindex);
    mg_param.smoother[idnum]=static_cast<QudaInverterType>(xmlputQLEnum("SmootherType",
                 {QUDA_CA_GCR_INVERTER,QUDA_CA_CG_INVERTER},ivalues,ivalindex));
    mg_param.smoother_solve_type[idnum]=static_cast<QudaSolveType>(xmlputQLEnum("SmootherSolveType",
                 {QUDA_DIRECT_PC_SOLVE,QUDA_NORMOP_PC_SOLVE},ivalues,ivalindex));
    mg_param.smoother_tol[idnum]=xmlputQLReal("SmootherTolerance",rvalues,rvalindex);
    mg_param.nu_pre[idnum]=xmlputQLInt("NumPreSmooth",ivalues,ivalindex);
    mg_param.nu_post[idnum]=xmlputQLInt("NumPostSmooth",ivalues,ivalindex);
    mg_param.omega[idnum]=xmlputQLReal("OverRelaxOmega",rvalues,rvalindex);
    mg_param.cycle_type[idnum]=static_cast<QudaMultigridCycleType>(xmlputQLEnum("CycleType",
                 {QUDA_MG_CYCLE_RECURSIVE,QUDA_MG_CYCLE_VCYCLE,QUDA_MG_CYCLE_FCYCLE,
                  QUDA_MG_CYCLE_WCYCLE},ivalues,ivalindex));
    mg_param.smoother_schwarz_cycle[idnum]=xmlputQLInt("SmootherSchwarzCycle",ivalues,ivalindex);   

    mg_param.location[idnum]=QUDA_CUDA_FIELD_LOCATION;
    mg_param.smoother_schwarz_type[idnum]=QUDA_INVALID_SCHWARZ;

    mg_param.precision_null[idnum]=QUDA_HALF_PRECISION;
    mg_param.verbosity[idnum]=QUDA_SILENT;
    mg_param.global_reduction[idnum]=QUDA_BOOLEAN_YES;   // unless using schwarz smoother

    }

 mg_param.compute_null_vector=xmlputQLBool("ComputeNullVectors",ivalues,ivalindex)?QUDA_COMPUTE_NULL_VECTOR_YES:QUDA_COMPUTE_NULL_VECTOR_NO;
 mg_param.generate_all_levels=xmlputQLBool("GenerateAllLevels",ivalues,ivalindex)?QUDA_BOOLEAN_TRUE:QUDA_BOOLEAN_FALSE;
 mg_param.run_verify=xmlputQLBool("RunVerify",ivalues,ivalindex)?QUDA_BOOLEAN_TRUE:QUDA_BOOLEAN_FALSE;

// QudaInvertParam mg_inv_param = newQudaInvertParam();
// const int QUDA_MAX_MG_LEVEL=4;
// QudaEigParam mg_eig_param[QUDA_MAX_MG_LEVEL];
// QudaEigParam eig_param = newQudaEigParam();
// bool use_split_grid = false;



/ *

    / ** Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size * /
    size_t struct_size;

    QudaInvertParam *invert_param;



    / ** The Gflops rate of the multigrid solver setup * /
    double gflops;

    / **< The time taken by the multigrid solver setup * /
    double secs;



  } QudaMultigridParam;
* /

/ *
  if (inv_multigrid) {
    setQudaMgSolveTypes();
    setMultigridInvertParam(inv_param);
    // Set sub structures
    mg_param.invert_param = &mg_inv_param;
    for (int i = 0; i < mg_levels; i++) {
      if (mg_eig[i]) {
        mg_eig_param[i] = newQudaEigParam();
        setMultigridEigParam(mg_eig_param[i], i);
        mg_param.eig_param[i] = &mg_eig_param[i];
      } else {
        mg_param.eig_param[i] = nullptr;
      }
    }
    // Set MG
    setMultigridParam(mg_param);
  } else {
    setInvertParam(inv_param);
  }

  if (inv_deflate) {
    setEigParam(eig_param);
    inv_param.eig_param = &eig_param;
  } else {
    inv_param.eig_param = nullptr;
  }
* /





























/ *


void setQudaDefaultMgTestParams()
{
  // We give here some default values
  for (int i = 0; i < QUDA_MAX_MG_LEVEL; i++) {
    mg_verbosity[i] = QUDA_SUMMARIZE;
#ifdef QUDA_MMA_AVAILABLE
    mg_setup_use_mma[i] = true;
#else
    mg_setup_use_mma[i] = false;
#endif
    mg_dslash_use_mma[i] = false;
    setup_inv[i] = QUDA_BICGSTAB_INVERTER;
    num_setup_iter[i] = 1;
    setup_tol[i] = 5e-6;
    setup_maxiter[i] = 500;
    setup_maxiter_refresh[i] = 20;
    mu_factor[i] = 1.;
    coarse_solve_type[i] = QUDA_INVALID_SOLVE;
    smoother_solve_type[i] = QUDA_INVALID_SOLVE;
    mg_schwarz_type[i] = QUDA_INVALID_SCHWARZ;
    mg_schwarz_cycle[i] = 1;
    smoother_type[i] = QUDA_MR_INVERTER;
    smoother_tol[i] = 0.25;
    coarse_solver[i] = QUDA_GCR_INVERTER;
    coarse_solver_tol[i] = 0.25;
    coarse_solver_maxiter[i] = 100;
    solver_location[i] = QUDA_CUDA_FIELD_LOCATION;
    setup_location[i] = QUDA_CUDA_FIELD_LOCATION;
    nu_pre[i] = 2;
    nu_post[i] = 2;
    n_block_ortho[i] = 1;
    block_ortho_two_pass[i] = true;

    // Default eigensolver params
    mg_eig[i] = false;
    mg_eig_tol[i] = 1e-3;
    mg_eig_n_ev[i] = nvec[i];
    mg_eig_n_kr[i] = 3 * nvec[i];
    mg_eig_require_convergence[i] = QUDA_BOOLEAN_TRUE;
    mg_eig_type[i] = QUDA_EIG_TR_LANCZOS;
    mg_eig_spectrum[i] = QUDA_SPECTRUM_SR_EIG;
    mg_eig_check_interval[i] = 5;
    mg_eig_max_restarts[i] = 100;
    mg_eig_max_ortho_attempts[i] = 10;
    mg_eig_use_normop[i] = QUDA_BOOLEAN_FALSE;
    mg_eig_use_dagger[i] = QUDA_BOOLEAN_FALSE;
    mg_eig_use_poly_acc[i] = QUDA_BOOLEAN_TRUE;
    mg_eig_poly_deg[i] = 100;
    mg_eig_amin[i] = 1.0;
    mg_eig_amax[i] = -1.0; // use power iterations
    mg_eig_save_prec[i] = QUDA_DOUBLE_PRECISION;

    setup_ca_basis[i] = QUDA_POWER_BASIS;
    setup_ca_basis_size[i] = 4;
    setup_ca_lambda_min[i] = 0.0;
    setup_ca_lambda_max[i] = -1.0; // use power iterations

    coarse_solver_ca_basis[i] = QUDA_POWER_BASIS;
    coarse_solver_ca_basis_size[i] = 4;
    coarse_solver_ca_lambda_min[i] = 0.0;
    coarse_solver_ca_lambda_max[i] = -1.0;

    smoother_solver_ca_basis[i] = QUDA_POWER_BASIS;
    smoother_solver_ca_lambda_min[i] = 0.0;
    smoother_solver_ca_lambda_max[i] = -1.0; // use power iterations
  }
}

  // Set QUDA's internal parameters
  gauge_param = newQudaGaugeParam();
  setWilsonGaugeParam(gauge_param);

  inv_param = newQudaInvertParam();
  mg_param = newQudaMultigridParam();
  mg_inv_param = newQudaInvertParam();
  eig_param = newQudaEigParam();

  if (inv_multigrid) {
    setQudaMgSolveTypes();
    setMultigridInvertParam(inv_param);
    // Set sub structures
    mg_param.invert_param = &mg_inv_param;
    for (int i = 0; i < mg_levels; i++) {
      if (mg_eig[i]) {
        mg_eig_param[i] = newQudaEigParam();
        setMultigridEigParam(mg_eig_param[i], i);
        mg_param.eig_param[i] = &mg_eig_param[i];
      } else {
        mg_param.eig_param[i] = nullptr;
      }
    }
    // Set MG
    setMultigridParam(mg_param);
  } else {
    setInvertParam(inv_param);
  }

  if (inv_deflate) {
    setEigParam(eig_param);
    inv_param.eig_param = &eig_param;
  } else {
    inv_param.eig_param = nullptr;
  }

  // set parameters for the reference Dslash, and prepare fields to be loaded
  if (dslash_type == QUDA_DOMAIN_WALL_DSLASH || dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH
      || dslash_type == QUDA_MOBIUS_DWF_DSLASH || dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH) {
    dw_setDims(gauge_param.X, inv_param.Ls);
  } else {
    setDims(gauge_param.X);
  }

  // Allocate host side memory for the gauge field.
  //----------------------------------------------------------------------------
  gauge_.resize(4 * V * gauge_site_size * host_gauge_data_type_size);
  for (int i = 0; i < 4; i++) gauge[i] = gauge_.data() + i * V * gauge_site_size * host_gauge_data_type_size;
  constructHostGaugeField(gauge.data(), gauge_param, argc, argv);

  // Allocate host side memory for clover terms if needed.
  if (dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
    clover.resize(V * clover_site_size * host_clover_data_type_size);
    clover_inv.resize(V * clover_site_size * host_spinor_data_type_size);
    constructHostCloverField(clover.data(), clover_inv.data(), inv_param);
  }

  if (!enable_testing) {
    // Load the gauge field to the device
    loadGaugeQuda(gauge.data(), &gauge_param);

    if (dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
      // Load the clover terms to the device
      loadCloverQuda(clover.data(), clover_inv.data(), &inv_param);
    }

    // Compute plaquette as a sanity check
    double plaq[3];
    plaqQuda(plaq);
    printfQuda("Computed plaquette is %e (spatial = %e, temporal = %e)\n", plaq[0], plaq[1], plaq[2]);
  }
}

std::vector<std::array<double, 2>> solve(test_t param)
{
  inv_param.cuda_prec = ::testing::get<0>(param);
  inv_param.clover_cuda_prec = ::testing::get<0>(param);
  inv_param.cuda_prec_sloppy = ::testing::get<1>(param);
  inv_param.cuda_prec_refinement_sloppy = ::testing::get<1>(param);
  inv_param.cuda_prec_eigensolver = ::testing::get<1>(param);
  inv_param.clover_cuda_prec_sloppy = ::testing::get<1>(param);
  inv_param.clover_cuda_prec_refinement_sloppy = ::testing::get<1>(param);
  inv_param.clover_cuda_prec_eigensolver = ::testing::get<1>(param);
  inv_param.inv_type = ::testing::get<2>(param);
  inv_param.solution_type = ::testing::get<3>(param);
  inv_param.solve_type = ::testing::get<4>(param);
  multishift = ::testing::get<5>(param);
  inv_param.solution_accumulator_pipeline = ::testing::get<6>(param);

  // schwarz parameters
  auto schwarz_param = ::testing::get<7>(param);
  inv_param.schwarz_type           = ::testing::get<0>(schwarz_param);
  inv_param.inv_type_precondition  = ::testing::get<1>(schwarz_param);
  inv_param.cuda_prec_precondition = ::testing::get<2>(schwarz_param);
  inv_param.clover_cuda_prec_precondition = ::testing::get<2>(schwarz_param);

  inv_param.residual_type = ::testing::get<8>(param);

  // reset lambda_max if we're doing a testing loop to ensure correct lambma_max
  if (enable_testing) inv_param.ca_lambda_max = -1.0;

  logQuda(QUDA_SUMMARIZE, "Solution = %s, Solve = %s, Solver = %s, Precision = %s, Sloppy precision = %s\n",
          get_solution_str(inv_param.solution_type), get_solve_str(inv_param.solve_type),
          get_solver_str(inv_param.inv_type), get_prec_str(inv_param.cuda_prec),
          get_prec_str(inv_param.cuda_prec_sloppy));

  // params corresponds to split grid
  for (int i = 0; i < 4; i++) inv_param.split_grid[i] = grid_partition[i];
  int num_sub_partition = grid_partition[0] * grid_partition[1] * grid_partition[2] * grid_partition[3];
  use_split_grid = num_sub_partition > 1;

  // Now QUDA is initialised and the fields are loaded, we may setup the preconditioner
  void *mg_preconditioner = nullptr;
  if (inv_multigrid) {
    if (use_split_grid) { errorQuda("Split grid does not work with MG yet."); }
    mg_preconditioner = newMultigridQuda(&mg_param);
    inv_param.preconditioner = mg_preconditioner;

    printfQuda("MG Setup Done: %g secs, %g Gflops\n", mg_param.secs, mg_param.gflops / mg_param.secs);
  }

  // Vector construct START
  //-----------------------------------------------------------------------------------
  std::vector<quda::ColorSpinorField> in(Nsrc);
  std::vector<quda::ColorSpinorField> out(Nsrc);
  std::vector<quda::ColorSpinorField> out_multishift(multishift * Nsrc);
  quda::ColorSpinorField check;
  quda::ColorSpinorParam cs_param;
  constructWilsonTestSpinorParam(&cs_param, &inv_param, &gauge_param);
  check = quda::ColorSpinorField(cs_param);
  std::vector<std::vector<void *>> _hp_multi_x(Nsrc, std::vector<void *>(multishift));

  // QUDA host array for internal checks and malloc
  // Vector construct END
  //-----------------------------------------------------------------------------------

  // Shifts
  std::vector<double> shifts(multishift);

  // QUDA invert test BEGIN
  //----------------------------------------------------------------------------
  if (multishift > 1) {
    if (use_split_grid) { errorQuda("Split grid does not work with multishift yet."); }
    inv_param.num_offset = multishift;

    // Consistency check for shifts, tols, tols_hq size if we're setting custom values
    if (multishift_masses.size() != 0)
      errorQuda("Multishift masses are not supported for Wilson-type fermions");
    if (multishift_shifts.size() != 0 && multishift_shifts.size() != static_cast<unsigned long>(multishift))
      errorQuda("Multishift shift count %d does not agree with number of shifts passed in %lu\n", multishift, multishift_shifts.size());
    if (multishift_tols.size() != 0 && multishift_tols.size() != static_cast<unsigned long>(multishift))
      errorQuda("Multishift tolerance count %d does not agree with number of masses passed in %lu\n", multishift, multishift_tols.size());
    if (multishift_tols_hq.size() != 0 && multishift_tols_hq.size() != static_cast<unsigned long>(multishift))
      errorQuda("Multishift hq tolerance count %d does not agree with number of masses passed in %lu\n", multishift, multishift_tols_hq.size());

    // Copy offsets and tolerances into inv_param; copy data pointers
    for (int i = 0; i < multishift; i++) {
      shifts[i] = (multishift_shifts.size() == 0 ? (i * i * 0.01) : multishift_shifts[i]);
      inv_param.offset[i] = shifts[i];
      inv_param.tol_offset[i] = (multishift_tols.size() == 0 ? inv_param.tol : multishift_tols[i]);
      inv_param.tol_hq_offset[i] = (multishift_tols_hq.size() == 0 ? inv_param.tol_hq : multishift_tols_hq[i]);

      // Allocate memory and set pointers
      for (int n = 0; n < Nsrc; n++) {
        out_multishift[n * multishift + i] = quda::ColorSpinorField(cs_param);
        _hp_multi_x[n][i] = out_multishift[n * multishift + i].data();
      }
    }
  }

  std::vector<double> time(Nsrc);
  std::vector<double> gflops(Nsrc);
  std::vector<int> iter(Nsrc);

  quda::RNG rng(check, 1234);

  for (int i = 0; i < Nsrc; i++) {
    // Populate the host spinor with random numbers.
    in[i] = quda::ColorSpinorField(cs_param);
    spinorNoise(in[i], rng, QUDA_NOISE_GAUSS);
    out[i] = quda::ColorSpinorField(cs_param);
  }

  if (distance_pc_alpha0 != 0.0 && distance_pc_t0 >= 0) {
    inv_param.distance_pc_alpha0 = distance_pc_alpha0;
    inv_param.distance_pc_t0 = distance_pc_t0;
    verifySpinorDistanceReweight(in[0], distance_pc_alpha0, distance_pc_t0);
  }

  if (!use_split_grid) {

    for (int i = 0; i < Nsrc; i++) {
      // If deflating, preserve the deflation space between solves
      if (inv_deflate) eig_param.preserve_deflation = i < Nsrc - 1 ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
      // Perform QUDA inversions
      if (multishift > 1) {
        invertMultiShiftQuda(_hp_multi_x[i].data(), in[i].data(), &inv_param);
      } else {
        invertQuda(out[i].data(), in[i].data(), &inv_param);
      }

      time[i] = inv_param.secs;
      gflops[i] = inv_param.gflops / inv_param.secs;
      iter[i] = inv_param.iter;
      printfQuda("Done: %i iter / %g secs = %g Gflops\n", inv_param.iter, inv_param.secs,
                 inv_param.gflops / inv_param.secs);
    }
  } else {

    inv_param.num_src = Nsrc;
    inv_param.num_src_per_sub_partition = Nsrc / num_sub_partition;
    // Host arrays for solutions, sources, and check
    std::vector<void *> _hp_x(Nsrc);
    std::vector<void *> _hp_b(Nsrc);
    for (int i = 0; i < Nsrc; i++) {
      _hp_x[i] = out[i].data();
      _hp_b[i] = in[i].data();
    }
    
		// Run split grid
    invertMultiSrcQuda(_hp_x.data(), _hp_b.data(), &inv_param);

    quda::comm_allreduce_int(inv_param.iter);
    inv_param.iter /= quda::comm_size() / num_sub_partition;
    quda::comm_allreduce_sum(inv_param.gflops);
    inv_param.gflops /= quda::comm_size() / num_sub_partition;
    quda::comm_allreduce_max(inv_param.secs);
    printfQuda("Done: %d sub-partitions - %i iter / %g secs = %g Gflops\n", num_sub_partition, inv_param.iter,
               inv_param.secs, inv_param.gflops / inv_param.secs);
  }

  // QUDA invert test COMPLETE
  //----------------------------------------------------------------------------

  // free the multigrid solver
  if (inv_multigrid) destroyMultigridQuda(mg_preconditioner);

  // Compute performance statistics
  if (Nsrc > 1 && !use_split_grid) performanceStats(time, gflops, iter);

  std::vector<std::array<double, 2>> res(Nsrc);
  // Perform host side verification of inversion if requested
  if (verify_results) {
    for (int i = 0; i < Nsrc; i++) {
      res[i] = verifyInversion(out[i].data(), _hp_multi_x[i].data(), in[i].data(), check.data(), gauge_param, inv_param,
                               gauge.data(), clover.data(), clover_inv.data());
    }
  }
  return res;
}
* /




 mg_param.struct_size = sizeof(mg_param);
*/
}




// *******************************************************************
}
