#ifndef LAPH_INVERTER_INFO_H
#define LAPH_INVERTER_INFO_H

#include "quark_action_info.h"
#include "quda.h"
#include "xml_handler.h"
#include <memory>

namespace LaphEnv {

// **********************************************************************
// *                                                                    *
// *  Class "InverterInfo" holds information about the inverter.        *
// *  The role of the inverters is to solve the linear system of        *
// *  equations   M*x=y  for given source y.  Conjugate gradient        *
// *  methods actually solve M^dagger M x = M^dagger y, but at the      *
// *  cost of a worse condition number.  The solution methods are       *
// *  iterative, so one must specify the maximum number of iterations   *
// *  before giving up and the relative residuum for adequate           *
// *  convergence to a solution.  Convergence is usually specified      *
// *  by   r = Mx-y  with  norm(r) <= eps^2 * norm(y).  This            *
// *  convergence criterion is applied to the **pre-conditioned**       *
// *  system (pre-conditioning is usually done, based on the            *
// *  fermion action specified).  The top-level solver reports the      *
// *  obtained residual (the relative residual is divided by norm(y)).  *
// *  Note that the relative residuum for the original                  *
// *  not-preconditioned system may be a little larger (usually         *
// *  a factor of 2-3) than for the preconditioned system.              *
// *                                                                    *
// *  An important member is                                            *
// *                                                                    *
// *    setQudaInvertParam(QudaInvertParam& invParam,                   *
// *                       const QuarkActionInfo& qactioninfo) const;   *
// *                                                                    *
// *  This sets up the struct "invParam" which quda needs to perform    *
// *  the inversion.                                                    *
// *                                                                    *
// *  Sample XML inputs are given below. All tags except <Name> are     *
// *  optional, and the default values are given.                       *
// *                                                                    *
// *    Conjugate-Gradient on the Normal Equations:                     *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>CGNR</Name>                                              *
// *     <Tolerance>1.0e-10</Tolerance>                                 *
// *     <MaxIterations>10000</MaxIterations>                           *
// *   </InvertInfo>                                                    *
// *                                                                    *
// *    Biconjugate Gradient Stabilized:                                *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>BICGSTAB</Name>                                          *
// *     <Tolerance>1.0e-10</Tolerance>                                 *
// *     <MaxIterations>6000</MaxIterations>                            *
// *   </InvertInfo>                                                    *
// *                                                                    *
// *    Generalized Conjugate Residual (GCR)                            *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>GCR</Name>                                               *
// *     <Tolerance>1.0e-10</Tolerance>                                 *
// *     <MaxIterations>5000</MaxIterations>                            *
// *     <NKrylov>16</NKrylov>                                          *
// *   </InvertInfo>                                                    *
// *                                                                    *
// *    Generalized Conjugate Residual (GCR) with Adaptive Multigrid    *
// *      preconditioning                                               *
// *                                                                    *
// *   <InvertInfo>                                                     *
// *     <Name>GCR_MULTIGRID</Name>                                     *
// *     <Tolerance>1.0e-12</Tolerance>                                 *
// *     <MaxIterations>5000</MaxIterations>                            *
// *     <NKrylov>24</NKrylov>                                          *
// *     .... (see below)                                               *
// *   </InvertInfo>                                                    *
// *                                                                    *
// **********************************************************************

// *******************************************************************************
// * *
// *    GCR Inverter with Adaptive Multigrid Preconditioning *
// * *
// *      Limitations: - only QUDA_MG_CYCLE_RECURSIVE cycle type allowed *
// *                   - deflation only allowed on coarsest grid *
// *                   - inverse iteration used for null vector generation *
// *                   - null space set up solver is CGNR *
// *                   - coarse solvers: GCR except CA_GCR coarsest level *
// *                   - smoother solvers: CA_GCR *
// *                   - any use of communications avoiding inverters *
// *                       uses Chebyshev basis with basis size 4 *
// * *
// *   Preconditioning:  The goal is to solve M x = b for x, where M is *
// *   the Dirac clover matrix.  The idea behind preconditioning is to *
// *   find a matrix K that is an approximation to M^(-1) but is easy to *
// *   apply onto an arbitrary vector.  Then one solves K M x = K b. *
// *   The matrix K M should be close to the identity matrix, so its *
// *   condition number is much smaller than that of M, and an iterative *
// *   solver convergences much faster. *
// * *
// *   Adaptive Multigrid (MG) Preconditioning is a very complicated method *
// *   of forming a preconditioning matrix.  This method is most helpful *
// *   for large lattices and very light quark masses since it combats *
// *   critical slowing down very well.  Note that MG forms a *
// *   preconditioning matrix for EACH iteration of the outer inverter *
// *   (GCR here). So the preconditioning matrix can change iteration by *
// *   iteration, making the construction adaptive. An outer solver *
// *   that can use variable preconditioning (flexible) must be used. *
// * *
// *   MG preconditioning is characterized by forming the preconditioner *
// *   matrix in a sequence of constructions on coarser and coarser *
// *   grids.  Construction "transfers" to coarser lattices (restriction) *
// *   and back to the finer lattices (prolongation).  The use of near *
// *   null space vectors to define the transfer operations in an *
// *   adaptive procedure is a key feature, and "Nlevels" is a key *
// *   parameter.  Solution "smoothing" before and after coarsening *
// *   is also an important feature. *
// * *
// *   Set up step: *
// * *
// *      Construct the prolongation matrix P[l] on level l=0,1,...,Nlevels-2. *
// *      The restriction matrix is R[l]=P[l]^dagger.  The coarse grid matrix *
// *      is given by K[l+1] = R[l] K[l] P[l]  so P[l] is an n[l] X n[l+1] *
// *      matrix where n[l+1] < n[l]. *
// * *
// *      In the adaptive algorithm, one lets the MG process itself define the *
// *      appropriate prolongator by an iterative procedure. A set *
// *      of near null vectors on each level except the last is generated. *
// *      The method employed (here) is inverse iteration.  For several *
// *      random vectors eta, the homogeneous system M eta = 0 is solved *
// *      for some number of iterations. CGNR is used here.  Iterative *
// *      solvers converge most slowly for the low modes, so the resulting *
// *      vectors should mostly contain low modes.  A block orthonormalization *
// *      of the vectors and a lattice blocking scheme is used to define the *
// *      prolongation and restriction operators.  Parameters (for each level *
// *      except last) important for this step are *
// *           <NullSpaceDim>S</NullSpaceDim>  (S=small 24, L=large 32) *
// *           <NullSpaceSolveSteps>600<NullSpaceSolveSteps> *
// * *
// *   MG cycle: *
// * *
// *      The outer solver starts with an initial guess vector, then iterates *
// *      modifying the vector until convergence to a small residual achieves. *
// *      At the beginning of each iteration of the outer solver, the MG *
// *      preconditioner is called to form a new preconditioning matrix, *
// *      which it does in a recursive cycle of MG steps, transferring to a *
// *      coarser grid and calling itself. *
// * *
// *   MG step: *
// * *
// *      An MG step in level "L" is the following sequence of steps: *
// *         -- presmoothing some number of times using an inverter *
// *         -- restrict to coarse grid, call MG step for level "L+1" *
// *              to correct *
// *         -- prolongate back to this grid *
// *         -- postsmoothing some number of types using an inverter *
// *      The coarse grid solve corrects directions that the smoother *
// *      struggles to fix. *
// * *
// *   Coarsest level: *
// * *
// *      Only a solve is done at the coarsest level, and low-mode *
// *      deflation is applied here. *
// * *
// *   XML input should have the following form: *
// *      (all tags except initial <Name> optional; default values shown, *
// *       except as discussed below) *
// * *
// *   <InvertInfo> *
// *     <Name>GCR_MULTIGRID</Name> *
// *     <Tolerance>1.0e-11</Tolerance> *
// *     <MaxIterations>200</MaxIterations> *
// *     <NKrylov>24</NKrylov> *
// *     <MGPreconditioner> *
// *        <NumLevels>2</NumLevels>    2,3,4 (2 or 3 usually best) *
// *        <Level0> *
// *           <XYZTBlockExtents>4 4 4 4</XYZTBlockExtents> (see below) *
// *           <NullSpaceDim>S</NullSpaceDim>  (S=small 24, L=large 32) *
// *           <NullSpaceSolveSteps>500</NullSpaceSolveSteps> *
// *           <CoarseSolverTolerance>0.25</CoarseSolverTolerance> *
// *           <CoarseSolverMaxIterations>50</CoarseSolverMaxIterations> *
// *           <NumPreSmooth>0</NumPreSmooth> *
// *           <NumPostSmooth>4</NumPostSmooth> *
// *        </Level0> *
// *        <Level1> ...similar to level 0 </Level1> *
// *        <Level2>  (the coarsest level) *
// *           <CoarseSolverTolerance>0.25</CoarseSolverTolerance> *
// *           <CoarseSolverMaxIterations>50</CoarseSolverMaxIterations> *
// *           <TRLMDeflation>...</TRLMDeflation>  (default: absent) *
// *        </Level2> *
// *        <LoadNullVectors>false</LoadNullVectors>     compute or reload *
// *        <GenerateAllLevels>true</GenerateAllLevels>  all or level 0 only *
// *        <PreOrthoNullVectors>true</PreOrthoNullVectors> *
// *        <PostOrthoNullVectors>true</PostOrthoNullVectors> *
// *        <SetupMinimizeMemory>false<SetupMinimizeMemory> *
// *        <RunVerify>false</RunVerify>   (use true for initial runs) *
// *     </MGPreconditioner> *
// *   </InvertInfo> *
// * *
// *  Input XML for the deflation on coarsest level (if present) *
// * *
// *   <TRLMDeflation> *
// *      <NumVectors> 128 </NumVectors> *
// *      <Tolerance> 1e-6 </Tolerance> *
// *      <MaxIterations> 200 </MaxIterations> *
// *      <ChebyshevOrder> 8 </ChebyshevOrder> *
// *      <CutoffEigenvalue> 0.1 </CutoffEigenvalue> *
// *   </TRLMDeflation> *
// * *
// *  On an L^3 x T lattice: *
// *     - default value for <NumLevels> is 2 *
// *     - default <XYZTBlockExtents> values determined based on L,T *
// * *
// *  Tips: *
// *     - block extents are the size of each block; for example, if there are *
// *       24 sites in a direction and the block extent is 3, the coarse *
// *       lattice will have 8 sites in that direction, each site including *
// *       3 sites of the original lattice in that direction *
// *     - block extents on a level must be divisors of fine or coarse local *
// *       (each mpi rank) lattice of previous level *
// *     - block extents should be as large as possible without the product of *
// *       the four block extents exceeding quda limit 1024 *
// *     - NullSpaceSolveSteps should be around 500; want the vectors created *
// *       to have significant overlap with null space; solver tolerance *
// *       should be low to force solver to reach number of steps requested *
// *     - CoarseSolverTolerance should be kept large, around 0.25 since *
// *       only approximate solutions needed for preconditioning matrix *
// *     - NumPreSmooth can be set to 0, NumPostSmooth is important and should *
// *       be around 4 to 8. *
// *     - Deflation on the coarsest level will almost always improve solution *
// *       time.  The number of vectors should be reasonably large, such as *
// *       100 to 200. *
// * *
// *  Typical XML input is below. Deflation is not turned on by default, so *
// *  your XML input should turn it on.  You might wish to control the blocking
// *
// *  yourself. *
// * *
// *   <InvertInfo> *
// *     <Name>GCR_MULTIGRID</Name> *
// *     <Tolerance>1.0e-11</Tolerance> *
// *     <MaxIterations>200</MaxIterations> *
// *     <NKrylov>24</NKrylov> *
// *     <MGPreconditioner> *
// *        <NumLevels>2</NumLevels> *
// *        <Level0> *
// *           <XYZTBlockExtents>4 4 4 4</XYZTBlockExtents> *
// *        </Level0> *
// *        <Level1> *
// *           <TRLMDeflation/> *
// *        </Level1> *
// *     </MGPreconditioner> *
// *   </InvertInfo> *
// * *
// *******************************************************************************

class InverterInfo {

  std::vector<std::string> svalues;

  std::vector<int> ivalues;

  std::vector<double> rvalues;

public:
  InverterInfo(const XMLHandler &xmlin);

  InverterInfo(const InverterInfo &rhs);

  InverterInfo &operator=(const InverterInfo &rhs);

  ~InverterInfo();

  std::string getName() const { return svalues[0]; }

  uint getMaxIterations() const { return ivalues[0]; }

  void checkEqual(const InverterInfo &rhs) const;

  bool operator==(const InverterInfo &rhs) const;

  std::string output(int indent = 0) const;

  void output(XMLHandler &xmlout) const;

  void setQudaInvertParam(QudaInvertParam &invParam,
                          const QuarkActionInfo &qactioninfo) const;

private:
  void set_info_cgnr(XMLHandler &xmlr);

  void output_cgnr(XMLHandler &xmlout) const;

  void setQudaInvertParam_cgnr(QudaInvertParam &invParam) const;

  void set_info_bicgstab(XMLHandler &xmlr);

  void output_bicgstab(XMLHandler &xmlout) const;

  void setQudaInvertParam_bicgstab(QudaInvertParam &invParam) const;

  void set_info_gcr(XMLHandler &xmlr);

  void output_gcr(XMLHandler &xmlout) const;

  void setQudaInvertParam_gcr(QudaInvertParam &invParam) const;

  void set_info_gcr_multigrid(XMLHandler &xmlr);

  void output_gcr_multigrid(XMLHandler &xmlout) const;

  void
  setQudaInvertParam_gcr_multigrid(QudaInvertParam &invParam,
                                   const QuarkActionInfo &qactioninfo) const;

  void initDeflationParam_gcr_multigrid() const;

  // auxiliary structures needed as part of an invParam
  // for the complicated multigrid inverter

  struct QudaMGInfo {

    QudaInvertParam mg_inv_param;

    QudaMultigridParam mg_param;

    QudaEigParam coarse_deflate_param;

    QudaMGInfo() {
      mg_inv_param = newQudaInvertParam();
      mg_param = newQudaMultigridParam();
      coarse_deflate_param = newQudaEigParam();
    }

    ~QudaMGInfo() {}
  };

  mutable std::unique_ptr<QudaMGInfo> QudaMGInfoPtr;

  void factorize(uint ivalue, std::list<uint> &factors) const;

  void do_blocking(const std::vector<int> &extents,
                   std::vector<int> &block_extents) const;

  bool check_blocking(std::vector<int> &block_extents,
                      std::vector<int> &extents, std::string &errmessage) const;

  const uint mg_blockprodmin = 4;
  const uint mg_blockprodmax = 1024;

  friend class QuarkHandler;
  friend class PerambulatorHandler;
};

// *****************************************************************
} // namespace LaphEnv
#endif
