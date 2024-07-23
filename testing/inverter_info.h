#ifndef LAPH_INVERTER_INFO_H
#define LAPH_INVERTER_INFO_H

#include "xml_handler.h"
#include "quark_action_info.h"
#include "quda.h"
#include <memory>

namespace LaphEnv {


/*
  CG  (Conjugate Gradient: use MdagM)
  CGNE (Conjugate Gradient on the Normal Equations)
  CGNR (Conjugate Gradient on the Normal Equations
  BICGSTAB (Biconjugate Gradients, Stabilized)
  BICGSTABL
  GCR (Generalized Conjugate Residual)
  MR  (Minimum Residual: use MdagM)
  GMRESDR (Generalized Minimum Residual with Deflated Restarts)

  CA_CG (Communication Avoiding Conjugate Gradient: use MdagM)
  CA_CGNE (Communication Avoiding Conjugate Gradient on the Normal Equations)
  CA_CGNR
  CA_GCR (Communication Avoiding Generalized Conjugate Residual)


  SD
  PCG
  EIGCG
  INC_EIGCG
  GMRESDR_PROJ
  GMRESDR_SH
  FGMRESDR
  MG
  CG3
  CG3NE
  CG3NR

Power method A simple iteration for approximating the principal eigenvector and eigen-
value of a matrix A, by means of repeated sparse matrix-vector multiplications.


Communication Avoiding
Messages between processors, in a distributed-memory system
Cache coherency traffic, in a shared-memory system
Data transfers between coprocessors linked by a bus, such as between a CPU and a GPU

s-step Krylov methods: goal is to take s steps of a Krylov Space Method for the 
same communication cost as 1 step for standard method.

Two basis choices for communication avoiding:
  QUDA_POWER_BASIS,
  QUDA_CHEBYSHEV_BASIS,



CG

Method of Conjugate Gradients for solving symmetric positive definite linear systems
from M. R. Hestenes and E. Stiefel, Methods of conjugate gradients for solving linear
systems, Journal of Research of the National Bureau of Standards, 49 (1952).


CGNE: M Mdag y = b is solved; x = Mdag y is returned as solution.
CGNR: Mdag M x = Mdag b is solved.

GMRESDR
The Generalized Minimum Residual method of Saad and Schultz [209], a Krylov
subspace method for solving nonsymmetric linear systems.



*/




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
// **********************************************************************


class InverterInfo
{

    std::vector<std::string> svalues;
    
    std::vector<int> ivalues;
    
    std::vector<double> rvalues;


  public:
  
    InverterInfo(const XMLHandler& xmlin);

    InverterInfo(const InverterInfo& rhs);

    InverterInfo& operator=(const InverterInfo& rhs);

    ~InverterInfo();
    
    std::string getName() const {return svalues[0];}
    
    uint getMaxIterations() const {return ivalues[0];}

    void checkEqual(const InverterInfo& rhs) const;

    bool operator==(const InverterInfo& rhs) const;

    std::string output(int indent = 0) const;

    void output(XMLHandler& xmlout) const;

    void setQudaInvertParam(QudaInvertParam& invParam, 
                            const QuarkActionInfo& qactioninfo) const;


  private:

    void set_info_cgnr(XMLHandler& xmlr);

    void output_cgnr(XMLHandler& xmlout) const;
 
    void setQudaInvertParam_cgnr(QudaInvertParam& invParam) const;

    void set_info_bicgstab(XMLHandler& xmlr);

    void output_bicgstab(XMLHandler& xmlout) const;
 
    void setQudaInvertParam_bicgstab(QudaInvertParam& invParam) const;

    void set_info_gcr(XMLHandler& xmlr);

    void output_gcr(XMLHandler& xmlout) const;
 
    void setQudaInvertParam_gcr(QudaInvertParam& invParam) const;

    void set_info_gcr_multigrid(XMLHandler& xmlr);

    void output_gcr_multigrid(XMLHandler& xmlout) const;
 
    void setQudaInvertParam_gcr_multigrid(QudaInvertParam& invParam,
                                          const QuarkActionInfo& qactioninfo) const;

       // auxiliary structures needed as part of an invParam
       // for the complicated multigrid inverter

    struct QudaMGInfo {

       QudaInvertParam mg_inv_param;

       QudaMultigridParam mg_param;

       std::vector<QudaEigParam> mg_eig_param;

       QudaMGInfo()
       {mg_inv_param = newQudaInvertParam();
        mg_param = newQudaMultigridParam();}

       ~QudaMGInfo(){}

    };

    mutable std::unique_ptr<QudaMGInfo> QudaMGInfoPtr;
    
    friend class QuarkHandler;

};


// *****************************************************************
}
#endif
