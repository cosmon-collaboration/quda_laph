#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include "latt_field.h"
#include "gauge_configuration_info.h"
#include "quark_action_info.h"
#include "momenta.h"
#include <complex>


namespace LaphEnv {


// **************************************************************************
// *                                                                        *
// *   This file contains routines for performing various lattice field     *
// *   operations on the CPU.  Typically, these routines are used for       *
// *   performing checks on solutions obtained by QUDA.  These routines     *
// *   are not particularly slow, but they are not particularly fast.       *
// *   They should only by used for occasional checks.                      *
// *                                                                        *
// **************************************************************************

    //  This routine applies an inner product 
    //         conj(leftfield).rightfield     if "left_conj" is true,
    //              leftfield .rightfield     if "left_conj" is false,
    //  where leftfield and rightfield are fields of the same site type, but the
    //  inner product is taken over the time slices.  The ntime inner
    //  products are returned.  

std::vector<std::complex<double>> getTimeSlicedInnerProducts(const LattField& leftfield, 
                                                             const LattField& rightfield, 
                                                             bool left_conj=true);


void setConstantField(LattField& field, const std::complex<double>& zconst);

void setUnitField(LattField& field);

void setZeroField(LattField& field);

   // Sets the zeroth element at a site to exp(I*sum(site[k]*momfactors[k],k))
   // and each subsequent element at the site gets multiplied by another
   // factor of exp(I*local_phase)

void setVariablePhaseField(LattField& field, const std::vector<double>& momfactors,
                           const double& local_phase);

void compare_latt_fields(const LattField& src1, const LattField& src2);

void convertSpinBasis(LattField& outferm, const LattField& inferm, const std::string& from_to);

   //   Applies the clover Dirac operation to "infield", returning
   //   the result in "outfield".  This operation is
   //
   //    outfield =  [ 1/(2*kappa) - (1/2) Dterm  + CFterm ] infield
   //
   //   where
   //
   //        Dterm = sum_mu [ (1-gamma_mu) U  + (1+gamma_mu) * U^dag ]
   //
   //        CFterm = csw (i/4)  sigma[mu,nu] F[mu,nu]
   //
   //            sigma[mu,nu] (i/2) [gamma_mu, gamma_nu]
   //
   //            F[mu,nu] = (1/8) ( Q[mu,nu]-Q[nu,mu] )
   //
   //            Q[mu,nu] = U(mu,nu,-mu,-nu) + U(nu,-mu,-nu,mu)
   //                     + U(-mu,-nu,mu,nu) + U(-nu,mu,nu,-mu)
   //
   //   Lattice shifts must take the fermion temporal boundary
   //   conditions into account.  This routine assumes the Dirac-Pauli
   //   convention for the Dirac gamma-matrices.
   

void applyCloverDirac(LattField& outfield, const LattField& infield,
                      const std::vector<LattField>& gauge_field,
                      const GaugeConfigurationInfo& gaction,
                      const QuarkActionInfo& qaction);

    //  Applies the 3d spatial Laplacian with a smeared gauge field
    //  onto "infield", returning result in "outfield".  These fields
    //  must be color vectors

void applyMinusSpatialLaplacian(LattField& outfield, const LattField& infield,
                                const std::vector<LattField>& smeared_gauge_field);



    //  Evaluates sum_a conj( qbar[a](x,t) ) q[a](x,t)   a = color index

void doColorContract(LattField& result, const LattField& qbar, const LattField& q);

    //  Evaluates   sum_a,b,c epsilson(a,b,c)  q1[a](x,t) q2[b](x,t) q3[c](x,t) 

void doColorContract(LattField& result, const LattField& q1, const LattField& q2, 
                     const LattField& q3);

    //  Evaluates   sum_a,b epsilson(a,b,c)  q1[a](x,t) q2[b](x,t)

void doColorCrossProduct(LattField& result, const LattField& q1, const LattField& q2);

    //  Evaluates sum_a  q1[a](x,t) q2[a](x,t)   a = color index

void doColorVectorContract(LattField& result, const LattField& q1, const LattField& q2);

    //  Returns the field   exp( -I * pvec.x )  in  "phases"

void makeMomentumPhaseField(LattField& phases, const Momentum& pvec);


// **************************************************************************
}
#endif
