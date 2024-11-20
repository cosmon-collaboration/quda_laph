#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include "gauge_configuration_info.h"
#include "latt_field.h"
#include "quark_action_info.h"
#include <complex>
#include <vector>

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

//  This routine applies an inner product conj(leftfield).rightfield,
//  where leftfield and rightfield are color-vector fields, but the
//  inner product is taken over the time slices.  The ntime inner
//  products are returned.

std::vector<std::complex<double>>
getTimeSlicedInnerProducts(const LattField &leftfield,
                           const LattField &rightfield);

void setConstantField(LattField &field, const std::complex<double> &zconst);

void setUnitField(LattField &field);

void setZeroField(LattField &field);

void compare_latt_fields(const LattField &src1, const LattField &src2);

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
//   conditions into account.

void applyCloverDirac(LattField &outfield, const LattField &infield,
                      const std::vector<LattField> &gauge_field,
                      const GaugeConfigurationInfo &gaction,
                      const QuarkActionInfo &qaction);

//  Applies the 3d spatial Laplacian with a smeared gauge field
//  onto "infield", returning result in "outfield".  These fields
//  must be color vectors

void applyMinusSpatialLaplacian(
    LattField &outfield, const LattField &infield,
    const std::vector<LattField> &smeared_gauge_field);

} // namespace LaphEnv
#endif
