#ifndef LAPH_DWF_H
#define LAPH_DWF_H

#ifdef LAPH_DOMAIN_WALL
#include "QudaLaphIncludes.h"
namespace LaphEnv {

// code from witnesser sortof
template <class T>
static void
DirtyWilson4D( const size_t fermsize,
	       const QudaInvertParam inv_param,
	       std::complex<T> sol[] )
{
  std::complex<T> *tmp = new std::complex<T>[ fermsize ] ;
  memcpy( tmp , sol , fermsize*sizeof(std::complex<T>) ) ;
  
  QudaInvertParam Wil = newQudaInvertParam() ;
  Wil = inv_param ;
  
  Wil.dslash_type = QUDA_WILSON_DSLASH ;
  Wil.kappa = 1/(2.*(4.+inv_param.m5)) ;
  Wil.mass = inv_param.m5 ;

  void *in  = (void*)sol , *out = (void*)tmp ;
  
  MatQuda( out , in , &Wil ) ;

  const std::complex<T> cfac1 = 0.5*(inv_param.c_5[0])/Wil.kappa ;
  const T cmul = cfac1.real() ;
  
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < fermsize; i++ ) {
    sol[i] = cmul*tmp[i] - sol[i] ;
  }
  delete [] tmp ;
}

// does 0.5*(1+\gamma_5) in the NR basis that we seemingly have to use
template <class T>
static void
chiralProjectPlus( const size_t fermsize,
		   const void *in,
		   std::complex<T> tmp[] )
{
  std::complex<T> *ptin = (std::complex<T>*)in ;
  T mul = (T)0.5 ;
  size_t i ;
#pragma omp parallel for private(i) collapse(2)
  for( i = 0 ; i < fermsize/(12) ; i++ ) {
    for( int c = 0 ; c < 3 ; c++ ) {
      tmp[ c +     12*i ] = mul*( ptin[ c +     12*i ] - ptin[ c + 6 + 12*i ] ) ;
      tmp[ c + 3 + 12*i ] = mul*( ptin[ c + 3 + 12*i ] - ptin[ c + 9 + 12*i ] ) ;
      tmp[ c + 6 + 12*i ] = mul*( ptin[ c + 6 + 12*i ] - ptin[ c + 0 + 12*i ] ) ;
      tmp[ c + 9 + 12*i ] = mul*( ptin[ c + 9 + 12*i ] - ptin[ c + 3 + 12*i ] ) ;
    }
  }
}

// chiminus = 0.5*(1-\gamma_5)
template <class T>
static void
chiralProjectMinus( const size_t fermsize,
                    const void *in,
		    std::complex<T> tmp[] )
{
  std::complex<T> *ptin = (std::complex<T>*)in ;
  T mul = (T)0.5 ;
  size_t i ;
#pragma omp parallel for private(i) collapse(2)
  for( i = 0 ; i < fermsize/(12) ; i++ ) {
    for( int c = 0 ; c < 3 ; c++ ) {
      tmp[ c +     12*i ] = mul*( ptin[ c +     12*i ] + ptin[ c + 6 + 12*i ] ) ;
      tmp[ c + 3 + 12*i ] = mul*( ptin[ c + 3 + 12*i ] + ptin[ c + 9 + 12*i ] ) ;
      tmp[ c + 6 + 12*i ] = mul*( ptin[ c + 6 + 12*i ] + ptin[ c + 0 + 12*i ] ) ;
      tmp[ c + 9 + 12*i ] = mul*( ptin[ c + 9 + 12*i ] + ptin[ c + 3 + 12*i ] ) ;
    }
  }
}

template <class T>
static void
solve_DWF( void *out ,
	   void *in ,
	   QudaInvertParam *inv_param)
{
  const int Ls = inv_param -> Ls ;
  const int fullsize = LayoutInfo::getRankLatticeNumSites() ;
  const size_t fermsize = fullsize*3*4 ;
  const size_t half_fermsize = fermsize/2 ;
  std::complex<T> *spinorIn  = new std::complex<T>[ Ls*fermsize ] ;
  std::complex<T> *spinorOut = new std::complex<T>[ Ls*fermsize ] ;

  memset(reinterpret_cast<char*>(spinorIn) , 0, fermsize*Ls*sizeof(std::complex<T>));
  memset(reinterpret_cast<char*>(spinorOut), 0, fermsize*Ls*sizeof(std::complex<T>));
  
  // chiral project
  std::complex<T> *tmp1 = new std::complex<T>[ fermsize ] ;

  chiralProjectPlus<T>( fermsize, in , tmp1 ) ;
  DirtyWilson4D<T>( fermsize , *inv_param , tmp1 ) ;
  
  std::complex<T> *tmp2 = new std::complex<T>[ fermsize ] ;

  chiralProjectMinus<T>( fermsize, in , tmp2 ) ;
  DirtyWilson4D<T>( fermsize , *inv_param , tmp2 ) ;
  
  // copies for the checkerboarded sources
  memcpy(reinterpret_cast<char*>(&spinorIn[0]),
	 reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp1[0]))),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(&spinorIn[half_fermsize*Ls]),
	 reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp1[half_fermsize]))),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(&spinorIn[half_fermsize*(Ls-1)]),
	 reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp2[0]))),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(&spinorIn[half_fermsize*(2*Ls-1)]),
	 reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp2[half_fermsize]))),
	 half_fermsize*sizeof(std::complex<T>));

  invertQuda( reinterpret_cast<void*>(spinorOut),
	      reinterpret_cast<void*>(spinorIn), inv_param ) ;
  delete [] spinorIn ;
  
  // undo the solve
  memcpy(reinterpret_cast<char*>(const_cast<std::complex<T>*>(&tmp1[0])),
	 reinterpret_cast<char*>(&spinorOut[0]),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp1[half_fermsize]))),
	 reinterpret_cast<char*>(&spinorOut[half_fermsize*Ls]),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp2[0]))),
	 reinterpret_cast<char*>(&spinorOut[half_fermsize*(Ls-1)]),
	 half_fermsize*sizeof(std::complex<T>));
  memcpy(reinterpret_cast<char*>(const_cast<std::complex<T>*>(&(tmp2[half_fermsize]))),
	 reinterpret_cast<char*>(&spinorOut[half_fermsize*(2*Ls-1)]),
	 half_fermsize*sizeof(std::complex<T>));
  delete [] spinorOut ;

  const std::complex<T> invTwoKappaB = inv_param -> b_5[0]*( 4 + inv_param -> m5 ) + 1.0;
  const T twoKappaB = 1/(2*invTwoKappaB.real());
  
  // chiral projections
  std::complex<T> *outptr = (std::complex<T>*)out ;
  size_t i ;
  // out = chiralprojecminus(tmp1) + chiralprojectplus(tmp2)
#pragma omp parallel for private(i) collapse(2)
  for( i = 0 ; i < fermsize/12 ; i++ ) {
    for( int c = 0 ; c < 3 ; c++ ) {
      outptr[ c + 0 + 12*i ] = (twoKappaB)*(   tmp1[ c + 0 + 12*i ] + tmp1[ c + 6 + 12*i ]
					     + tmp2[ c + 0 + 12*i ] - tmp2[ c + 6 + 12*i ] ) ;
      outptr[ c + 3 + 12*i ] = (twoKappaB)*(   tmp1[ c + 3 + 12*i ] + tmp1[ c + 9 + 12*i ]
					     + tmp2[ c + 3 + 12*i ] - tmp2[ c + 9 + 12*i ] ) ;
      outptr[ c + 6 + 12*i ] = (twoKappaB)*(   tmp1[ c + 6 + 12*i ] + tmp1[ c + 0 + 12*i ]
					     + tmp2[ c + 6 + 12*i ] - tmp2[ c + 0 + 12*i ] ) ;
      outptr[ c + 9 + 12*i ] = (twoKappaB)*(   tmp1[ c + 9 + 12*i ] + tmp1[ c + 3 + 12*i ]
					     + tmp2[ c + 9 + 12*i ] - tmp2[ c + 3 + 12*i ] ) ;
    }
  }
  delete [] tmp1 ;
  delete [] tmp2 ;
}
}
#endif  
#endif
