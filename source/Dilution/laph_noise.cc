#include "laph_noise.h"
#include "laph_stdio.h"
#include <cstdint>

namespace LaphEnv {

UniformDeviate32::UniformDeviate32() {
  Reseed(0); // choose seed based on clock time
}

UniformDeviate32::UniformDeviate32(const std::uint32_t seed) { Reseed(seed); }

// Set initial "state" using "seed".  If seed==0, then
// initialization is done using the system clock.
// NOTES on hexidecimal numbers:  a constant beginning
// with "0x" indicates a hexidecimal number (base 16).
// The "digits" of a hexidecimal number are
// 0-9,a,b,c,d,e,f (or upppercase).  A "u" or "U" then
// indicates it is unsigned, and an "L" indicates it is
// a long integer.

void UniformDeviate32::Reseed(const std::uint32_t seed) {
  std::uint32_t s = seed;
  int j;
  if (seed == 0) {
    printLaph("WARNING: Setting MT seed using system time");
    s = time(NULL);
    s ^= s << 8;
    s ^= s << 8;
    s ^= s << 8;
  }
  state[0] = s & 0xffffffffUL; // 2^32
  for (j = 1; j < N; j++) {
    state[j] = (1812433253UL * (state[j - 1] ^ (state[j - 1] >> 30)) + j);
    state[j] &= 0xffffffffUL; // for >32 bit machines
  }
  left = 0;
}

// works for arbitrary zn group to the same precision as the old code
static inline std::complex<double>
computeZ( const int idx , const int zn )
{
  long double cn = 0 , sn = 0 ;
  if( (2*idx)%zn == 0 ) {
    return std::complex( (double)(1-2*(2*idx)/zn) , 0.0 ) ;
  }
  if( (4*idx)%zn == 0 ) {
    return std::complex( 0.0 , (double)(2-(4*idx)/zn) ) ;
  }
  // generic case cast up to long double and back down to recover better precision
  sincosl( 2*idx*(M_PIl/(long double)zn) , &sn , &cn ) ;
  return std::complex( (double)cn , (double)sn ) ;
}

LaphZnNoise::LaphZnNoise(const int zn, const std::uint32_t &seed) : rng(seed), values(zn) {
  // generic version
  for( int i = 0 ; i < zn ; i++ ) {
    values[i] = computeZ( i , zn ) ;
  }
  switch( zn ) {
  case 1  :
    printLaph("Warning: ZN group set to N=1 for debugging ONLY");
    genptr = &LaphZnNoise::z1_generate  ;
    break ;
  case 4  : genptr = &LaphZnNoise::z4_generate  ; break ;
  case 8  : genptr = &LaphZnNoise::z8_generate  ; break ;
  case 32 : genptr = &LaphZnNoise::z32_generate ; break ;
  default :
    errorLaph("Unsupported Zn group in LaphZnNoise: only Z4, Z8, Z32");
    break ;
  }
  // call Mersenne Twister 12 times to start things
  for (int i = 0; i < 12; i++) {
    current = rng.generate();
  }
  znGroup = zn;
  count = 1;
};

// Must generate noise in **exactly** the same way every time.
// Store results in array to ensure optimizations do not omit terms in
// sequence.  This routine generates the full source, which will be used
// to compute a quark sink.

Array<std::complex<double>> LaphZnNoise::generateLapHQuarkSourceForSink(const int Textent,
							 const int Nspin,
                                                         const int nEigs) {
  Array<std::complex<double>> laph_noise(Textent, Nspin, nEigs);
  for (int t = 0; t < Textent; t++)
    for (int s = 0; s < Nspin; s++)
      for (int v = 0; v < nEigs; v++)
        laph_noise(t, s, v) = generate(); // same on all compute nodes
  return laph_noise;
}

// Must generate noise in **exactly** the same way every time.
// Store results in array to ensure optimizations do not omit terms in
// sequence.  This routine generates just a source: loop can stop when
// temporal loops hits timeValue since only noise(timeValue,spin,v)
// will be needed.

Array<std::complex<double>> LaphZnNoise::generateLapHQuarkSource(const int timeValue,
						  const int Textent,
                                                  const int Nspin,
						  const int nEigs) {
  Array<std::complex<double>> laph_noise(timeValue + 1, Nspin, nEigs);
  for (int t = 0; t <= timeValue; t++)
    for (int s = 0; s < Nspin; s++)
      for (int v = 0; v < nEigs; v++)
        laph_noise(t, s, v) = generate(); // same on all compute nodes
  return laph_noise;
}

} // namespace LaphEnv
