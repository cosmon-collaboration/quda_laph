#ifndef LAPH_NOISE_H
#define LAPH_NOISE_H

#include "array.h"
#include <cstdint>

namespace LaphEnv {

//                                                               *
//    Random number generator:  the 32-bit Mersenne Twister      
//    Using objects of class "UniformDeviate32":                 
//                                                               
//    UniformDeviate32 Q;        => defines Q (time-based seed)  
//    UniformDeviate32 Q(seed);  => or use explicit seed         
//    Q.Reseed(seed);  => explicit reseeding                     
//    Q.Generate();    => returns random 32-bit unsigned int     
//                                                               

// The celebrated Mersenne Twister.
// Mersenne Twister random number generator.
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v0.8  24 March 2002  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.
// It was designed with consideration of the flaws in various other
// generators. The period, 2^19937-1, and the order of equidistribution,
// 623 dimensions, are far greater.  The generator is also fast; it
// avoids multiplication and division, and it benefits from caches and
// pipelines.

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

class UniformDeviate32 {

private:
  static const int N = 624;     // length of state vector
  static const int M = 397;     // period parameter
  static const std::uint32_t MATRIX_A = 0x9908b0dfUL; // constant vector a
  static const std::uint32_t UMASK = 0x80000000UL;    // most significant w-r bits
  static const std::uint32_t LMASK =
      0x7fffffffUL; // least significant r bits
                    // hexadecimal constants start with 0x, UL=unsigned long

  std::uint32_t state[N]; // internal state
  std::uint32_t *next;    // next value to get from state
  int left;        // number of values left before next_state needed

public:
  UniformDeviate32(); // Default constructor, seeded by time.
  UniformDeviate32(
      std::uint32_t seed);      // Constructor with explicit seed (0 <= seed < 2^32).
  ~UniformDeviate32() {} // Destructor.
  void Reseed(std::uint32_t seed); // Explicitly re-seeds the generator (0 <= seed <
                            // 2^32   or  seed < 0 to seed by time).

  std::uint32_t generate(); // generates random 32 bit unsigned int

private:
  std::uint32_t twist(const std::uint32_t &m, const std::uint32_t &s0, const std::uint32_t &s1) const {
    return m ^ (((s0 & UMASK) | (s1 & LMASK)) >> 1) ^ (-(s1 & 1UL) & MATRIX_A);
  }
};

  inline std::uint32_t UniformDeviate32::generate() {
  if (left == 0) {
    int i;
    std::uint32_t *p = state;
    left = N;
    next = state;
    for (i = N - M + 1; --i; p++)
      *p = twist(p[M], p[0], p[1]);
    for (i = M; --i; p++)
      *p = twist(p[M - N], p[0], p[1]);
    *p = twist(p[M - N], p[0], state[0]);
  }
  left--;
  std::uint32_t s1 = *next++;
  s1 ^= (s1 >> 11);
  s1 ^= (s1 << 7) & 0x9d2c5680UL;
  s1 ^= (s1 << 15) & 0xefc60000UL;
  s1 ^= (s1 >> 18);
  return s1;
}

//                                                    *
//    Generates the Z_N noise in the Laph subspace.   *
//    Only Z[4], Z[8], and Z[32] are supported.       *
//    Usage:                                          *
//                                                    *
//        uint32 seed = 1293962291; ->  0 .. 2^32-1   *
//        int zngroup = 4;          -> 4, 8, or 32    *
//        LaphZnNoise rho(zngroup,seed);              *
//        rho.generate();  -> returns complex Z[N]    *
//                                                    *
//    Using N=1 is allowed and produces noise that    *
//    is always unity (for debugging purposes only).  *
//                                                    *

class LaphZnNoise {
  int znGroup;
  UniformDeviate32 rng;
  std::vector<std::complex<double>> values;
  std::uint32_t current;
  int count;
  std::complex<double> (LaphZnNoise::*genptr)();

public:
  LaphZnNoise(int zn, const std::uint32_t &seed);

  std::complex<double> generate() { return (this->*genptr)(); }

  Array<std::complex<double>> generateLapHQuarkSourceForSink(const int Textent,
							     const int Nspin,
							     const int nEigs);

  Array<std::complex<double>> generateLapHQuarkSource(const int timeValue,
						      const int Textent,
						      const int Nspin,
						      const int nEigs);

  int getZnGroup() const { return znGroup; }

private:
  // prevent copying
  LaphZnNoise(const LaphZnNoise &);
  LaphZnNoise &operator=(const LaphZnNoise &);

  std::complex<double> z1_generate(); // for debugging purposes
  std::complex<double> z4_generate();
  std::complex<double> z8_generate();
  std::complex<double> z32_generate();
};

inline std::complex<double> LaphZnNoise::z1_generate() { return values[0]; }

inline std::complex<double> LaphZnNoise::z4_generate() {
  if (count == 16) {
    current = rng.generate();
    count = 1;
  } else {
    current >>= 2;
    count++;
  }
  return values[current & 0x3UL];
}

inline std::complex<double> LaphZnNoise::z8_generate() {
  if (count == 10) {
    current = rng.generate();
    count = 1;
  } else {
    current >>= 3;
    count++;
  }
  return values[current & 0x7UL];
}

inline std::complex<double> LaphZnNoise::z32_generate() {
  if (count == 6) {
    current = rng.generate();
    count = 1;
  } else {
    current >>= 5;
    count++;
  }
  return values[current & 0x1FUL];
}
} // namespace LaphEnv
#endif
