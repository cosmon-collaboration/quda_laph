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

LaphZnNoise::LaphZnNoise(const int zn, const std::uint32_t &seed) : rng(seed), values(zn) {
  const double one = 1.0, zero = 0.0;
  if (zn == 4) {
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(zero, one);
    values[2] = std::complex<double>(-one, zero);
    values[3] = std::complex<double>(zero, -one);
    genptr = &LaphZnNoise::z4_generate;
  } else if (zn == 8) {
    const double sqrthalf = 0.70710678118654752440;
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(sqrthalf, sqrthalf);
    values[2] = std::complex<double>(zero, one);
    values[3] = std::complex<double>(-sqrthalf, sqrthalf);
    values[4] = std::complex<double>(-one, zero);
    values[5] = std::complex<double>(-sqrthalf, -sqrthalf);
    values[6] = std::complex<double>(zero, -one);
    values[7] = std::complex<double>(sqrthalf, -sqrthalf);
    genptr = &LaphZnNoise::z8_generate;
  } else if (zn == 32) {
    const double c0 = 0.98078528040323044912; // cos(Pi/16)
    const double c1 = 0.92387953251128675613; // cos(Pi/8)
    const double c2 = 0.83146961230254523707; // cos(3*Pi/16)
    const double c3 = 0.70710678118654752440; // sqrt(1/2);
    const double c4 = 0.55557023301960222475; // cos(5*Pi/16)
    const double c5 = 0.38268343236508977173; // cos(3*Pi/8)
    const double c6 = 0.19509032201612826785; // cos(7*Pi/16)
    values[0] = std::complex<double>(one, zero);
    values[1] = std::complex<double>(c0, c6);
    values[2] = std::complex<double>(c1, c5);
    values[3] = std::complex<double>(c2, c4);
    values[4] = std::complex<double>(c3, c3);
    values[5] = std::complex<double>(c4, c2);
    values[6] = std::complex<double>(c5, c1);
    values[7] = std::complex<double>(c6, c0);
    values[8] = std::complex<double>(zero, one);
    values[9] = std::complex<double>(-c6, c0);
    values[10] = std::complex<double>(-c5, c1);
    values[11] = std::complex<double>(-c4, c2);
    values[12] = std::complex<double>(-c3, c3);
    values[13] = std::complex<double>(-c2, c4);
    values[14] = std::complex<double>(-c1, c5);
    values[15] = std::complex<double>(-c0, c6);
    values[16] = std::complex<double>(-one, zero);
    values[17] = std::complex<double>(-c0, -c6);
    values[18] = std::complex<double>(-c1, -c5);
    values[19] = std::complex<double>(-c2, -c4);
    values[20] = std::complex<double>(-c3, -c3);
    values[21] = std::complex<double>(-c4, -c2);
    values[22] = std::complex<double>(-c5, -c1);
    values[23] = std::complex<double>(-c6, -c0);
    values[24] = std::complex<double>(zero, -one);
    values[25] = std::complex<double>(c6, -c0);
    values[26] = std::complex<double>(c5, -c1);
    values[27] = std::complex<double>(c4, -c2);
    values[28] = std::complex<double>(c3, -c3);
    values[29] = std::complex<double>(c2, -c4);
    values[30] = std::complex<double>(c1, -c5);
    values[31] = std::complex<double>(c0, -c6);
    genptr = &LaphZnNoise::z32_generate;
  } else if (zn == 1) {
    printLaph("Warning: ZN group set to N=1 for debugging ONLY");
    values[0] = std::complex<double>(one, zero);
    genptr = &LaphZnNoise::z1_generate;
  } else {
    errorLaph("Unsupported Zn group in LaphZnNoise: only Z4, Z8, Z32");
  }

  // call Mersenne Twister 12 times to start things
  for (int i = 0; i < 11; i++)
    rng.generate();
  current = rng.generate();
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
