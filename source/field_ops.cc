#include "field_ops.h"
#include "laph_stdio.h"
#include "utils.h"
#include <cstring>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

#include "QudaLaphBlas.h"

using namespace std;

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

std::vector<complex<double>>
getTimeSlicedInnerProducts(const LattField &leftfield,
                           const LattField &rightfield) {
  if ((leftfield.getFieldSiteType() != FieldSiteType::ColorVector) ||
      (rightfield.getFieldSiteType() != FieldSiteType::ColorVector)) {
    errorLaph("getTimeSlicedInnerProducts only supports ColorVector fields");
  }
  if (leftfield.bytesPerWord() != rightfield.bytesPerWord()) {
    errorLaph(
        "getTimeSlicedInnerProducts requires both fields have same precision");
  }
  bool dp = (leftfield.bytesPerWord() == sizeof(std::complex<double>));
  int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  int loc_npsites = loc_nsites / 2;
  int start_parity = LayoutInfo::getMyStartParity();
  int nloctime = LayoutInfo::getRankLattExtents()[3];
  int tstride = LayoutInfo::getRankLattExtents()[0] *
                LayoutInfo::getRankLattExtents()[1] *
                LayoutInfo::getRankLattExtents()[2];
  int cbytes =
      (dp) ? sizeof(std::complex<double>) : sizeof(std::complex<float>);
  int bps = leftfield.bytesPerSite();

  vector<char> iprods1(nloctime * cbytes);
  vector<char> iprods2(nloctime * cbytes);
  const char *lp = reinterpret_cast<const char *>(leftfield.getDataConstPtr());
  const char *rp = reinterpret_cast<const char *>(rightfield.getDataConstPtr());
  char *ip1 = reinterpret_cast<char *>(iprods1.data());
  char *ip2 = reinterpret_cast<char *>(iprods2.data());
  for (int tloc = 0; tloc < nloctime; ++tloc, ip1 += cbytes, ip2 += cbytes) {
    int parshift = loc_npsites * ((start_parity + tloc) % 2);
    int start1 = ((tstride * tloc) / 2) + parshift;
    int stop1 = ((1 + tstride * (tloc + 1)) / 2) + parshift;
    int n1 = FieldNcolor * (stop1 - start1);
    parshift = loc_npsites * ((start_parity + 1 + tloc) % 2);
    int start2 = ((1 + tstride * tloc) / 2) + parshift;
    int stop2 = ((tstride * (tloc + 1)) / 2) + parshift;
    int n2 = FieldNcolor * (stop2 - start2);
    const char *x1 = lp + bps * start1;
    const char *x2 = lp + bps * start2;
    const char *y1 = rp + bps * start1;
    const char *y2 = rp + bps * start2;
    if (dp) {
      cblas_zdotc_sub(n1, x1, 1, y1, 1, ip1);
      cblas_zdotc_sub(n2, x2, 1, y2, 1, ip2);
    } else {
      cblas_cdotc_sub(n1, x1, 1, y1, 1, ip1);
      cblas_cdotc_sub(n2, x2, 1, y2, 1, ip2);
    }
  }

  int ntime = LayoutInfo::getLattExtents()[3];
  vector<complex<double>> iprods(ntime);
  for (int t = 0; t < ntime; ++t) {
    iprods[t] = complex<double>(0.0, 0.0);
  }
  int mytmin =
      LayoutInfo::getMyCommCoords()[3] * LayoutInfo::getRankLattExtents()[3];
  if (dp) {
    complex<double> *z1 = reinterpret_cast<complex<double> *>(iprods1.data());
    complex<double> *z2 = reinterpret_cast<complex<double> *>(iprods2.data());
    for (int tloc = 0; tloc < nloctime; ++tloc, ++z1, ++z2) {
      iprods[tloc + mytmin] += (*z1) + (*z2);
    }
  } else {
    float *f1 = reinterpret_cast<float *>(iprods1.data());
    float *f2 = reinterpret_cast<float *>(iprods2.data());
    for (int tloc = 0; tloc < nloctime; ++tloc, ++f1, ++f2) {
      double zr = (*f1) + (*f2);
      ++f1;
      ++f2;
      double zi = (*f1) + (*f2);
      iprods[tloc + mytmin] += complex<double>(zr, zi);
    }
  }

#ifdef ARCH_PARALLEL
  vector<complex<double>> results(ntime);
  int status = MPI_Allreduce(iprods.data(), results.data(), 2 * ntime,
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in getTimeSlicedInnerProducts");
  }
  return results;
#else
  return iprods;
#endif
}

void setConstantField(LattField &field, const std::complex<double> &zconst) {
  bool dp = (field.bytesPerWord() == sizeof(complex<double>));
  int n = field.elemsPerSite() * LayoutInfo::getRankLatticeNumSites();
  if (dp) {
    complex<double> *z =
        reinterpret_cast<complex<double> *>(field.getDataPtr());
    for (int k = 0; k < n; ++k, ++z) {
      *z = zconst;
    }
  } else {
    complex<float> zfconst(float(real(zconst)), float(imag(zconst)));
    complex<float> *z = reinterpret_cast<complex<float> *>(field.getDataPtr());
    for (int k = 0; k < n; ++k, ++z) {
      *z = zfconst;
    }
  }
}

void setUnitField(LattField &field) {
  setConstantField(field, complex<double>(1.0, 0.0));
}

void setZeroField(LattField &field) {
  setConstantField(field, complex<double>(0.0, 0.0));
}

void compare_latt_fields(const LattField &src1, const LattField &src2) {
  bool flag = true;
  int n1 = src1.elemsPerSite();
  int n2 = src2.elemsPerSite();
  if (n1 != n2) {
    flag = false;
  } else {
    int s1 = src1.get_cpu_prec_bytes();
    int s2 = src2.get_cpu_prec_bytes();
    if (s1 != s2) {
      flag = false;
    } else {
      if (s1 == sizeof(complex<double>)) {
        double eps = 1e-9;
        const complex<double> *z1 =
            reinterpret_cast<const complex<double> *>(src1.getDataConstPtr());
        const complex<double> *z2 =
            reinterpret_cast<const complex<double> *>(src2.getDataConstPtr());
        int ncomp = n1 * LayoutInfo::getRankLatticeNumSites();
        for (int k = 0; k < ncomp; ++k, ++z1, ++z2) {
          if (std::abs((*z1) - (*z2)) > eps) {
            flag = false;
            break;
          }
        }
      } else {
        float eps = 1e-5;
        const complex<float> *z1 =
            reinterpret_cast<const complex<float> *>(src1.getDataConstPtr());
        const complex<float> *z2 =
            reinterpret_cast<const complex<float> *>(src2.getDataConstPtr());
        int ncomp = n1 * LayoutInfo::getRankLatticeNumSites();
        for (int k = 0; k < ncomp; ++k, ++z1, ++z2) {
          if (std::abs((*z1) - (*z2)) > eps) {
            flag = false;
            break;
          }
        }
      }
    }
  }
  if (globalAnd(flag)) {
    printLaph("Fields AGREE");
  } else
    printLaph("Fields DISAGREE");
}

template <typename T>
void su3color_mult(complex<T> *prod, const complex<T> *cdata1,
                   const complex<T> *cdata2, int krange, int kstride,
                   int istride, int jstride) {
  for (int k = 0; k < krange; ++k)
    for (int i = 0; i < FieldNcolor; ++i) {
      complex<T> z(0.0, 0.0);
      for (int j = 0; j < FieldNcolor; ++j) {
        z += cdata1[FieldNcolor * i + j] * cdata2[jstride * j + kstride * k];
      }
      prod[istride * i + kstride * k] = z;
    }
}

template <typename T>
void su3color_adjmult(complex<T> *prod, const complex<T> *cdata1,
                      const complex<T> *cdata2, int krange, int kstride,
                      int istride, int jstride) {
  for (int k = 0; k < krange; ++k)
    for (int i = 0; i < FieldNcolor; ++i) {
      complex<T> z(0.0, 0.0);
      for (int j = 0; j < FieldNcolor; ++j) {
        z += std::conj(cdata1[FieldNcolor * j + i]) *
             cdata2[jstride * j + kstride * k];
      }
      prod[istride * i + kstride * k] = z;
    }
}

//  Lattice-site-wise color-matrix multiplies outfield = fieldL * fieldR.
//  fieldL must be of type color-matrix. Any spin indices go along untouched.

void su3color_multiplier(LattField &outfield, const LattField &fieldL,
                         const LattField &fieldR, char Lmat) {
  if (fieldL.getFieldSiteType() != FieldSiteType::ColorMatrix) {
    errorLaph("left field must be ColorMatrix in su3color_mult");
  }
  if (fieldR.getFieldSiteType() == FieldSiteType::Complex) {
    errorLaph(
        "right field in su3color_mult cannot be complex (non-color) field");
  }
  outfield.reset(fieldR.getFieldSiteType());
  if (fieldR.bytesPerWord() != fieldL.bytesPerWord()) {
    errorLaph(
        "su3color_mult requires same precision between left and right fields");
  }
  int oinc = outfield.elemsPerSite();
  int rinc = fieldR.elemsPerSite();
  int linc = fieldL.elemsPerSite();
  int nsites = LayoutInfo::getRankLatticeNumSites();
  int kextent = 0, kstride = 0, jstride = 0, istride = 0;
  if (fieldR.getFieldSiteType() == FieldSiteType::ColorMatrix) {
    kstride = 1;
    kextent = 3;
    istride = 3;
    jstride = 3;
  } else if (fieldR.getFieldSiteType() == FieldSiteType::ColorVector) {
    kstride = 1;
    kextent = 1;
    istride = 1;
    jstride = 1;
  } else if (fieldR.getFieldSiteType() == FieldSiteType::ColorSpinVector) {
    kstride = 3;
    kextent = 4;
    istride = 1;
    jstride = 1;
  }
  if (fieldR.bytesPerWord() == sizeof(complex<double>)) {
    complex<double> *op =
        reinterpret_cast<complex<double> *>(outfield.getDataPtr());
    const complex<double> *rp =
        reinterpret_cast<const complex<double> *>(fieldR.getDataConstPtr());
    const complex<double> *lp =
        reinterpret_cast<const complex<double> *>(fieldL.getDataConstPtr());
    void (*multfunc)(complex<double> *, const complex<double> *,
                     const complex<double> *, int, int, int, int) =
        (Lmat == 'm') ? &su3color_mult<double> : &su3color_adjmult<double>;
    for (int ind = 0; ind < nsites; ++ind, op += oinc, rp += rinc, lp += linc) {
      multfunc(op, lp, rp, kextent, kstride, istride, jstride);
    }
  } else {
    complex<float> *op =
        reinterpret_cast<complex<float> *>(outfield.getDataPtr());
    const complex<float> *rp =
        reinterpret_cast<const complex<float> *>(fieldR.getDataConstPtr());
    const complex<float> *lp =
        reinterpret_cast<const complex<float> *>(fieldL.getDataConstPtr());
    void (*multfunc)(complex<float> *, const complex<float> *,
                     const complex<float> *, int, int, int, int) =
        (Lmat == 'm') ? &su3color_mult<float> : &su3color_adjmult<float>;
    for (int ind = 0; ind < nsites; ++ind, op += oinc, rp += rinc, lp += linc) {
      multfunc(op, lp, rp, kextent, kstride, istride, jstride);
    }
  }
}

//  Lattice-site-wise color-matrix multiplies outfield = fieldL * fieldR.
//  fieldL must be of type color-matrix. Any spin indices go along untouched.

void su3color_mult(LattField &outfield, const LattField &fieldL,
                   const LattField &fieldR) {
  su3color_multiplier(outfield, fieldL, fieldR, 'm');
}

//  Lattice-site-wise color-matrix multiplies outfield = colorAdj(fieldL) *
//  fieldR. fieldL must be of type color-matrix.  Any spin indices go along
//  untouched.

void su3color_adjmult(LattField &outfield, const LattField &fieldL,
                      const LattField &fieldR) {
  su3color_multiplier(outfield, fieldL, fieldR, 'a');
}

// assigns y[d] from x[d]; if d==dir, y[d]=(x[d]+1) % N[d]

void index_increaser(int dir, int d, vector<int> &y, const vector<int> &x,
                     const vector<int> &N) {
  if (dir == d) {
    y[d] = (x[d] == (N[d] - 1)) ? 0 : x[d] + 1;
  } else {
    y[d] = x[d];
  }
}

// assigns y[d] from x[d]; if d==dir, y[d]=(x[d]-1+N[d]) % N[d]

void index_decreaser(int dir, int d, vector<int> &y, const vector<int> &x,
                     const vector<int> &N) {
  if (dir == d) {
    y[d] = (x[d] == 0) ? N[d] - 1 : x[d] - 1;
  } else {
    y[d] = x[d];
  }
}

void flipsign(char *sitedata, int bps, bool dp) {
  if (dp) {
    int nelem = bps / sizeof(complex<double>);
    complex<double> *op = reinterpret_cast<complex<double> *>(sitedata);
    for (int k = 0; k < nelem; ++k, ++op) {
      *op = -(*op);
    }
  } else {
    int nelem = bps / sizeof(complex<float>);
    complex<float> *op = reinterpret_cast<complex<float> *>(sitedata);
    for (int k = 0; k < nelem; ++k, ++op) {
      *op = -(*op);
    }
  }
}

#ifdef ARCH_SERIAL

//  This routine is slow: meant only for debugging, testing
//   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
//   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

// linear index for site (x,y,z,t) is
//      (x+Nx*(y+Ny*(z+Nz*t)))/2 + (nsites/2)*((x+y+z+t)%2);

void latt_shifter(LattField &outfield, const LattField &infield, int dir,
                  void (*indexfunc)(int, int, vector<int> &,
                                    const vector<int> &, const vector<int> &),
                  char fwd_or_bwd, bool apply_antiperiodic_fermbc = false) {
  if ((dir < 0) || (dir > 3)) {
    errorLaph("Invalid direction in lattice field shift");
  }
  outfield.reset(infield.getFieldSiteType());
  const vector<int> &N = LayoutInfo::getLattExtents();
  int nsites = LayoutInfo::getLatticeNumSites();
  int npsites = nsites / 2;
  int bps = infield.bytesPerSite();
  char *oute = outfield.getDataPtr();
  char *outo = oute + npsites * bps;
  const char *ine = infield.getDataConstPtr();
  const char *ino = ine + npsites * bps;
  int fliptime = -1;
  if ((infield.getFieldSiteType() == FieldSiteType::ColorSpinVector) &&
      (apply_antiperiodic_fermbc) && (dir == 3)) {
    fliptime = (fwd_or_bwd == 'B') ? 0 : N[3] - 1;
  }
  bool dp = (infield.bytesPerWord() == sizeof(complex<double>));
  vector<int> x(LayoutInfo::Ndim);
  vector<int> y(LayoutInfo::Ndim);
  for (x[3] = 0; x[3] < N[3]; ++x[3]) {
    indexfunc(dir, 3, y, x, N);
    bool signflip = (x[3] == fliptime);
    for (x[2] = 0; x[2] < N[2]; ++x[2]) {
      indexfunc(dir, 2, y, x, N);
      for (x[1] = 0; x[1] < N[1]; ++x[1]) {
        indexfunc(dir, 1, y, x, N);
        for (x[0] = 0; x[0] < N[0]; ++x[0]) {
          indexfunc(dir, 0, y, x, N);
          int shift =
              bps * ((y[0] + N[0] * (y[1] + N[1] * (y[2] + N[2] * y[3]))) / 2);
          if ((x[0] + x[1] + x[2] + x[3]) % 2) {
            const char *inp = ine + shift;
            std::memcpy(outo, inp, bps);
            if (signflip) {
              flipsign(outo, bps, dp);
            }
            outo += bps;
          } else {
            const char *inp = ino + shift;
            std::memcpy(oute, inp, bps);
            if (signflip) {
              flipsign(oute, bps, dp);
            }
            oute += bps;
          }
        }
      }
    }
  }
}

//  This routine is slow: meant only for debugging, testing
//   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
//   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

void lattice_shift(LattField &outfield, const LattField &infield, int dir,
                   char fwd_or_bwd, bool apply_antiperiodic_fermbc = false) {
  if (fwd_or_bwd == 'F') {
    latt_shifter(outfield, infield, dir, &index_increaser, fwd_or_bwd,
                 apply_antiperiodic_fermbc);
  } else if (fwd_or_bwd == 'B') {
    latt_shifter(outfield, infield, dir, &index_decreaser, fwd_or_bwd,
                 apply_antiperiodic_fermbc);
  } else {
    errorLaph("lattice shift needs F or B for forward/backward");
  }
}

#else // parallel version now

bool fwd_stay_local(int dir, const vector<int> &x, const vector<int> &N,
                    bool nocomm) {
  if (x[dir] < (N[dir] - 1))
    return true;
  return nocomm;
}

bool bwd_stay_local(int dir, const vector<int> &x, const vector<int> &N,
                    bool nocomm) {
  if (x[dir] > 0)
    return true;
  return nocomm;
}

void get_fwd_comm_to_from(int dir, int &send_to, int &recv_from) {
  vector<int> rank_coord(
      LayoutInfo::getMyCommCoords()); // rank coords of this node
                                      // get mpi rank where to send
  if (rank_coord[dir] > 0)
    rank_coord[dir]--;
  else
    rank_coord[dir] = LayoutInfo::getCommNumPartitions()[dir] - 1;
  send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
  // get mpi rank where to receive from
  rank_coord = LayoutInfo::getMyCommCoords();
  rank_coord[dir]++;
  if (rank_coord[dir] == LayoutInfo::getCommNumPartitions()[dir])
    rank_coord[dir] = 0;
  recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);
}

void get_bwd_comm_to_from(int dir, int &send_to, int &recv_from) {
  vector<int> rank_coord(
      LayoutInfo::getMyCommCoords()); // rank coords of this node
                                      // get mpi rank where to receive from
  if (rank_coord[dir] > 0)
    rank_coord[dir]--;
  else
    rank_coord[dir] = LayoutInfo::getCommNumPartitions()[dir] - 1;
  recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);
  // get mpi rank where to send
  rank_coord = LayoutInfo::getMyCommCoords();
  rank_coord[dir]++;
  if (rank_coord[dir] == LayoutInfo::getCommNumPartitions()[dir])
    rank_coord[dir] = 0;
  send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
}

//  This routine is slow: meant only for debugging, testing
//   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
//   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

// linear index for site (x,y,z,t) is
//      (x+Nx*(y+Ny*(z+Nz*t)))/2 + (nsites/2)*((x+y+z+t)%2);

void latt_shifter(LattField &outfield, const LattField &infield, int dir,
                  void (*indexfunc)(int, int, vector<int> &,
                                    const vector<int> &, const vector<int> &),
                  bool (*staylocal)(int, const vector<int> &,
                                    const vector<int> &, bool),
                  void (*neighbors)(int, int &, int &), char fwd_or_bwd,
                  bool apply_antiperiodic_fermbc = false) {
  int start_parity = LayoutInfo::getMyStartParity();
  if ((dir < 0) || (dir > 3)) {
    errorLaph("Invalid direction in lattice field shift");
  }
  outfield.reset(infield.getFieldSiteType());
  const vector<int> &N = LayoutInfo::getRankLattExtents();
  int nsites = LayoutInfo::getRankLatticeNumSites();
  int npsites = nsites / 2;
  int bps = infield.bytesPerSite();
  char *oute = outfield.getDataPtr();
  char *outo = oute + npsites * bps;
  const char *ine = infield.getDataConstPtr();
  const char *ino = ine + npsites * bps;
  vector<int> x(LayoutInfo::Ndim);
  vector<int> y(LayoutInfo::Ndim);
  bool needcomm = (LayoutInfo::getCommNumPartitions()[dir] > 1);
  bool nocm = !needcomm;
  vector<char> sendbuffer, recvbuffer;
  char *snde = 0;
  char *sndo = 0;
  int nbuf = 0;
  bool nospflip = ((N[dir] % 2) == 0);
  if (needcomm) {
    nbuf = nsites / N[dir];
    nbuf += nbuf % 2;
    sendbuffer.resize(nbuf * bps);
    recvbuffer.resize(nbuf * bps);
    snde = sendbuffer.data();
    sndo = snde + (nbuf / 2) * bps;
    if (!nospflip) {
      sndo = sendbuffer.data();
      snde = sndo + (nbuf / 2) * bps;
    }
  }
  int fliptime = -1;
  if ((infield.getFieldSiteType() == FieldSiteType::ColorSpinVector) &&
      (apply_antiperiodic_fermbc) && (dir == 3)) {
    if ((fwd_or_bwd == 'B') && (LayoutInfo::getMyCommCoords()[3] ==
                                (LayoutInfo::getCommNumPartitions()[3] - 1))) {
      fliptime = 0;
    } else if ((fwd_or_bwd == 'F') && (LayoutInfo::getMyCommCoords()[3] == 0)) {
      fliptime = N[3] - 1;
    }
  }
  bool dp = (infield.bytesPerWord() == sizeof(complex<double>));
  // make local changes and put data in send buffer
  for (x[3] = 0; x[3] < N[3]; ++x[3]) {
    indexfunc(dir, 3, y, x, N);
    bool signflip = (x[3] == fliptime);
    for (x[2] = 0; x[2] < N[2]; ++x[2]) {
      indexfunc(dir, 2, y, x, N);
      for (x[1] = 0; x[1] < N[1]; ++x[1]) {
        indexfunc(dir, 1, y, x, N);
        for (x[0] = 0; x[0] < N[0]; ++x[0]) {
          indexfunc(dir, 0, y, x, N);
          int site_parity = (start_parity + x[0] + x[1] + x[2] + x[3]) % 2;
          int shift =
              bps * ((y[0] + N[0] * (y[1] + N[1] * (y[2] + N[2] * y[3]))) / 2);
          if (site_parity) {
            if (staylocal(dir, x, N, nocm)) {
              const char *inp = ine + shift;
              std::memcpy(outo, inp, bps);
              if (signflip) {
                flipsign(outo, bps, dp);
              }
            } else {
              const char *inp = (nospflip ? ine : ino) + shift;
              std::memcpy(sndo, inp, bps);
              if (signflip) {
                flipsign(sndo, bps, dp);
              }
              sndo += bps;
            }
            outo += bps;
          } else {
            if (staylocal(dir, x, N, nocm)) {
              const char *inp = ino + shift;
              std::memcpy(oute, inp, bps);
              if (signflip) {
                flipsign(oute, bps, dp);
              }
            } else {
              const char *inp = (nospflip ? ino : ine) + shift;
              std::memcpy(snde, inp, bps);
              if (signflip) {
                flipsign(snde, bps, dp);
              }
              snde += bps;
            }
            oute += bps;
          }
        }
      }
    }
  }
  if (!needcomm) {
    return;
  }

  // do the communication (send and receive)
  int send_to = 0;
  int recv_from = 0;
  neighbors(dir, send_to, recv_from);
  MPI_Status mpistatus;
  int status =
      MPI_Sendrecv(sendbuffer.data(), nbuf * bps, MPI_BYTE, send_to,
                   LayoutInfo::getMyRank(), recvbuffer.data(), nbuf * bps,
                   MPI_BYTE, recv_from, recv_from, MPI_COMM_WORLD, &mpistatus);
  if (status != MPI_SUCCESS) {
    errorLaph("Error in communication while doing shift");
  }

  // put data received into local memory appropriately
  sendbuffer.clear();
  oute = outfield.getDataPtr();
  outo = oute + npsites * bps;
  const char *rcve = recvbuffer.data();
  const char *rcvo = rcve + (nbuf / 2) * bps;
  for (x[3] = 0; x[3] < N[3]; ++x[3]) {
    for (x[2] = 0; x[2] < N[2]; ++x[2]) {
      for (x[1] = 0; x[1] < N[1]; ++x[1]) {
        for (x[0] = 0; x[0] < N[0]; ++x[0]) {
          if ((start_parity + x[0] + x[1] + x[2] + x[3]) % 2) {
            if (!staylocal(dir, x, N, nocm)) {
              std::memcpy(outo, rcvo, bps);
              rcvo += bps;
            }
            outo += bps;
          } else {
            if (!staylocal(dir, x, N, nocm)) {
              std::memcpy(oute, rcve, bps);
              rcve += bps;
            }
            oute += bps;
          }
        }
      }
    }
  }
}

//  This routine is slow: meant only for debugging, testing
//   outfield(x) <=  infield(x+dir)  if fwd_or_bwd=='F'
//   outfield(x) <=  infield(x-dir)  if fwd_or_bwd=='B'

void lattice_shift(LattField &outfield, const LattField &infield, int dir,
                   char fwd_or_bwd, bool apply_antiperiodic_fermbc = false) {
  if (fwd_or_bwd == 'F') {
    latt_shifter(outfield, infield, dir, &index_increaser, &fwd_stay_local,
                 &get_fwd_comm_to_from, fwd_or_bwd, apply_antiperiodic_fermbc);
  } else if (fwd_or_bwd == 'B') {
    latt_shifter(outfield, infield, dir, &index_decreaser, &bwd_stay_local,
                 &get_bwd_comm_to_from, fwd_or_bwd, apply_antiperiodic_fermbc);
  } else {
    errorLaph("lattice shift needs F or B for forward/backward");
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif
//   outfield(x) <=  U_dir(x) infield(x+mu)            if fwd_or_bwd=='F'
//   outfield(x) <=  U^dag_dir(x-mu) infield(x-mu)     if fwd_or_bwd=='B'

//   for 'F', su3mult( U[dir], shift(infield, mu, 'F') )
//   for 'B', shift(  su3mult( adj(U[dir]), infield ), mu, 'B')

void lattice_cov_shift(LattField &outfield, const LattField &infield,
                       const vector<LattField> &gauge_field, int dir,
                       char fwd_or_bwd,
                       bool apply_antiperiodic_fermbc = false) {
  if (int(gauge_field.size()) != LayoutInfo::Ndim) {
    errorLaph("invalid gauge field in lattice_cov_shift");
  }
  for (uint dir = 0; dir < gauge_field.size(); ++dir) {
    if (gauge_field[dir].getFieldSiteType() != FieldSiteType::ColorMatrix) {
      errorLaph("invalid gauge field in lattice_cov_shift");
    }
  }
  if (infield.getFieldSiteType() == FieldSiteType::Complex) {
    errorLaph("cannot covariantly shift a complex (non-color) field");
  }
  LattField tmp;
  if (fwd_or_bwd == 'F') {
    lattice_shift(tmp, infield, dir, 'F', apply_antiperiodic_fermbc);
    su3color_mult(outfield, gauge_field[dir], tmp);
  } else if (fwd_or_bwd == 'B') {
    su3color_adjmult(tmp, gauge_field[dir], infield);
    lattice_shift(outfield, tmp, dir, 'B', apply_antiperiodic_fermbc);
  } else {
    errorLaph("lattice shift needs F or B for forward/backward");
  }
}

//  covariant shift by a path of directions: each direction is an integer
//  having values    1,2,3,4 (forward directions) and -1,-2,-3,-4 (backwards)
//  where 1 means x, 2 means y, 3 means z, and 4 means t

void lattice_cov_shift(LattField &outfield, const LattField &infield,
                       const vector<LattField> &gauge_field,
                       vector<int> &path) {
  LattField tmp;
  LattField *p1 = 0, *p2 = 0, *sw = 0;
  const LattField *pc;
  if (path.size() % 2) {
    p2 = &outfield;
    pc = &infield;
  } else {
    p2 = &tmp;
    pc = &infield;
  }
  for (int seg = int(path.size()) - 1; seg >= 0; --seg) {
    lattice_cov_shift(*p2, *pc, gauge_field, abs(path[seg]) - 1,
                      path[seg] > 0 ? 'F' : 'B');
    if (seg == int(path.size()) - 1) {
      if (path.size() % 2) {
        p1 = &tmp;
      } else {
        p1 = &outfield;
      }
    }
    sw = p1;
    p1 = p2;
    p2 = sw;
    pc = p1;
  }
}

template <typename T>
void su3color_adjcopy(complex<T> *out, const complex<T> *in) {
  for (int i = 0; i < FieldNcolor; ++i)
    for (int j = 0; j < FieldNcolor; ++j) {
      out[FieldNcolor * i + j] = std::conj(in[FieldNcolor * j + i]);
    }
}

void su3color_adjcopy(LattField &outfield, const LattField &infield) {
  if (infield.getFieldSiteType() != FieldSiteType::ColorMatrix) {
    errorLaph("field must be ColorMatrix in su3color_adjcopy");
  }
  outfield.reset(FieldSiteType::ColorMatrix);
  int inc = infield.elemsPerSite();
  int nsites = LayoutInfo::getRankLatticeNumSites();
  if (infield.bytesPerWord() == sizeof(complex<double>)) {
    complex<double> *op =
        reinterpret_cast<complex<double> *>(outfield.getDataPtr());
    const complex<double> *ip =
        reinterpret_cast<const complex<double> *>(infield.getDataConstPtr());
    for (int ind = 0; ind < nsites; ++ind, op += inc, ip += inc) {
      su3color_adjcopy(op, ip);
    }
  } else {
    complex<float> *op =
        reinterpret_cast<complex<float> *>(outfield.getDataPtr());
    const complex<float> *ip =
        reinterpret_cast<const complex<float> *>(infield.getDataConstPtr());
    for (int ind = 0; ind < nsites; ++ind, op += inc, ip += inc) {
      su3color_adjcopy(op, ip);
    }
  }
}

void lattice_link_path(LattField &outfield, const vector<LattField> &gaugefield,
                       vector<int> &dir_path) {
  uint nshifts = dir_path.size();
  if (nshifts == 0) {
    errorLaph("lattice_link_path requires a non-trivial path");
  }
  LattField temp;
  vector<int>::const_reverse_iterator seg = dir_path.rbegin();
  bool doprod = false;
  for (uint k = 0; k < dir_path.size(); ++k) {
    int dir = *seg;
    if (dir > 0) {
      if (doprod) {
        lattice_shift(temp, outfield, dir - 1, 'F');
        su3color_mult(outfield, gaugefield[dir - 1], temp);
      } else {
        outfield = gaugefield[dir - 1];
        doprod = true;
      }
    } else {
      if (doprod) {
        su3color_adjmult(temp, gaugefield[-dir - 1], outfield);
      } else {
        su3color_adjcopy(temp, gaugefield[-dir - 1]);
        doprod = true;
      }
      lattice_shift(outfield, temp, -dir - 1, 'B');
    }
    ++seg;
  }
}

void lattice_addto(LattField &outfield, const LattField &infield,
                   const complex<double> &zcoef = complex<double>(1.0, 0.0)) {
  if (outfield.getFieldSiteType() != infield.getFieldSiteType()) {
    errorLaph("fields must be of same type to add");
  }
  if (outfield.bytesPerWord() != infield.bytesPerWord()) {
    errorLaph("addition requires same precision of fields");
  }
  int nelem = LayoutInfo::getRankLatticeNumSites() * infield.elemsPerSite();
  if (infield.bytesPerWord() == sizeof(complex<double>)) {
    complex<double> *op =
        reinterpret_cast<complex<double> *>(outfield.getDataPtr());
    const complex<double> *ip =
        reinterpret_cast<const complex<double> *>(infield.getDataConstPtr());
    for (int k = 0; k < nelem; ++k, ++op, ++ip) {
      *op += zcoef * (*ip);
    }
  } else {
    complex<float> *op =
        reinterpret_cast<complex<float> *>(outfield.getDataPtr());
    const complex<float> *ip =
        reinterpret_cast<const complex<float> *>(infield.getDataConstPtr());
    complex<float> zf(float(zcoef.real()), float(zcoef.imag()));
    for (int k = 0; k < nelem; ++k, ++op, ++ip) {
      *op += zf * (*ip);
    }
  }
}

//  *op = Gamma(spin_matrix_index) * (*ip)  site-wise
//
//        spin_matrix_index     matrix  (DeGrand-Rossi basis)
//               1              gamma[1]
//               2              gamma[2]     sigma[mu,nu] = i/2 [gamma[mu],
//               gamma[nu]] 3              gamma[3] 4              gamma[4] 5
//               gamma[5] = gamma[4]*gamma[1]*gamma[2]*gamma[3] 6 sigma[1,2] 7
//               sigma[1,3] 8              sigma[1,4] 9              sigma[2,3]
//              10              sigma[2,4]
//              11              sigma[3,4]

template <typename T>
void spin_mult(complex<T> *op, const complex<T> *ip,
               const vector<pair<int, complex<T>>> &spinmat) {
  int spinstride = FieldNcolor;
  complex<T> *opp = op;
  const complex<T> *ipp = ip;
  for (int color = 0; color < FieldNcolor; ++color, ++opp, ++ipp) {
    for (int outspin = 0; outspin < FieldNspin; ++outspin) {
      *(opp + outspin * spinstride) =
          *(ipp + spinmat[outspin].first * spinstride) *
          spinmat[outspin].second;
    }
  }
}

template <typename T>
void assign_spin_matrix(int spin_matrix_index,
                        vector<pair<int, complex<T>>> &spin_mat) {
  spin_mat.resize(4);
  complex<T> I(0.0, 1.0);
  complex<T> one(1.0, 0.0);
  if (spin_matrix_index == 1) {
    spin_mat[0] = pair<int, complex<T>>(3, I);
    spin_mat[1] = pair<int, complex<T>>(2, I);
    spin_mat[2] = pair<int, complex<T>>(1, -I);
    spin_mat[3] = pair<int, complex<T>>(0, -I);
  } else if (spin_matrix_index == 2) {
    spin_mat[0] = pair<int, complex<T>>(3, -one);
    spin_mat[1] = pair<int, complex<T>>(2, one);
    spin_mat[2] = pair<int, complex<T>>(1, one);
    spin_mat[3] = pair<int, complex<T>>(0, -one);
  } else if (spin_matrix_index == 3) {
    spin_mat[0] = pair<int, complex<T>>(2, I);
    spin_mat[1] = pair<int, complex<T>>(3, -I);
    spin_mat[2] = pair<int, complex<T>>(0, -I);
    spin_mat[3] = pair<int, complex<T>>(1, I);
  } else if (spin_matrix_index == 4) {
    spin_mat[0] = pair<int, complex<T>>(2, one);
    spin_mat[1] = pair<int, complex<T>>(3, one);
    spin_mat[2] = pair<int, complex<T>>(0, one);
    spin_mat[3] = pair<int, complex<T>>(1, one);
  } else if (spin_matrix_index == 5) {
    spin_mat[0] = pair<int, complex<T>>(0, -one);
    spin_mat[1] = pair<int, complex<T>>(1, -one);
    spin_mat[2] = pair<int, complex<T>>(2, one);
    spin_mat[3] = pair<int, complex<T>>(3, one);
  } else if (spin_matrix_index == 6) {
    spin_mat[0] = pair<int, complex<T>>(0, one);
    spin_mat[1] = pair<int, complex<T>>(1, -one);
    spin_mat[2] = pair<int, complex<T>>(2, one);
    spin_mat[3] = pair<int, complex<T>>(3, -one);
  } else if (spin_matrix_index == 7) {
    spin_mat[0] = pair<int, complex<T>>(1, -I);
    spin_mat[1] = pair<int, complex<T>>(0, I);
    spin_mat[2] = pair<int, complex<T>>(3, -I);
    spin_mat[3] = pair<int, complex<T>>(2, I);
  } else if (spin_matrix_index == 8) {
    spin_mat[0] = pair<int, complex<T>>(1, -one);
    spin_mat[1] = pair<int, complex<T>>(0, -one);
    spin_mat[2] = pair<int, complex<T>>(3, one);
    spin_mat[3] = pair<int, complex<T>>(2, one);
  } else if (spin_matrix_index == 9) {
    spin_mat[0] = pair<int, complex<T>>(1, one);
    spin_mat[1] = pair<int, complex<T>>(0, one);
    spin_mat[2] = pair<int, complex<T>>(3, one);
    spin_mat[3] = pair<int, complex<T>>(2, one);
  } else if (spin_matrix_index == 10) {
    spin_mat[0] = pair<int, complex<T>>(1, -I);
    spin_mat[1] = pair<int, complex<T>>(0, I);
    spin_mat[2] = pair<int, complex<T>>(3, I);
    spin_mat[3] = pair<int, complex<T>>(2, -I);
  } else if (spin_matrix_index == 11) {
    spin_mat[0] = pair<int, complex<T>>(0, -one);
    spin_mat[1] = pair<int, complex<T>>(1, one);
    spin_mat[2] = pair<int, complex<T>>(2, one);
    spin_mat[3] = pair<int, complex<T>>(3, -one);
  } else {
    errorLaph("Unsupported Dirac gamma spin index");
  }
}

void lattice_spin_multiply(LattField &outfield, const LattField &infield,
                           int spin_matrix_index) {
  if (infield.getFieldSiteType() != FieldSiteType::ColorSpinVector) {
    errorLaph("spin_multiply can only be done on a ColorSpinVector");
  }
  outfield.reset(FieldSiteType::ColorSpinVector);
  int nsites = LayoutInfo::getRankLatticeNumSites();
  int nelem = infield.elemsPerSite();
  if (infield.bytesPerWord() == sizeof(complex<double>)) {
    vector<pair<int, complex<double>>> spin_mat(4);
    assign_spin_matrix<double>(spin_matrix_index, spin_mat);
    complex<double> *op =
        reinterpret_cast<complex<double> *>(outfield.getDataPtr());
    const complex<double> *ip =
        reinterpret_cast<const complex<double> *>(infield.getDataConstPtr());
    for (int k = 0; k < nsites; ++k, op += nelem, ip += nelem) {
      spin_mult<double>(op, ip, spin_mat);
    }
  } else {
    vector<pair<int, complex<float>>> spin_mat(4);
    assign_spin_matrix<float>(spin_matrix_index, spin_mat);
    complex<float> *op =
        reinterpret_cast<complex<float> *>(outfield.getDataPtr());
    const complex<float> *ip =
        reinterpret_cast<const complex<float> *>(infield.getDataConstPtr());
    for (int k = 0; k < nsites; ++k, op += nelem, ip += nelem) {
      spin_mult<float>(op, ip, spin_mat);
    }
  }
}

void calcCloverLeaves(LattField &cloverleaf,
                      const vector<LattField> &gaugefield, int dir1, int dir2) {
  LattField Utmp;
  vector<int> path(4);
  // upper right leaf
  path[0] = dir1;
  path[1] = dir2;
  path[2] = -dir1;
  path[3] = -dir2;
  lattice_link_path(cloverleaf, gaugefield, path);
  // upper left leaf
  path[0] = dir2;
  path[1] = -dir1;
  path[2] = -dir2;
  path[3] = dir1;
  lattice_link_path(Utmp, gaugefield, path);
  lattice_addto(cloverleaf, Utmp);
  // lower left leaf
  path[0] = -dir1;
  path[1] = -dir2;
  path[2] = dir1;
  path[3] = dir2;
  lattice_link_path(Utmp, gaugefield, path);
  lattice_addto(cloverleaf, Utmp);
  // lower right leaf
  path[0] = -dir2;
  path[1] = dir1;
  path[2] = dir2;
  path[3] = -dir1;
  lattice_link_path(Utmp, gaugefield, path);
  lattice_addto(cloverleaf, Utmp);
}

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
                      const vector<LattField> &gauge_field,
                      const GaugeConfigurationInfo &gaction,
                      const QuarkActionInfo &qaction) {
  if (int(gauge_field.size()) != LayoutInfo::Ndim) {
    errorLaph("invalid gauge field in applyCloverDirac");
  }
  for (uint dir = 0; dir < gauge_field.size(); ++dir) {
    if (gauge_field[dir].getFieldSiteType() != FieldSiteType::ColorMatrix) {
      errorLaph("invalid gauge field in applyCloverDirac");
    }
  }
  if (infield.getFieldSiteType() != FieldSiteType::ColorSpinVector) {
    errorLaph(
        "can apply clover Dirac operator only to a color-spin vector field");
  }
  bool tbc1 = qaction.isFermionTimeBCAntiPeriodic();
  bool tbc2 = gaction.isFermionTimeBCAntiPeriodic();
  if (tbc1 != tbc2) {
    errorLaph("Inconsistent fermion time boundary conditions in "
              "QuarkActionInfo and GaugeConfigurationInfo",
              true);
  }
  outfield.reset(FieldSiteType::ColorSpinVector);

  if (qaction.getName() != "WILSON_CLOVER") {
    errorLaph("Can only applyCloverDirac if action name is WILSON_CLOVER");
  }
  bool timebc_antiperiodic = qaction.isFermionTimeBCAntiPeriodic();
  double kappa = qaction.getRValues()[2];
  double csw = qaction.getRValues()[4];
  double anisotropy = qaction.getRValues()[3];
  if (std::abs(anisotropy - 1.0) > 1e-12) {
    errorLaph("Current applyCloverDirac only applies for isotropic actions");
  }

  //  Now for the clover Dirac operator acting on a color-spin field
  LattField cloverleaf(FieldSiteType::ColorMatrix);
  LattField cloverleaf2(FieldSiteType::ColorMatrix);
  LattField phi(FieldSiteType::ColorSpinVector);
  LattField phi2(FieldSiteType::ColorSpinVector);
  LattField CFterm(FieldSiteType::ColorSpinVector);
  LattField Dterm(FieldSiteType::ColorSpinVector);
  setZeroField(CFterm);
  setZeroField(Dterm);
  int spin_index = 6;
  complex<double> cfclover(0.0, csw / 16.0);
  for (int dir1 = 1; dir1 <= LayoutInfo::Ndim; ++dir1)
    for (int dir2 = dir1 + 1; dir2 <= LayoutInfo::Ndim; ++dir2, ++spin_index) {
      calcCloverLeaves(cloverleaf, gauge_field, dir1, dir2);
      calcCloverLeaves(cloverleaf2, gauge_field, dir2, dir1);
      lattice_addto(cloverleaf, cloverleaf2, complex<double>(-1.0, 0.0));
      su3color_mult(phi, cloverleaf, infield);
      lattice_spin_multiply(phi2, phi, spin_index);
      lattice_addto(CFterm, phi2, cfclover);
    }

  //  Now the leading derivative and the Wilson term
  for (int dir = 1; dir <= 4; ++dir) {
    lattice_cov_shift(phi, infield, gauge_field, dir - 1, 'F',
                      timebc_antiperiodic);
    lattice_addto(Dterm, phi);
    lattice_spin_multiply(phi2, phi, dir);
    lattice_addto(Dterm, phi2, complex<double>(-1.0, 0.0));
    lattice_cov_shift(phi, infield, gauge_field, dir - 1, 'B',
                      timebc_antiperiodic);
    lattice_addto(Dterm, phi);
    lattice_spin_multiply(phi2, phi, dir);
    lattice_addto(Dterm, phi2);
  }

  setZeroField(outfield);
  lattice_addto(outfield, infield, complex<double>(1.0 / (2.0 * kappa), 0.0));
  lattice_addto(outfield, Dterm, complex<double>(-0.5, 0.0));

  // vector<vector<int>> sites;
  // sites.push_back(vector<int>{4,7,3,12});
  // sites.push_back(vector<int>{11,15,7,19});
  // printField(outfield,"Dirac term",sites);

  lattice_addto(outfield, CFterm);
  // printField(CFterm,"clover term",sites);
}

//  Applies the 3d spatial Laplacian with a smeared gauge field
//  onto "infield", returning result in "outfield".  These fields
//  must be color vectors

void applyMinusSpatialLaplacian(LattField &outfield, const LattField &infield,
                                const vector<LattField> &smeared_gauge_field) {
  if (int(smeared_gauge_field.size()) != LayoutInfo::Ndim) {
    errorLaph("invalid number of components in smeared gauge field in "
              "applySpatialLaplacian");
  }
  for (uint dir = 0; dir < smeared_gauge_field.size(); ++dir) {
    if (smeared_gauge_field[dir].getFieldSiteType() !=
        FieldSiteType::ColorMatrix) {
      errorLaph(make_strf("invalid field type in smeared gauge field[%d] in "
                          "applySpatialLaplacian",
                          dir));
    }
  }
  if (infield.getFieldSiteType() != FieldSiteType::ColorVector) {
    errorLaph(
        "can apply spatial laplacian operator only to a color vector field");
  }
  outfield.reset(FieldSiteType::ColorVector);

  setZeroField(outfield);
  lattice_addto(outfield, infield, complex<double>(6.0, 0.0));
  LattField phi(FieldSiteType::ColorVector);

  //  Apply the spatial covariant shifts
  for (int dir = 1; dir <= 3; ++dir) {
    lattice_cov_shift(phi, infield, smeared_gauge_field, dir - 1, 'F');
    lattice_addto(outfield, phi, complex<double>(-1.0, 0.0));
    lattice_cov_shift(phi, infield, smeared_gauge_field, dir - 1, 'B');
    lattice_addto(outfield, phi, complex<double>(-1.0, 0.0));
  }
}

// **************************************************************************
} // namespace LaphEnv
