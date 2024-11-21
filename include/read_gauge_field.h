#ifndef READ_GAUGE_H
#define READ_GAUGE_H

#include "gauge_configuration_info.h"
#include "latt_field.h"
#include <cstring>
#include <string>
#include <vector>

namespace LaphEnv {

// ********************************************************************************
// * *
// *    An object of this class reads a 4-dimensional lattice SU3 gauge *
// *    configuration from file and puts the field into "U". *
// * *
// ********************************************************************************

class GaugeConfigReader {

public:
  GaugeConfigReader() {}

  GaugeConfigReader(const GaugeConfigReader &in) = delete;

  GaugeConfigReader &operator=(const GaugeConfigReader &in) = delete;

  ~GaugeConfigReader() {}

  bool read(std::vector<LattField> &U, XMLHandler &gauge_xmlinfo,
            const GaugeConfigurationInfo &ginfo);
};

// ********************************************************************************
// * *
// *  An object of this class reads a 4-dimensional lattice SU3 gauge *
// *  configuration in CERN format from a file named "cfg_file" and puts *
// *  the field into "u". *
// * *
// *  Details about CERN format: *
// *     - little endian, ints are 4 bytes, doubles are 8 bytes *
// *     - 4 ints NT,NX,NY,NZ *
// *     - 1 double for average plaquette *
// *     - links as SU3 matrices (row major, double precision complex) *
// *         in the following order: *
// *            - the 8 links in directions +0,-0,...,+3,-3 at the first *
// *              odd point, the second odd point, and so on. *
// *            - for the negative directions, this only refers to the *
// *              starting sites, but does NOT refer to the direction of the *
// *              parallel transport, so NO Hermitian conjugate is implied *
// *            - direction 0 refers to time, 1 to x, 2 to y, and 3 to z *
// *            - The order of the point (x0,x1,x2,x3) with *
// *               Cartesian coordinates in the range *
// *                 0<=x0<N0,...,0<=x3<N3 is determined by the index *
// *                   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0, *
// *               where N0,N1,N2,N3 are the global lattice sizes *
// *  Note that N0,N1,N2,N3 must all be EVEN for format to work. In summary, *
// *  the site ordering is (z,y,x,t) with left (z) fastest. *
// * *
// *  QDP lexicographic ordering is assumed in "u".  The links in the
// x-direction *
// *  will be returned in u[0], the y-dir links in u[1], the z-dir links in
// u[2], *
// *  and the temporal links will be in u[3].  In each, the site ordering is *
// *  (x,y,z,t) with x fastest varying.  At each site is an SU3 matrix in row *
// *  major format. *
// * *
// ********************************************************************************

class GaugeCERNConfigReader {

public:
  GaugeCERNConfigReader() {}

  GaugeCERNConfigReader(const GaugeCERNConfigReader &in) = delete;

  GaugeCERNConfigReader &operator=(const GaugeCERNConfigReader &in) = delete;

  ~GaugeCERNConfigReader() {}

  bool read(std::vector<LattField> &u, const std::string &cfg_file);

private:
  void su3_copy(double *dest, const double *src, int su3dble);

  void lexico_to_evenodd(std::vector<LattField> &u);

#ifdef ARCH_SERIAL

  void cern_to_qdp_lexico(const double *cern, double *xl, double *yl,
                          double *zl, double *tl, int NX, int NY, int NZ,
                          int NT, int su3dble);
#else

  void local_cern_to_qdp_lexico(const double *cernlinks, double *xl, double *yl,
                                double *zl, double *tl, double *xsend,
                                double *ysend, double *zsend, double *tsend,
                                int locNX, int locNY, int locNZ, int locNT,
                                bool start_parity, int su3dble);

  void local_cern_to_qdp_lexico_edge(const double *xedge, const double *yedge,
                                     const double *zedge, const double *tedge,
                                     double *xl, double *yl, double *zl,
                                     double *tl, int locNX, int locNY,
                                     int locNZ, int locNT, bool start_parity,
                                     int su3dble);

  void get_file_view(const std::vector<int> &global_sizes,
                     const std::vector<int> &local_sizes,
                     const std::vector<int> &rank_coords,
                     std::vector<int> &displacements, std::vector<int> &lengths,
                     int &start_parity);

#endif
};

inline void GaugeCERNConfigReader::su3_copy(double *dest, const double *src,
                                            int su3dble) {
  std::memcpy(dest, src, su3dble * sizeof(double));
}

// *************************************************************************
} // namespace LaphEnv
#endif
