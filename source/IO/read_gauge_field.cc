#include "QudaLaphIncludes.h"

#include "read_gauge_field.h"
#include "byte_handler.h"
#include "gauge_configuration_info.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include "util_quda.h"

using namespace std;
using namespace quda;

namespace LaphEnv {

// *************************************************************************

// The general reader: reads from file, assigns to "U',
// and extracts the gauge XML into "gauge_xmlinfo"

bool GaugeConfigReader::read(std::vector<LattField> &U,
                             XMLHandler &gauge_xmlinfo,
                             const GaugeConfigurationInfo &ginfo) {
  // if gauge config type is CERN
  if ((ginfo.file_format == "CERN") || (ginfo.file_format == "CLS")) {
    GaugeCERNConfigReader GCCR;
    GCCR.read(U, ginfo.getFileName());
  }

  return true;
}

// ************************************************************
// *                                                          *
// *                                                          *
// *                  CERN format read                        *
// *                                                          *
// *                                                          *
// ************************************************************

#ifdef ARCH_SERIAL

//  This routine loops through the links in CERN format (z(fastest),y,x,t) odd,
//  and rearranges them into QDP lexicographic order (x(fastest),y,z,t) all,
//  putting the links in the different directions into different arrays.  The
//  calling routine is responsible for allocating sufficient memory.  The CERN
//  links must be input in the array at location "cern".  The temporal links
//  will go into the array at "tl", and so on. "NX","NY","NZ","NT" are the sizes
//  of the lattice.  Each of the NX, NY, NZ, NT is assumed to be EVEN.

void GaugeCERNConfigReader::cern_to_qdp_lexico(
    const double *cern, double *xl, double *yl, double *zl, double *tl,
    const int NX, const int NY, const int NZ, const int NT, const int su3dble) {
  const int dx = 1 - NY * NZ;
  const int dy = NX - NZ;
  const int dz = NX * NY - 1;
  const int tstride = NX * NY * NZ;
  int t;
  // this could be collapsed but we would have to watch out for cindex increment
#pragma omp parallel for private(t)
  for (t = 0; t < NT; ++t) {
    int cindex =
        t * tstride; // calculate these so each loop independent for threading
    const double *cptr = cern + 4 * t * tstride * su3dble;
    const int td = (t > 0) ? (t - 1) : (NT - 1);
    const int tshift = (td - t) * NX * NY * NZ;
    for (int x = 0; x < NX; ++x) {
      const int xd = (x > 0) ? (x - 1) : (NX - 1);
      const int xshift = xd - x;
      for (int y = 0; y < NY; ++y) {
        const int yd = (y > 0) ? (y - 1) : (NY - 1);
        const int yshift = (yd - y) * NX;
        bool oddparity = (((t + x + y) % 2) == 1);
        for (int z = 0; z < NZ; ++z, ++cindex, oddparity = !oddparity) {
          if (oddparity) {
            const int zd = (z > 0) ? (z - 1) : (NZ - 1);
            const int zshift = (zd - z) * NX * NY;
            const int qindex = cindex + x * dx + y * dy + z * dz;
            const int qq = qindex * su3dble;
            su3_copy(tl + qq, cptr, su3dble);
            cptr += su3dble;
            su3_copy(tl + (qindex + tshift) * su3dble, cptr, su3dble);
            cptr += su3dble;
            su3_copy(xl + qq, cptr, su3dble);
            cptr += su3dble;
            su3_copy(xl + (qindex + xshift) * su3dble, cptr, su3dble);
            cptr += su3dble;
            su3_copy(yl + qq, cptr, su3dble);
            cptr += su3dble;
            su3_copy(yl + (qindex + yshift) * su3dble, cptr, su3dble);
            cptr += su3dble;
            su3_copy(zl + qq, cptr, su3dble);
            cptr += su3dble;
            su3_copy(zl + (qindex + zshift) * su3dble, cptr, su3dble);
            cptr += su3dble;
          }
        }
      }
    }
  }
}

bool GaugeCERNConfigReader::read(vector<LattField> &u,
                                 const std::string &in_cfg_file) {
  string cfg_file = tidyString(in_cfg_file);
  if (cfg_file.empty()) {
    errorLaph("Empty file name in GaugeCERNConfigReader::read");
  }
  if ((sizeof(int) != 4) || (sizeof(double) != 8)) {
    errorLaph("CERN files contain 4-byte ints, 8-byte doubles");
  }
  if (FieldNcolor != 3) {
    errorLaph("readCERN only supports Nc=3");
  }
  if (LayoutInfo::Ndim != 4) {
    errorLaph("readCERN only supported for 4 space-time dimensions");
  }

  const int NT = LayoutInfo::getLattExtents()[3];
  const int NZ = LayoutInfo::getLattExtents()[2]; // these should all be even
  const int NY = LayoutInfo::getLattExtents()[1];
  const int NX = LayoutInfo::getLattExtents()[0];
  if ((NT % 2) || (NX % 2) || (NY % 2) || (NZ % 2)) {
    errorLaph("Each lattice size in each direction must be even");
  }
  printLaph(
      make_strf("\nBeginning read of CERN gauge field (%d x %d x %d) x %d\n",
                NX, NY, NZ, NT));
  StopWatch rtimer;
  rtimer.start();
  StopWatch iotimer;

  //  CERN field always stored as little endian
  ByteHandler BH;
  const bool endian_convert = (BH.big_endian()) ? true : false;
  // open file for reading
  iotimer.start();
  ifstream fin(cfg_file.c_str(), std::ios::binary | std::ios::in);
  iotimer.stop();
  if (!fin) {
    errorLaph(
        make_strf("Error open CERN gauge configuration file %s\n", cfg_file));
  }
  //  read the lattice sizes
  int Ncern[4];
  iotimer.start();
  fin.read(reinterpret_cast<char *>(&Ncern[0]),
           4 * sizeof(int)); // assumes little endian, no checksums
  iotimer.stop();
  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(Ncern), sizeof(int), 4);
  }
  if ((Ncern[0] != NT) || (Ncern[1] != NX) || (Ncern[2] != NY) ||
      (Ncern[3] != NZ)) {
    errorLaph(make_strf("readCERN: Lattice size mismatch\nread %d %d %d "
                        "%d\nQudaLapH wants %d %d %d %d\n",
                        Ncern[1], Ncern[2], Ncern[3], Ncern[0], NX, NY, NZ,
                        NT));
  }

  // read the average plaquette
  double avgplaq = 0.0;
  iotimer.start();
  fin.read(reinterpret_cast<char *>(&avgplaq), sizeof(double));
  iotimer.stop();
  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(&avgplaq), sizeof(double), 1);
  }
  avgplaq /= 3.0;
  // now read the links
  const size_t linkdbles = 2 * FieldNcolor * FieldNcolor;
  const size_t ndir = LayoutInfo::Ndim;
  const size_t nelemsite = ndir * linkdbles;
  const size_t nelem = NX * NY * NZ * NT * nelemsite;
  vector<double> buffer(nelem);
  iotimer.start();
  fin.read(reinterpret_cast<char *>(buffer.data()), sizeof(double) * nelem);
  iotimer.stop();
  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(buffer.data()), sizeof(double),
                 nelem);
  }
  // now convert from CERN to QDP LEXICO format
  u.resize(LayoutInfo::Ndim);
  for (int dir = 0; dir < LayoutInfo::Ndim; ++dir) {
    u[dir].reset_by_precision(FieldSiteType::ColorMatrix, 'D');
  }
  double *xlinks = reinterpret_cast<double *>(u[0].getDataPtr());
  double *ylinks = reinterpret_cast<double *>(u[1].getDataPtr());
  double *zlinks = reinterpret_cast<double *>(u[2].getDataPtr());
  double *tlinks = reinterpret_cast<double *>(u[3].getDataPtr());

  cern_to_qdp_lexico(buffer.data(), xlinks, ylinks, zlinks, tlinks, NX, NY, NZ,
                     NT, linkdbles);

  // lastly, convert from lexico to even-odd, then to requested precision
  lexico_to_evenodd(u);

  rtimer.stop();
  printLaph(make_strf("      readCERN: plaq read: %16.12f\n", avgplaq));
  printLaph(make_strf("Read of CERN gauge field done in %g seconds\n",
                      rtimer.getTimeInSeconds()));
  printLaph(make_strf("Time of IO operations = %g seconds\n\n",
                      iotimer.getTimeInSeconds()));
  return true;
}

#else

//  This routine is called right after MPI-IO has been used to read from file
//  the part of the lattice assigned to this MPI process.  This routine loops
//  through the links of this sublattice in CERN format (z(fastest),y,x,t) odd,
//  and rearranges them into QDP lexicographic order (x(fastest),y,z,t), putting
//  the links in the different directions into different arrays.  The calling
//  routine is responsible for allocating sufficient memory.  The CERN links
//  must be input in the array at location "cernlinks".  The temporal links will
//  go into the array at "tl", and so on. "locNX","locNY","locNZ","locNT" are
//  the sizes of the sublattice on this MPI rank.  Edge case links which belong
//  to other MPI processes and must be sent to another MPI process are copied in
//  the send arrays pointed to by "xsend","ysend", "zsend","tsend".  If
//  "locNX==globalNX", then "xsend"=0 should be input; similar for other
//  directions. "start_parity" is the parity of the local (0,0,0,0) site on the
//  global lattice. "su3dble" is the number of doubles in an SU3 matrix.

void GaugeCERNConfigReader::local_cern_to_qdp_lexico(
    const double *cernlinks, double *xl, double *yl, double *zl, double *tl,
    double *xsend, double *ysend, double *zsend, double *tsend, const int locNX,
    const int locNY, const int locNZ, const int locNT, const bool start_parity,
    const int su3dble) {
  const int dx = 1 - locNY * locNZ;
  const int dy = locNX - locNZ;
  const int dz = locNX * locNY - 1;
  const int tstride = locNX * locNY * locNZ;
  double *ts = tsend;
  int t;
#pragma omp parallel for private(t)
  for (t = 0; t < locNT; ++t) {

    // for threading, each t must be independent
    int cindex =
        t * tstride; // calculate these so each loop independent for threading
    const double *cptr =
        cernlinks + 8 * ((t * tstride + start_parity) / 2) * su3dble;
    const int td = ((t > 0) || (tsend != 0))
                       ? (t - 1)
                       : (locNT - 1); // tsend==0 means locNT==globalNT
    const int tshift = (td - t) * locNX * locNY * locNZ;
    double *xs = xsend + t * ((locNY * locNZ + 1) / 2) * su3dble;
    double *ys = ysend + t * ((locNX * locNZ + 1) / 2) * su3dble;
    double *zs = zsend + t * ((locNY * locNX + 1) / 2) * su3dble;
    double *tdn = 0;
    double *zdn = 0;
    double *ydn = 0;
    double *xdn = 0;
    for (int x = 0; x < locNX; ++x) {
      const int xd = ((x > 0) || (xsend != 0))
                         ? (x - 1)
                         : (locNX - 1); // xsend==0 means locNX==globalNX
      const int xshift = xd - x;
      for (int y = 0; y < locNY; ++y) {
        const int yd = ((y > 0) || (ysend != 0))
                           ? (y - 1)
                           : (locNY - 1); // ysend==0 means locNY==globalNY
        const int yshift = (yd - y) * locNX;
        bool oddparity = (((t + x + y + start_parity) % 2) == 1);
        for (int z = 0; z < locNZ; ++z, ++cindex, oddparity = !oddparity) {
          if (oddparity) {
            const int zd = ((z > 0) || (zsend != 0))
                               ? (z - 1)
                               : (locNZ - 1); // zsend==0 means locNZ==globalNZ
            const int zshift = (zd - z) * locNX * locNY;
            const int qindex = cindex + x * dx + y * dy + z * dz;
            const int qq = qindex * su3dble;

            su3_copy(tl + qq, cptr, su3dble);
            cptr += su3dble;

            if (td >= 0) {
              tdn = tl + (qindex + tshift) * su3dble;
            } else {
              tdn = ts;
              ts += su3dble;
            }
            su3_copy(tdn, cptr, su3dble);
            cptr += su3dble;

            su3_copy(xl + qq, cptr, su3dble);
            cptr += su3dble;

            if (xd >= 0) {
              xdn = xl + (qindex + xshift) * su3dble;
            } else {
              xdn = xs;
              xs += su3dble;
            }
            su3_copy(xdn, cptr, su3dble);
            cptr += su3dble;

            su3_copy(yl + qq, cptr, su3dble);
            cptr += su3dble;

            if (yd >= 0) {
              ydn = yl + (qindex + yshift) * su3dble;
            } else {
              ydn = ys;
              ys += su3dble;
            }
            su3_copy(ydn, cptr, su3dble);
            cptr += su3dble;

            su3_copy(zl + qq, cptr, su3dble);
            cptr += su3dble;

            if (zd >= 0) {
              zdn = zl + (qindex + zshift) * su3dble;
            } else {
              zdn = zs;
              zs += su3dble;
            }
            su3_copy(zdn, cptr, su3dble);
            cptr += su3dble;
          }
        }
      }
    }
  }
}

//  This routine is the last step of reading the CERN config file and putting
//  the links into appropriate memory.  Before calling this routine, the
//  sublattices assigned to each MPI process have been read into memory on all
//  MPI processes, rearranged into QDP format, with edge links stored in
//  separate arrays.  These edge links have then been sent to their appropriate
//  MPI processes. This routine now loops through these received arrays and
//  completes the reordering of the links into QDP format.  The received links
//  should be available in the arrays pointed to by "xedge", "yedge", "zedge",
//  and "tedge", and the rearrangement is done into the arrays pointed to by
//  "xl", "yl", "zl", "tl".  "locNX","locNY", "locNZ","locNT" are the sizes of
//  the sublattice on this MPI rank.  "start_parity" is the parity of the local
//  (0,0,0,0) site on the global lattice. "su3dble" is the number of doubles in
//  an SU3 matrix.

void GaugeCERNConfigReader::local_cern_to_qdp_lexico_edge(
    const double *xedge, const double *yedge, const double *zedge,
    const double *tedge, double *xl, double *yl, double *zl, double *tl,
    const int locNX, const int locNY, const int locNZ, const int locNT,
    const bool start_parity, const int su3dble) {
  const int tpar = locNT % 2;
  const int xpar = locNX % 2;
  const int ypar = locNY % 2;
  const int zpar = locNZ % 2;
  const int dx = 1 - locNY * locNZ;
  const int dy = locNX - locNZ;
  const int dz = locNX * locNY - 1;
  const int tstride = locNX * locNY * locNZ;
  const double *tr = tedge;
  int t;
#pragma omp parallel for private(t)
  for (t = 0; t < locNT; ++t) {

    // for threading, each t must be independent
    int cindex =
        t * tstride; // calculate these so each loop independent for threading
    const int td = ((t == 0) && (tedge != 0)) ? (locNT - 1) : -1;
    const int tshift = (td - t) * locNX * locNY * locNZ;
    const double *xr = xedge + t * ((locNY * locNZ + 1) / 2) * su3dble;
    const double *yr = yedge + t * ((locNX * locNZ + 1) / 2) * su3dble;
    const double *zr = zedge + t * ((locNY * locNX + 1) / 2) * su3dble;
    double *tdn = 0;
    double *zdn = 0;
    double *ydn = 0;
    double *xdn = 0;
    for (int x = 0; x < locNX; ++x) {
      const int xd = ((x == 0) && (xedge != 0)) ? (locNX - 1) : -1;
      const int xshift = xd - x;
      for (int y = 0; y < locNY; ++y) {
        const int yd = ((y == 0) && (yedge != 0)) ? (locNY - 1) : -1;
        const int yshift = (yd - y) * locNX;
        bool oddparity = (((t + x + y + start_parity) % 2) == 1);
        for (int z = 0; z < locNZ; ++z, ++cindex, oddparity = !oddparity) {
          const int zd = ((z == 0) && (zedge != 0)) ? (locNZ - 1) : -1;
          const int zshift = (zd - z) * locNX * locNY;
          const int qindex = cindex + x * dx + y * dy + z * dz;

          if ((td >= 0) && ((oddparity + tpar) % 2)) {
            tdn = tl + (qindex + tshift) * su3dble;
            su3_copy(tdn, tr, su3dble);
            tr += su3dble;
          }

          if ((xd >= 0) && ((oddparity + xpar) % 2)) {
            xdn = xl + (qindex + xshift) * su3dble;
            su3_copy(xdn, xr, su3dble);
            xr += su3dble;
          }

          if ((yd >= 0) && ((oddparity + ypar) % 2)) {
            ydn = yl + (qindex + yshift) * su3dble;
            su3_copy(ydn, yr, su3dble);
            yr += su3dble;
          }

          if ((zd >= 0) && ((oddparity + zpar) % 2)) {
            zdn = zl + (qindex + zshift) * su3dble;
            su3_copy(zdn, zr, su3dble);
            zr += su3dble;
          }
        }
      }
    }
  }
}

//  This routine creates data need to initialize an MPI-IO file view of a CERN
//  config file for an MPI process having rank coordinates given in
//  "rank_coords". The global sizes of the entire lattice must be input in
//  "global_sizes", and the local sublattice sizes for this MPI process must be
//  input in "local_sizes".  This routine then returns "displacements",
//  "lengths", and "start_parity".  The file view is a sequences of blocks, the
//  blocks starting at serial indices specified in the vector "displacements",
//  and the blocks have lengths given in "lengths".  The values in
//  these two arrays are all in terms of CERN sites.  A CERN site is a
//  quantity having 8 su3 matrices in double precision complex.

void GaugeCERNConfigReader::get_file_view(const vector<int> &global_sizes,
                                          const vector<int> &local_sizes,
                                          const vector<int> &rank_coords,
                                          vector<int> &displacements,
                                          vector<int> &lengths,
                                          int &start_parity) {
  const int ndisp = local_sizes[1] * local_sizes[2] * local_sizes[3];
  displacements.resize(ndisp);
  lengths.resize(ndisp);
  vector<int> start_indices(4);
  for (int dir = 0; dir < 4; ++dir) {
    start_indices[dir] = rank_coords[dir] * local_sizes[dir];
  }
  start_parity = (start_indices[0] + start_indices[1] + start_indices[2] +
                  start_indices[3]) %
                 2;
  int count = 0;
  for (int g3 = start_indices[3]; g3 < (start_indices[3] + local_sizes[3]);
       ++g3)
    for (int g2 = start_indices[2]; g2 < (start_indices[2] + local_sizes[2]);
         ++g2)
      for (int g1 = start_indices[1]; g1 < (start_indices[1] + local_sizes[1]);
           ++g1, ++count) {
        int g0start = start_indices[0];
        int loc0start = 0;
        if (((g0start + g1 + g2 + g3) % 2) == 0) {
          ++g0start;
          loc0start = 1;
        }
        displacements[count] =
            (g0start +
             global_sizes[0] *
                 (g1 + global_sizes[1] * (g2 + global_sizes[2] * g3))) /
            2;
        lengths[count] = (local_sizes[0] + 1 - loc0start) / 2;
      }
}

//  This routine reads a 4-dimensional lattice SU3 gauge configuration in
//  CERN format from a file named "cfg_file" and puts the field into "u".
//
//  Details about CERN format:
//     - little endian, ints are 4 bytes, doubles are 8 bytes
//     - 4 ints NT,NX,NY,NZ
//     - 1 double for average plaquette
//     - links as SU3 matrices (row major, double precision complex)
//         in the following order:
//            - the 8 links in directions +0,-0,...,+3,-3 at the first
//              odd point, the second odd point, and so on.
//            - for the negative directions, this only refers to the
//              starting sites, but does NOT refer to the direction of the
//              parallel transport, so NO Hermitian conjugate is implied
//            - direction 0 refers to time, 1 to x, 2 to y, and 3 to z
//            - The order of the point (x0,x1,x2,x3) with
//               Cartesian coordinates in the range
//                 0<=x0<N0,...,0<=x3<N3 is determined by the index
//                   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
//               where N0,N1,N2,N3 are the global lattice sizes
//  Note that N0,N1,N2,N3 must all be EVEN for format to work. In summary,
//  the site ordering is (z,y,x,t) with left (z) fastest.
//
//  QDP lexicographic ordering is assumed in "u".  The links in the x-direction
//  will be returned in u[0], the y-dir links in u[1], the z-dir links in u[2],
//  and the temporal links will be in u[3].  In each, the site ordering is
//  (x,y,z,t) with x fastest varying.  At each site is an SU3 matrix in row
//  major format.

bool GaugeCERNConfigReader::read(vector<LattField> &u,
                                 const std::string &in_cfg_file) {
  string cfg_file = tidyString(in_cfg_file);
  if (cfg_file.empty()) {
    errorLaph("Empty file name in GaugeCERNConfigReader::read");
  }
  if ((sizeof(int) != 4) || (sizeof(double) != 8)) {
    errorLaph("CERN files contain 4-byte ints, 8-byte doubles");
  }
  if (FieldNcolor != 3) {
    errorLaph("readCERN only supports Nc=3");
  }
  if (LayoutInfo::Ndim != 4) {
    errorLaph("readCERN only supported for 4 space-time dimensions");
  }

  const int NT = LayoutInfo::getLattExtents()[3];
  const int NZ = LayoutInfo::getLattExtents()[2]; // these should all be even
  const int NY = LayoutInfo::getLattExtents()[1];
  const int NX = LayoutInfo::getLattExtents()[0];
  if ((NT % 2) || (NX % 2) || (NY % 2) || (NZ % 2)) {
    errorLaph("Each lattice size in each direction must be even");
  }
  if (LayoutInfo::getRankLatticeNumSites() % 2) {
    errorLaph(
        "Number of sites on each MPI rank must be even for checkerboarding\n");
  }
  printLaph(
      make_strf("\nBeginning read of CERN gauge field (%d x %d x %d) x %d", NX,
                NY, NZ, NT));
  StopWatch rtimer;
  rtimer.start();
  StopWatch iotimer;

  //  CERN field always stored as little endian
  ByteHandler BH;
  bool endian_convert = (BH.big_endian()) ? true : false;

  // open file for reading
  iotimer.start();
  MPI_File fh; // the MPI-IO file handler
  int status = MPI_File_open(MPI_COMM_WORLD, (const char *)cfg_file.c_str(),
                             MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  iotimer.stop();
  if (status != MPI_SUCCESS) {
    errorLaph(
        make_strf("Error open CERN gauge configuration file %s\n", cfg_file));
  }

  //  read the lattice sizes
  int Ncern[4];
  iotimer.start();
  MPI_Status mpistatus;
  status = MPI_File_read(fh, &Ncern[0], 4, MPI_INT, &mpistatus);
  iotimer.stop();
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "Bad read of lattice size from CERN gauge configuration file %s\n",
        cfg_file));
  }
  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(Ncern), sizeof(int), 4);
  }
  if ((Ncern[0] != NT) || (Ncern[1] != NX) || (Ncern[2] != NY) ||
      (Ncern[3] != NZ)) {
    errorLaph(make_strf("readCERN: Lattice size mismatch\nread %d %d %d "
                        "%d\nQudaLapH wants %d %d %d %d\n",
                        Ncern[1], Ncern[2], Ncern[3], Ncern[0], NX, NY, NZ,
                        NT));
  }

  // read the average plaquette
  double avgplaq = 0.0;
  iotimer.start();
  status = MPI_File_read(fh, &avgplaq, 1, MPI_DOUBLE, &mpistatus);
  iotimer.stop();
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "Bad read of lattice size from CERN gauge configuration file %s\n",
        cfg_file));
  }

  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(&avgplaq), sizeof(double), 1);
  }
  avgplaq /= 3.0;

  // now read the links in a partition of the lattice into memory
  const int ndir = LayoutInfo::Ndim;
  const int su3dble = 2 * FieldNcolor * FieldNcolor;
  const int ndblesite = ndir * su3dble;
  const int ndble = LayoutInfo::getRankLatticeNumSites() * ndblesite;
  vector<double> buffer(ndble);
  const int nbytes_per_site =
      2 * ndblesite * sizeof(double); // only odd sites, but 8 directions
  const int locNT = LayoutInfo::getRankLattExtents()[3];
  const int locNZ = LayoutInfo::getRankLattExtents()[2];
  const int locNY = LayoutInfo::getRankLattExtents()[1];
  const int locNX = LayoutInfo::getRankLattExtents()[0];
  vector<int> global_sizes(LayoutInfo::Ndim);
  vector<int> local_sizes(LayoutInfo::Ndim);
  vector<int> rank_coords(LayoutInfo::Ndim);
  for (int k = 0; k < LayoutInfo::Ndim; ++k) {
    global_sizes[k] = LayoutInfo::getLattExtents()[(LayoutInfo::Ndim + 2 - k) %
                                                   LayoutInfo::Ndim];
    local_sizes[k] =
        LayoutInfo::getRankLattExtents()[(LayoutInfo::Ndim + 2 - k) %
                                         LayoutInfo::Ndim];
    rank_coords[k] = LayoutInfo::getMyCommCoords()[(LayoutInfo::Ndim + 2 - k) %
                                                   LayoutInfo::Ndim];
  }
  int start_parity;
  vector<int> viewblockdisps;
  vector<int> viewblocklengths;
  get_file_view(global_sizes, local_sizes, rank_coords, viewblockdisps,
                viewblocklengths, start_parity);

  MPI_Datatype etype;
  MPI_Type_contiguous(nbytes_per_site, MPI_BYTE, &etype);
  MPI_Type_commit(&etype);
  MPI_Datatype ftype;
  status = MPI_Type_indexed(viewblockdisps.size(), viewblocklengths.data(),
                            viewblockdisps.data(), etype, &ftype);
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "Could not create needed filetype for reading CERN config file %s\n",
        cfg_file));
  }
  MPI_Type_commit(&ftype);
  MPI_Offset currdisp;
  status = MPI_File_get_position(fh, &currdisp);
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "File location error while reading CERN config file %s\n", cfg_file));
  }
  status = MPI_File_set_view(fh, currdisp, etype, ftype, (char *)"native",
                             MPI_INFO_NULL);
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "Failure during MPI_File_set_view while reading CERN config file %s\n",
        cfg_file));
  }
  const int readcount = LayoutInfo::getRankLatticeNumSites() / 2;
  iotimer.start();
  status = MPI_File_read_all(fh, buffer.data(), readcount, etype, &mpistatus);
  iotimer.stop();
  if (status == MPI_SUCCESS) {
    int count = 0;
    status = MPI_Get_count(&mpistatus, etype, &count);
    if (count != readcount)
      status = MPI_ERR_COUNT;
  }
  if (status != MPI_SUCCESS) {
    errorLaph(
        make_strf("Failure while reading CERN config file %s\n", cfg_file));
  }
  MPI_Type_free(&ftype);
  MPI_Type_free(&etype);
  status = MPI_File_close(&fh);
  if (status != MPI_SUCCESS) {
    errorLaph(make_strf(
        "Failure in closing the CERN gauge configuration file %s\n", cfg_file));
  }
  if (endian_convert) {
    BH.byte_swap(reinterpret_cast<char *>(buffer.data()), sizeof(double),
                 ndble);
  }

  // now convert from CERN to QDP LEXICO format
  u.resize(LayoutInfo::Ndim);
  for (int dir = 0; dir < LayoutInfo::Ndim; ++dir) {
    u[dir].reset_by_precision(FieldSiteType::ColorMatrix, 'D');
  }
  double *xlinks = reinterpret_cast<double *>(u[0].getDataPtr());
  double *ylinks = reinterpret_cast<double *>(u[1].getDataPtr());
  double *zlinks = reinterpret_cast<double *>(u[2].getDataPtr());
  double *tlinks = reinterpret_cast<double *>(u[3].getDataPtr());
  const int tsdim = ((locNX * locNY * locNZ + 1) / 2) * su3dble;
  const int xsdim = locNT * ((locNY * locNZ + 1) / 2) * su3dble;
  const int ysdim = locNT * ((locNX * locNZ + 1) / 2) * su3dble;
  const int zsdim = locNT * ((locNY * locNX + 1) / 2) * su3dble;
  vector<double> xsend;
  vector<double> ysend;
  vector<double> zsend;
  vector<double> tsend;
  double *xs = 0;
  double *ys = 0;
  double *zs = 0;
  double *ts = 0;
  if (locNT < NT) {
    tsend.resize(tsdim);
    ts = tsend.data();
  }
  if (locNX < NX) {
    xsend.resize(xsdim);
    xs = xsend.data();
  }
  if (locNY < NY) {
    ysend.resize(ysdim);
    ys = ysend.data();
  }
  if (locNZ < NZ) {
    zsend.resize(zsdim);
    zs = zsend.data();
  }

  local_cern_to_qdp_lexico(buffer.data(), xlinks, ylinks, zlinks, tlinks, xs,
                           ys, zs, ts, locNX, locNY, locNZ, locNT, start_parity,
                           su3dble);

  // now send edge links to the appropriate MPI processes
  double *xr = 0;
  double *yr = 0;
  double *zr = 0;
  double *tr = 0;
  buffer.clear();
  vector<double> trecv, xrecv, yrecv, zrecv;
  // send/recv edge case time links
  if (ts != 0) {
    trecv.resize(tsdim);
    tr = trecv.data();
    vector<int> rank_coord(
        LayoutInfo::getMyCommCoords()); // rank coords of this node
    const int ntproc = NT / locNT;
    const int tproc = rank_coord[3];
    // get mpi rank where to send
    rank_coord[3] = (tproc > 0) ? tproc - 1 : ntproc - 1;
    const int send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
    // get mpi rank where to receive from
    rank_coord[3] = (tproc < (ntproc - 1)) ? tproc + 1 : 0;
    const int recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);

    status = MPI_Sendrecv(ts, tsdim, MPI_DOUBLE, send_to,
                          LayoutInfo::getMyRank(), tr, tsdim, MPI_DOUBLE,
                          recv_from, recv_from, MPI_COMM_WORLD, &mpistatus);
    if (status != MPI_SUCCESS) {
      errorLaph(
          make_strf("Error in communcation while reading CERN config file %s\n",
                    cfg_file));
    }
  }
  tsend.clear();

  // send/recv edge case x-links
  if (xs != 0) {
    xrecv.resize(xsdim);
    xr = xrecv.data();
    vector<int> rank_coord(
        LayoutInfo::getMyCommCoords()); // rank coords of this node
    const int nxproc = NX / locNX;
    const int xproc = rank_coord[0];
    // get mpi rank where to send
    rank_coord[0] = (xproc > 0) ? xproc - 1 : nxproc - 1;
    const int send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
    // get mpi rank where to receive from
    rank_coord[0] = (xproc < (nxproc - 1)) ? xproc + 1 : 0;
    const int recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);

    status = MPI_Sendrecv(xs, xsdim, MPI_DOUBLE, send_to,
                          LayoutInfo::getMyRank(), xr, xsdim, MPI_DOUBLE,
                          recv_from, recv_from, MPI_COMM_WORLD, &mpistatus);
    if (status != MPI_SUCCESS) {
      errorLaph(
          make_strf("Error in communcation while reading CERN config file %s\n",
                    cfg_file));
    }
  }
  xsend.clear();

  // send/recv edge case y-links
  if (ys != 0) {
    yrecv.resize(ysdim);
    yr = yrecv.data();
    vector<int> rank_coord(
        LayoutInfo::getMyCommCoords()); // rank coords of this node
    const int nyproc = NY / locNY;
    const int yproc = rank_coord[1];
    // get mpi rank where to send
    rank_coord[1] = (yproc > 0) ? yproc - 1 : nyproc - 1;
    const int send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
    // get mpi rank where to receive from
    rank_coord[1] = (yproc < (nyproc - 1)) ? yproc + 1 : 0;
    const int recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);

    status = MPI_Sendrecv(ys, ysdim, MPI_DOUBLE, send_to,
                          LayoutInfo::getMyRank(), yr, ysdim, MPI_DOUBLE,
                          recv_from, recv_from, MPI_COMM_WORLD, &mpistatus);
    if (status != MPI_SUCCESS) {
      errorLaph(
          make_strf("Error in communcation while reading CERN config file %s\n",
                    cfg_file));
    }
  }
  ysend.clear();

  // send/recv edge case z-links
  if (zs != 0) {
    zrecv.resize(zsdim);
    zr = zrecv.data();
    vector<int> rank_coord(
        LayoutInfo::getMyCommCoords()); // rank coords of this node
    const int nzproc = NZ / locNZ;
    const int zproc = rank_coord[2];
    // get mpi rank where to send
    rank_coord[2] = (zproc > 0) ? zproc - 1 : nzproc - 1;
    const int send_to = LayoutInfo::getRankFromCommCoords(rank_coord);
    // get mpi rank where to receive from
    rank_coord[2] = (zproc < (nzproc - 1)) ? zproc + 1 : 0;
    const int recv_from = LayoutInfo::getRankFromCommCoords(rank_coord);

    status = MPI_Sendrecv(zs, zsdim, MPI_DOUBLE, send_to,
                          LayoutInfo::getMyRank(), zr, zsdim, MPI_DOUBLE,
                          recv_from, recv_from, MPI_COMM_WORLD, &mpistatus);
    if (status != MPI_SUCCESS) {
      errorLaph(
          make_strf("Error in communcation while reading CERN config file %s\n",
                    cfg_file));
    }
  }
  zsend.clear();

  // copy the edge case links to their final locations
  local_cern_to_qdp_lexico_edge(xr, yr, zr, tr, xlinks, ylinks, zlinks, tlinks,
                                locNX, locNY, locNZ, locNT, start_parity,
                                su3dble);

  // lastly, convert from lexico to even-odd, then to requested precision
  lexico_to_evenodd(u);

  rtimer.stop();
  printLaph(make_strf("      readCERN: plaq read: %16.12f\n", avgplaq));
  printLaph(make_strf("Read of CERN gauge field done in %g seconds\n",
                      rtimer.getTimeInSeconds()));
  printLaph(make_strf("Time of IO operations = %g seconds\n\n",
                      iotimer.getTimeInSeconds()));
  return true;
}

#endif

// convert from lexicographic site order to even-odd checkerboard as needed by
// quda

void GaugeCERNConfigReader::lexico_to_evenodd(vector<LattField> &u) {
  const int su3bytes = u[0].bytesPerSite();
  for (int dir = 0; dir < LayoutInfo::Ndim; ++dir) {
    char *links = (char *)(u[dir].getDataPtr());
    LattField temp;
    temp.reset_by_precision(FieldSiteType::ColorMatrix, 'D');
    LayoutInfo::lexico_to_evenodd((char *)(temp.getDataPtr()), links, su3bytes);
    u[dir].getDataRef() = std::move(temp.getDataRef());
    u[dir].to_quda_precision();
  }
}
} // namespace LaphEnv
