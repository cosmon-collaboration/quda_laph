#include "layout_info.h"
#include "laph_stdio.h"
#include "utils.h"
#include <algorithm>
#include <cstring>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

using namespace std;

namespace LaphEnv {

std::vector<int> LayoutInfo::latt_extents;

std::vector<int> LayoutInfo::num_partitions;

std::vector<int> LayoutInfo::rank_latt_extents;

std::vector<int> LayoutInfo::my_comm_coords;

int LayoutInfo::my_rank;

int LayoutInfo::latt_nsites;

int LayoutInfo::rank_latt_nsites;

const int LayoutInfo::Ndim = 4;

int LayoutInfo::num_ranks;

std::map<int, int> LayoutInfo::comm_map;

int LayoutInfo::my_start_parity;

// *  Expects the following XML:
// *         <XYZTExtents>24 24 24 96</XYZTExtents>
// *         <Precision>...</Precision>  double or single (default: double)
// *  This should be called AFTER QUDA has set up communications.

void LayoutInfo::init(const XMLHandler &xmlin,
                      const std::vector<int> npartitions, bool echo) {
  XMLHandler xml_in(xmlin);
  if (latt_extents.size() == Ndim) {
    throw(std::runtime_error(
        "LayoutInfo::init can only be run ONCE! Don't try it again!"));
  }
  xmlreadchild(xml_in, "XYZTExtents", latt_extents, "LayoutInfo");
  if (latt_extents.size() != Ndim) {
    throw(std::runtime_error("XYZTExtents in invalid"));
  }
  for (uint k = 0; k < Ndim; ++k) {
    if (latt_extents[k] < 1) {
      throw(std::runtime_error("XYZTExtents in invalid"));
    }
  }
  if (npartitions.size() == Ndim) {
    num_partitions = npartitions;
  } else if (npartitions.size() == 0) {
    num_partitions.resize(Ndim);
    for (int k = 0; k < Ndim; ++k)
      num_partitions[k] = 1;
  } else {
    throw(std::runtime_error("Invalid communication partitioning"));
  }
  for (uint k = 0; k < Ndim; ++k) {
    if (num_partitions[k] < 1) {
      throw(std::runtime_error("XYZTCommNumPartitions in invalid"));
    }
  }
  int checkranks = 1;
  for (uint k = 0; k < Ndim; ++k) {
    checkranks *= num_partitions[k];
  }
  rank_latt_extents.resize(Ndim);
  my_comm_coords.resize(Ndim);
#ifdef ARCH_PARALLEL
  my_rank = comm_rank();
  num_ranks = comm_size();
  for (int k = 0; k < Ndim; ++k)
    my_comm_coords[k] = commCoords(k);
#else
  my_rank = 0;
  num_ranks = 1;
  for (int k = 0; k < Ndim; ++k)
    my_comm_coords[k] = 0;
#endif
  if (num_ranks != checkranks) {
    throw(std::runtime_error(
        "Partitioning and number of MPI ranks inconsistent"));
  }
  latt_nsites = 1;
  rank_latt_nsites = 1;
  for (int k = 0; k < Ndim; ++k) {
    latt_nsites *= latt_extents[k];
    rank_latt_extents[k] = latt_extents[k] / num_partitions[k];
    if ((rank_latt_extents[k] * num_partitions[k]) != latt_extents[k]) {
      errorLaph("Lattice is not partitioned into identical sublattices");
    }
    rank_latt_nsites *= rank_latt_extents[k];
  }
  if ((rank_latt_nsites % 2) == 1) {
    errorLaph("Number of sites on each MPI rank must be even to accommodate "
              "checkerboarding");
  }
  set_up_comm_map();
  my_start_parity = 0;
  for (int k = 0; k < Ndim; ++k) {
    my_start_parity += my_comm_coords[k] * rank_latt_extents[k];
  }
  my_start_parity = my_start_parity % 2;
  if (echo) {
#ifdef ARCH_PARALLEL
    if (npartitions.size() == Ndim) {
      xml_in.seek_first_child();
      xml_in.put_sibling("XYZTCommNumPartitions", make_string(npartitions));
      xml_in.put_sibling("XYZTExtentsEachRank", make_string(rank_latt_extents));
    }
#endif
    printLaph(make_strf("Lattice Layout Information:\n%s\n", xml_in.output()));
  }
}

void LayoutInfo::get_comm_ranksite_coords(const std::vector<int> &latt_coords,
                                          std::vector<int> &comm_coords,
                                          std::vector<int> &rank_site_coords) {
  comm_coords.resize(Ndim);
  rank_site_coords.resize(Ndim);
  for (int dir = 0; dir < Ndim; ++dir) {
    comm_coords[dir] = latt_coords[dir] / rank_latt_extents[dir];
    rank_site_coords[dir] =
        latt_coords[dir] - comm_coords[dir] * rank_latt_extents[dir];
  }
}

std::vector<int> LayoutInfo::getCommCoordsFromLatticeCoords(
    const std::vector<int> &latt_coords) {
  vector<int> comm_coords, rank_site_coords;
  get_comm_ranksite_coords(latt_coords, comm_coords, rank_site_coords);
  return comm_coords;
}

int LayoutInfo::getRankFromLatticeCoords(const std::vector<int> &latt_coords) {
  vector<int> comm_coords = getCommCoordsFromLatticeCoords(latt_coords);
  return get_rank_from_comm_coords(comm_coords);
}

/*      // this routine is for lexicographic site ordering
void LayoutInfo::getCommInfoFromLatticeCoords(const std::vector<int>&
latt_coords, int& rank, int& rank_site_linear_index)
{
 vector<int> comm_coords, rank_site_coords;
 get_comm_ranksite_coords(latt_coords,comm_coords,rank_site_coords);
 rank=get_rank_from_comm_coords(comm_coords);
 rank_site_linear_index=rank_site_coords[0]+rank_latt_extents[0]*(rank_site_coords[1]
                       +rank_latt_extents[1]*(rank_site_coords[2]+rank_latt_extents[2]*rank_site_coords[3]));
}
*/

// given a lattice site in "latt_coords", this determines which MPI rank
// handles this site, and what the linear index of that site is on that
// rank, assuming even-odd checkerboard

void LayoutInfo::getCommInfoFromLatticeCoords(
    const std::vector<int> &latt_coords, int &rank,
    int &rank_site_linear_index) {
  vector<int> comm_coords, rank_site_coords;
  get_comm_ranksite_coords(latt_coords, comm_coords, rank_site_coords);
  rank = get_rank_from_comm_coords(comm_coords);
  int start_parity = 0;
  for (int k = 0; k < Ndim; ++k) {
    start_parity += comm_coords[k] * rank_latt_extents[k];
  }
  start_parity = start_parity % 2;
  rank_site_linear_index = get_linear_index(
      start_parity, rank_site_coords[0], rank_site_coords[1],
      rank_site_coords[2], rank_site_coords[3], rank_latt_extents[0],
      rank_latt_extents[1], rank_latt_extents[2], rank_latt_nsites);
}

// assumes only that number of sites on each MPI process is even
int LayoutInfo::get_linear_index(int start_parity, int x, int y, int z, int t,
                                 int Nx, int Ny, int Nz, int loc_nsites) {
  return ((x + Nx * (y + Ny * (z + Nz * t))) / 2) +
         (loc_nsites / 2) * ((start_parity + x + y + z + t) % 2);
}

int LayoutInfo::getRankFromLexicoLinearIndex(int lexico_lin_index) {
  int index = lexico_lin_index;
  vector<int> latt_coords(Ndim);
  for (int k = 0; k < Ndim; ++k) {
    latt_coords[k] = index % latt_extents[k];
    index -= latt_coords[k];
    index /= latt_extents[k];
  }
  return getRankFromLatticeCoords(latt_coords);
}

int LayoutInfo::getRankFromCommCoords(const std::vector<int> &comm_coords) {
  return get_rank_from_comm_coords(comm_coords);
}

//  The comm_map stores information about which rank stores which part of the
//  lattice. All mpi ranks have the information.  The comm_coords is encoded to
//  a linear index which is used as the key in the map.

void LayoutInfo::set_up_comm_map() {
  vector<int> comm_codes(num_ranks, 0);
  vector<int> gl_comm_codes(num_ranks);
  comm_codes[my_rank] = rank_encode(my_comm_coords);
#ifdef ARCH_SERIAL
  comm_map[comm_codes[my_rank]] = my_rank;
#elif defined(ARCH_PARALLEL)
  int status = MPI_Allreduce(comm_codes.data(), gl_comm_codes.data(), num_ranks,
                             MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (status != MPI_SUCCESS) {
    errorLaph("Problem occurred in setting up the comm_map in LayoutInfo");
  }
  for (int rank = 0; rank < num_ranks; ++rank) {
    comm_map[gl_comm_codes[rank]] = rank;
  }
#endif
  // do a sanity check
  vector<int> check;
  vector<int> comm_coords(Ndim);
  for (comm_coords[0] = 0; comm_coords[0] < num_partitions[0]; ++comm_coords[0])
    for (comm_coords[1] = 0; comm_coords[1] < num_partitions[1];
         ++comm_coords[1])
      for (comm_coords[2] = 0; comm_coords[2] < num_partitions[2];
           ++comm_coords[2])
        for (comm_coords[3] = 0; comm_coords[3] < num_partitions[3];
             ++comm_coords[3]) {
          check.push_back(comm_map[rank_encode(comm_coords)]);
        }
  std::sort(check.begin(), check.end());
  for (int rank = 0; rank < num_ranks; ++rank) {
    if (check[rank] != rank) {
      throw(std::runtime_error(
          "Sanity check did not pass in setting up comm_map in LayoutInfo"));
    }
  }
}

int LayoutInfo::get_rank_from_comm_coords(const std::vector<int> &comm_coords) {
  return comm_map[rank_encode(comm_coords)];
}

int LayoutInfo::rank_encode(const std::vector<int> &comm_coords) {
  return comm_coords[0] +
         num_partitions[0] *
             (comm_coords[1] +
              num_partitions[1] *
                  (comm_coords[2] + num_partitions[2] * comm_coords[3]));
}

// The input lattice quantity in "src" is assumed to be in lexicographic site
// order, then "dest" is output in even-odd site order.

void LayoutInfo::lexico_to_evenodd(char *dest, const char *src, int sitebytes) {
  int locNX = LayoutInfo::getRankLattExtents()[0];
  int locNY = LayoutInfo::getRankLattExtents()[1];
  int locNZ = LayoutInfo::getRankLattExtents()[2];
  int locNT = LayoutInfo::getRankLattExtents()[3];
  const char *lex = src;
  char *even = dest;
  char *odd = dest + ((locNX * locNY * locNZ * locNT) / 2) * sitebytes;
  for (int t = 0; t < locNT; ++t)
    for (int z = 0; z < locNZ; ++z)
      for (int y = 0; y < locNY; ++y) {
        bool oddparity = (((t + z + y + my_start_parity) % 2) == 1);
        for (int x = 0; x < locNX; ++x, oddparity = !oddparity) {
          if (oddparity) {
            std::memcpy(odd, lex, sitebytes);
            odd += sitebytes;
          } else {
            std::memcpy(even, lex, sitebytes);
            even += sitebytes;
          }
          lex += sitebytes;
        }
      }
}

// The input lattice quantity in "src" is assumed to be in even-odd site order,
// then "dest" is output in  lexicographic site order.

void LayoutInfo::evenodd_to_lexico(char *dest, const char *src, int sitebytes) {
  int locNX = LayoutInfo::getRankLattExtents()[0];
  int locNY = LayoutInfo::getRankLattExtents()[1];
  int locNZ = LayoutInfo::getRankLattExtents()[2];
  int locNT = LayoutInfo::getRankLattExtents()[3];
  char *lex = dest;
  const char *even = src;
  const char *odd = src + ((locNX * locNY * locNZ * locNT) / 2) * sitebytes;
  for (int t = 0; t < locNT; ++t)
    for (int z = 0; z < locNZ; ++z)
      for (int y = 0; y < locNY; ++y) {
        bool oddparity = (((t + z + y + my_start_parity) % 2) == 1);
        for (int x = 0; x < locNX; ++x, oddparity = !oddparity) {
          if (oddparity) {
            std::memcpy(lex, odd, sitebytes);
            odd += sitebytes;
          } else {
            std::memcpy(lex, even, sitebytes);
            even += sitebytes;
          }
          lex += sitebytes;
        }
      }
}

int LayoutInfo::linearSiteIndex_lexico(const vector<int> &loc_coords,
                                       const vector<int> &loc_extents) {
  return loc_coords[0] +
         loc_extents[0] * (loc_coords[1] +
                           loc_extents[1] * (loc_coords[2] +
                                             loc_extents[2] * loc_coords[3]));
}

// assumes only that number of sites on each MPI process is even
int LayoutInfo::linearSiteIndex_evenodd(const vector<int> &loc_coords,
                                        const vector<int> &loc_extents,
                                        int loc_nsites, int start_parity) {
  return (loc_coords[0] +
          loc_extents[0] *
              (loc_coords[1] +
               loc_extents[1] *
                   (loc_coords[2] + loc_extents[2] * loc_coords[3]))) /
             2 +
         loc_nsites / 2 *
             ((start_parity + loc_coords[0] + loc_coords[1] + loc_coords[2] +
               loc_coords[3]) %
              2);
}

// **********************************************************************
} // namespace LaphEnv
