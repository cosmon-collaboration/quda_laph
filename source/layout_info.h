#ifndef LAYOUT_INFO_H
#define LAYOUT_INFO_H

#include <vector>
#include <map>
#include "xml_handler.h"


namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Class "LayoutInfo" is a singleton which stores and makes      *
// *   available information about the lattice and how the fields    *
// *   are laid out the lattice on the hosts (cpus).  The init(...)  *
// *   member should be called near the beginning of program         *
// *   execution, but after MPI and Quda communications have been    *
// *   initialized. A 4 dimensional lattice is assumed (make the 4th *
// *   dimension have size 1 to use a 3 dimensional lattice). QDP    *
// *   lexicographic format is used. The lattice is dividing up on   *
// *   the MPI ranks in the same way on each rank.  The number of    *
// *   lattice sites on each MPI rank is the same.  The lattice size *
// *   in each direction must be divisible by the MPI partition      *
// *   size in each direction.                                       *
// *                                                                 *
// *   The input XML should have the following form:                 *
// *                                                                 *
// *       <LatticeLayoutInfo>                                       *
// *          <XYZTExtents>24 24 24 48</XYTZExtents>                 *
// *       </LatticeLayoutInfo>                                      *
// *                                                                 *
// *   This class can assist with getting and putting the data on    *
// *   a particular site into the lattice appropriately.             *
// *                                                                 *
// *   In order to accommodate red/black checkerboarding for an      *
// *   efficient solve of the Dirac equation, the number of sites    *
// *   on each MPI rank must be even. Lattice sites are ordered as   *
// *   even-odd checkerboard (even/odd sites have x+y+z+t even/odd), *
// *   with ordering (x,y,z,t) with t varying most slowly for each   *
// *   parity. "start_parity" is the parity of the local (0,0,0,0)   *
// *   site on the global lattice.                                   *
// *                                                                 *
// *                                                                 *
// *******************************************************************



class LayoutInfo
{

  LayoutInfo();   // All constructors private

  ~LayoutInfo();  // Private destructor

  static std::vector<int> latt_extents;

  static std::vector<int> num_partitions;
  
  static std::vector<int> rank_latt_extents;
  
  static std::vector<int> my_comm_coords;
  
  static int my_rank;
 
  static int latt_nsites;

  static int rank_latt_nsites;
  
  static int num_ranks;
  
  static std::map<int,int> comm_map;
  
  static int my_start_parity;


 public:

  LayoutInfo(const LayoutInfo &in) = delete;           // no copy constructor

  LayoutInfo(LayoutInfo &in) = delete;                 // no copy constructor

  LayoutInfo& operator=(const LayoutInfo &in) = delete; // not assignable

  LayoutInfo& operator=(LayoutInfo &in) = delete;       // not assignable
    
     
  static void init(const XMLHandler& xml_in, const std::vector<int> npartitions, bool echo=true);

  static const std::vector<int>& getLattExtents()
   {return latt_extents;}

  static const std::vector<int>& getCommNumPartitions()
   {return num_partitions;}

  static const std::vector<int>& getRankLattExtents()
   {return rank_latt_extents;}

  static const std::vector<int>& getMyCommCoords() 
   {return my_comm_coords;}

  static const int& getMyRank()
   {return my_rank;}

  static const int& getLatticeNumSites()
   {return latt_nsites;}

  static const int& getRankLatticeNumSites()
   {return rank_latt_nsites;}

  static const int& getNumRanks()
   {return num_ranks;}
  
  static const int& getMyStartParity()
   {return my_start_parity;}

  static const int Ndim;

  static std::vector<int> getCommCoordsFromLatticeCoords(const std::vector<int>& latt_coords);  

  static int getRankFromLatticeCoords(const std::vector<int>& latt_coords);
  
  static void getCommInfoFromLatticeCoords(const std::vector<int>& latt_coords, 
                                           int& rank, int& rank_site_linear_index);

  static int getRankFromLexicoLinearIndex(int lexico_lin_index);

  static int getRankFromCommCoords(const std::vector<int>& comm_coords);

    // The input lattice quantity in "src" is assumed to be in lexicographic site order,
    // then "dest" is output in even-odd site order.

  static void lexico_to_evenodd(char* dest, const char* src, int sitebytes);

    // The input lattice quantity in "src" is assumed to be in even-odd site order,
    // then "dest" is output in  lexicographic site order.

  static void evenodd_to_lexico(char* dest, const char* src, int sitebytes);


 private:

  static void set_up_comm_map();
  
  static int get_rank_from_comm_coords(const std::vector<int>& comm_coords);

  static void get_comm_ranksite_coords(const std::vector<int>& latt_coords,
               std::vector<int>& comm_coords, std::vector<int>& rank_site_coords);
  
  static int rank_encode(const std::vector<int>& comm_coords);

  static int get_linear_index(int start_parity, int x, int y, int z, int t, 
                              int Nx, int Ny, int Nz, int loc_nsites);

  static int linearSiteIndex_lexico(const std::vector<int>& loc_coords, 
                                    const std::vector<int>& loc_extents);

  static int linearSiteIndex_evenodd(const std::vector<int>& loc_coords, 
                                     const std::vector<int>& loc_extents,
                                     int loc_nsites, int start_parity);

};


// *************************************************************************
}
#endif
