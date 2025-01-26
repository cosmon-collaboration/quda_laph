#include "layout_info.h"
#include "utils.h"
#include "laph_stdio.h"
#include "verbosity_info.h"
#include "util_quda.h"
#include <algorithm>
#include <cstring>

#ifdef ARCH_PARALLEL
#include <mpi.h>
#endif

#define QUDA_RESTRICT

using namespace std;

// Lexicographic ordering of sites:   (all variables local to the MPI rank)
//
//   index = x+Nx*(y+Ny*(z+Nz*t))
//
// Even/odd checkerboard ordering of sites:
//
//   index = ((x+Nx*(y+Ny*(z+Nz*t)))/2) + (loc_nsites/2)*((start_parity+x+y+z+t)%2);
//
// "start_parity" is the parity of local site (0,0,0,0) on the global lattice. 
// "loc_nsites" should be equal to Nx*Ny*Nz*Nt. Even/odd parity sites have 
// global coordinates gx+gy+gz+gt even/odd. The above formula assumes only that 
// the number of sites on each MPI rank is even.  Note that the local extents
// Nx,Ny,Nz,Nt are the same on all MPI ranks.
//
// The above lexicographic ordering is straightforward to understand.
// The data can be viewed as an ordering of contiguous blocks.  Each
// block contains Nx sites, corresponding to x=0,1,2,...,Nx-1.  The
// blocks are ordered according to (y,z,t), and there are Ny*Nz*Nt blocks.
// The block index is y+Ny*(z+Nz*t), so the full index is
//     x + block_size * block_index = x + Nx * (y+Ny*(z+Nz*t)).
//
// The checkerboard formula is straightforward to understand if Nx is EVEN.
// For even lexicographic block size Nx, then Nx/2 sites in each lexico block
// will belong to one parity (x=0,2,4,...,Nx-2), and the other Nx/2 sites 
// in a lexico block will belong to the other parity (x=1,3,5,...Nx-1). The 
// sequential index for each x will be (x/2), remembering how integer division
// works. For x=0,2,4,6..., xseq = 0,1,2,3,..., and for x=1,3,5,... we also
// have xseq=0,1,2,...  Then each parity is a sequence of Ny*Nz*Nt checkerboard
// blocks, each of block size Nx/2.  For a given parity, the sequential index is 
//      xseq + block_size * block_index = (x/2) + (Nx/2) * (y+Ny*(z+Nz*t))
//           = (x+Nx*(y+Ny*(z+Nz*t)))/2
// For a given (x,y,z,t), the last step in determining its full index is to
// determine its parity. If even, the full index calculation is done.  If odd,
// then Nx*Ny*Nz*Nt/2 must be added.  Remember, Nx*Ny*Nz*Nt/2 sites are
// even parity, and Nx*Ny*Nz*Nt/2 sites are odd parity.  We assume that
// Nx*Ny*Nz*Nt is EVEN, which means at least one of the Nx, Ny, Nz, Nt
// must be even.
//
// The checkerboard formula is more difficult to understand if Nx is ODD.
// In this case, the checkerboard block size is not a single value. Now,
// the sequence within a lexico block for one parity is (x=0,2,4,...,Nx-1), 
// which contains (Nx/2)+1 sites.  We refer to these as "L blocks", and
// L_block_size = (Nx/2)+1.   Then for the other parity, the sequence 
// within a lexico block is (x=1,3,5,...,Nx-2), which contains (Nx/2) sites.
// We refer to these as "S blocks", and S_block_size = (Nx/2).  For either
// type of checkerboard block, the sequence index for each x is still 
// xseq = (x/2).  Note that L_block_size + S_block_size = Nx.  An example 
// helps to see this.  Let Nx = 7, then the L block is  0,2,4,6 
// (size (7/2)+1 = 4), and the S block is 1,3,5 (size (7/2) = 3 ).  
// For each parity, the data is now a sequence which alternates 
// L and S blocks.  To maintain the same parity, an L block for (y,z,t)
// must become an S block for the next block index.
// The sequential index for a given parity is now
//     xseq + L_block_size * N_L + S_block_size * N_S,
// where N_L is the number of L blocks, N_S is the number of S blocks,
// and N_L + N_S = block_index = y+Ny*(z+Nz*t).
// Case 1:
//   If block_index = y+Ny*(z+Nz*t) is EVEN, 
//   then N_L = N_S = block_index/2  and   
//   the index for a given parity is
//      xseq + L_block_size * N_L + S_block_size * N_S 
//       =  (x/2) + (L_block_size + S_block_size) * (block_index/2)
//       =  (x/2) + Nx *  (block_index/2) 
//       =  (x/2) + ( (Nx * block_index) /2 ) 
//       =  ( x + Nx * block_index ) / 2
//   since Nx*block_index is even, so the formula still works
//   for x odd or even.
// Case 2:
//   If block_index is ODD and N_L = N_S + 1,
//   then N_S = block_index/2  and x is in an S-block, so x is ODD.
//   The index for a given parity is
//      xseq + L_block_size * N_L + S_block_size * N_S 
//       =  (x/2) + L_block_size * (N_S + 1) + S_block_size * N_S 
//       =  (x/2) + (L_block_size + S_block_size) * N_S + L_block_size
//       =  (x/2) + Nx * (block_index/2) + (Nx/2) + 1
//   where x, Nx, block_index are all ODD.  Continuing,
//       =  (x/2) + Nx * ((block_index-1)/2) + ((Nx-1)/2) + (2/2)
//       =  (x/2) + (  Nx * (block_index-1) + (Nx-1) + 2  ) /2
//       =  (x/2) + (  Nx * block_index + 1  ) /2
//       =  (x +  Nx * block_index + 1  ) /2
//       =  ( x +  Nx * block_index  ) /2
//   since Nx*block_index and x are ODD, so x+Nx*block_index is even.
// Case 3:
//   If block_index is ODD and N_S = N_L + 1,
//   then N_L = block_index/2  and x is in an L-block, so x is EVEN
//   The index for a given parity is
//      xseq + L_block_size * N_L + S_block_size * N_S 
//       =  (x/2) + L_block_size * N_L + S_block_size * (N_L + 1) 
//       =  (x/2) + ( L_block_size + S_block_size ) * N_L  + S_block_size 
//       =  (x/2) + Nx * (block_index/2) + (Nx/2)
//   where Nx, block_index are ODD, but x is even.  Continuing,
//       =  (x/2) + Nx * ((block_index-1)/2) + ((Nx-1)/2)
//       =  (x/2) + (Nx * (block_index-1) + Nx-1 )/2
//       =  (x/2) + (Nx * block_index - 1 )/2
//       =  (x + Nx * block_index - 1 )/2
//       =  (x + Nx * block_index )/2
//   since x+Nx*block_index is odd. 
// Summary: the above formula works for Nx odd.  The only requirement is
// that the number of sites on each MPI rank must be even.


namespace LaphEnv {


std::vector<int> LayoutInfo::latt_extents;

std::vector<int> LayoutInfo::num_partitions;

std::vector<int> LayoutInfo::rank_latt_extents;
 
std::vector<int> LayoutInfo::my_comm_coords;

int LayoutInfo::my_rank;

int LayoutInfo::latt_nsites;

int LayoutInfo::rank_latt_nsites;

const int LayoutInfo::Ndim=4;

int LayoutInfo::num_ranks;

std::map<int,int> LayoutInfo::comm_map;

int LayoutInfo::my_start_parity;


// *  Expects the following XML:
// *         <XYZTExtents>24 24 24 96</XYZTExtents>
// *  This should be called AFTER QUDA has set up communications.

void LayoutInfo::init(const XMLHandler& xmlin, const std::vector<int> npartitions, bool echo)
{
 XMLHandler xml_in(xmlin);
 if (latt_extents.size()==Ndim){
    throw(std::runtime_error("LayoutInfo::init can only be run ONCE! Don't try it again!"));}
 xmlreadchild(xml_in,"XYZTExtents",latt_extents,"LayoutInfo");
 if (latt_extents.size()!=Ndim){
    throw(std::runtime_error("XYZTExtents in invalid"));}
 for (uint k=0;k<Ndim;++k){
    if (latt_extents[k]<1){
       throw(std::runtime_error("XYZTExtents in invalid"));}}
 if (npartitions.size()==Ndim){
    num_partitions=npartitions;}
 else if (npartitions.size()==0){
    num_partitions.resize(Ndim);
    for (int k=0;k<Ndim;++k) num_partitions[k]=1;}
 else{
    throw(std::runtime_error("Invalid communication partitioning"));}
 for (uint k=0;k<Ndim;++k){
    if (num_partitions[k]<1){
       throw(std::runtime_error("XYZTCommNumPartitions in invalid"));}}
 int checkranks=1;
 for (uint k=0;k<Ndim;++k){
    checkranks*=num_partitions[k];}
 rank_latt_extents.resize(Ndim);
 my_comm_coords.resize(Ndim);
#ifdef ARCH_PARALLEL
 my_rank=comm_rank();
 num_ranks=comm_size();
 for (int k=0;k<Ndim;++k) my_comm_coords[k]=commCoords(k);
#else
 my_rank=0;
 num_ranks=1;
 for (int k=0;k<Ndim;++k) my_comm_coords[k]=0;
#endif
 if (num_ranks!=checkranks){
    throw(std::runtime_error("Partitioning and number of MPI ranks inconsistent"));}
 latt_nsites=1;
 rank_latt_nsites=1;
 for (int k=0;k<Ndim;++k){
    latt_nsites*=latt_extents[k];
    rank_latt_extents[k]=latt_extents[k]/num_partitions[k];
    if ((rank_latt_extents[k]*num_partitions[k])!=latt_extents[k]){
       errorLaph("Lattice is not partitioned into identical sublattices");}
#ifdef QUDA_RESTRICT
    if (rank_latt_extents[k]<4){
       errorLaph("An extent in any direction on each MPI rank cannot be less than 4");}
    if (rank_latt_extents[k]%2){
       errorLaph("QUDA does not support an odd extent in any direction on each MPI rank");}
#endif
    rank_latt_nsites*=rank_latt_extents[k];}
 if ((rank_latt_nsites%2)==1){
    errorLaph("Number of sites on each MPI rank must be even to accommodate checkerboarding");}
 set_up_comm_map();
 my_start_parity=0;
 for (int k=0;k<Ndim;++k){
    my_start_parity+=my_comm_coords[k]*rank_latt_extents[k];}
 my_start_parity=my_start_parity%2;
 if (echo){
#ifdef ARCH_PARALLEL
    if (npartitions.size()==Ndim){
       xml_in.seek_first_child();
       xml_in.put_sibling("XYZTCommNumPartitions",make_string(npartitions));
       xml_in.put_sibling("XYZTExtentsEachRank",make_string(rank_latt_extents));}
#endif
    printLaph(make_strf("Lattice Layout Information:\n%s\n",xml_in.output()));}
}



void LayoutInfo::get_comm_ranksite_coords(const std::vector<int>& latt_coords,
                                          std::vector<int>& comm_coords, 
                                          std::vector<int>& rank_site_coords)
{
 comm_coords.resize(Ndim);
 rank_site_coords.resize(Ndim);
 for (int dir=0;dir<Ndim;++dir){
    comm_coords[dir]=latt_coords[dir]/rank_latt_extents[dir];
    rank_site_coords[dir]=latt_coords[dir]-comm_coords[dir]*rank_latt_extents[dir];}
}


std::vector<int> LayoutInfo::getCommCoordsFromLatticeCoords(const std::vector<int>& latt_coords)
{
 vector<int> comm_coords, rank_site_coords;
 get_comm_ranksite_coords(latt_coords,comm_coords,rank_site_coords);
 return comm_coords;
}


int LayoutInfo::getRankFromLatticeCoords(const std::vector<int>& latt_coords)
{
 vector<int> comm_coords=getCommCoordsFromLatticeCoords(latt_coords);
 return get_rank_from_comm_coords(comm_coords);
}
 
/*      // this routine is for lexicographic site ordering
void LayoutInfo::getCommInfoFromLatticeCoords(const std::vector<int>& latt_coords, 
                                              int& rank, int& rank_site_linear_index)
{
 vector<int> comm_coords, rank_site_coords;
 get_comm_ranksite_coords(latt_coords,comm_coords,rank_site_coords);
 rank=get_rank_from_comm_coords(comm_coords);
 rank_site_linear_index=rank_site_coords[0]+rank_latt_extents[0]*(rank_site_coords[1]
                       +rank_latt_extents[1]*(rank_site_coords[2]+rank_latt_extents[2]*rank_site_coords[3]));
}
*/

       // Given a lattice site in "latt_coords", this determines which MPI rank
       // handles this site, and what the linear index of that site is on that
       // rank, assuming even-odd checkerboard

void LayoutInfo::getCommInfoFromLatticeCoords(const std::vector<int>& latt_coords, 
                                              int& rank, int& rank_site_linear_index)
{
 vector<int> comm_coords, rank_site_coords;
 get_comm_ranksite_coords(latt_coords,comm_coords,rank_site_coords);
 rank=get_rank_from_comm_coords(comm_coords);
 int start_parity=0;
 for (int k=0;k<Ndim;++k){
    start_parity+=comm_coords[k]*rank_latt_extents[k];}
 start_parity=start_parity%2;
 rank_site_linear_index=get_linear_index(start_parity,
             rank_site_coords[0],rank_site_coords[1],rank_site_coords[2],rank_site_coords[3],
             rank_latt_extents[0],rank_latt_extents[1],rank_latt_extents[2],rank_latt_nsites);
}


     // Returns the linear index of the site (x,y,z,t) for this MPI rank,
     // where the coordinates are LOCAL to this MPI rank. Note that Nx,Ny,Nz,Nt
     // are the local extents for this MPI rank. start_parity is the parity of
     // local site (0,0,0,0) on the global lattice. loc_nsites should be equal to
     // Nx*Ny*Nz*Nt.  The formula used assumes only that the number of sites 
     // on each MPI process is even.

int LayoutInfo::get_linear_index(int start_parity, int x, int y, int z, int t, 
                                 int Nx, int Ny, int Nz, int loc_nsites)
{
 return ((x+Nx*(y+Ny*(z+Nz*t)))/2) + (loc_nsites/2)*((start_parity+x+y+z+t)%2);
}


int LayoutInfo::getRankFromLexicoLinearIndex(int lexico_lin_index)
{
 int index=lexico_lin_index;
 vector<int> latt_coords(Ndim);
 for (int k=0;k<Ndim;++k){
    latt_coords[k]=index % latt_extents[k];
    index-=latt_coords[k]; index/=latt_extents[k];}
 return getRankFromLatticeCoords(latt_coords);
}


int LayoutInfo::getRankFromCommCoords(const std::vector<int>& comm_coords)
{
 return get_rank_from_comm_coords(comm_coords);
}

   //  The comm_map stores information about which rank stores which part of the lattice.
   //  All mpi ranks have the information.  The comm_coords is encoded to a linear index
   //  which is used as the key in the map.

void LayoutInfo::set_up_comm_map()
{
 vector<int> comm_codes(num_ranks,0);
 vector<int> gl_comm_codes(num_ranks);
 comm_codes[my_rank]=rank_encode(my_comm_coords);
#ifdef ARCH_SERIAL
 comm_map[comm_codes[my_rank]]=my_rank;
#elif defined(ARCH_PARALLEL)
 int status=MPI_Allreduce(comm_codes.data(),gl_comm_codes.data(),num_ranks,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
 if (status!=MPI_SUCCESS){
    errorLaph("Problem occurred in setting up the comm_map in LayoutInfo");}
 for (int rank=0;rank<num_ranks;++rank){
    comm_map[gl_comm_codes[rank]]=rank;}
#endif
   // do a sanity check
 vector<int> check;
 vector<int> comm_coords(Ndim);
 for (comm_coords[0]=0;comm_coords[0]<num_partitions[0];++comm_coords[0])
 for (comm_coords[1]=0;comm_coords[1]<num_partitions[1];++comm_coords[1])
 for (comm_coords[2]=0;comm_coords[2]<num_partitions[2];++comm_coords[2])
 for (comm_coords[3]=0;comm_coords[3]<num_partitions[3];++comm_coords[3]){
    check.push_back(comm_map[rank_encode(comm_coords)]);}
 std::sort (check.begin(),check.end());
 for (int rank=0;rank<num_ranks;++rank){
    if (check[rank]!=rank){
        throw(std::runtime_error("Sanity check did not pass in setting up comm_map in LayoutInfo"));}}
   // output the comm map
 Verbosity v(getVerbosity());
 if (v.isMediumOrHigher()){
    for (comm_coords[0]=0;comm_coords[0]<num_partitions[0];++comm_coords[0])
    for (comm_coords[1]=0;comm_coords[1]<num_partitions[1];++comm_coords[1])
    for (comm_coords[2]=0;comm_coords[2]<num_partitions[2];++comm_coords[2])
    for (comm_coords[3]=0;comm_coords[3]<num_partitions[3];++comm_coords[3]){
       int rank=comm_map[rank_encode(comm_coords)];
       printLaph(make_strf("Sublattice at comm coord (%3d, %3d, %3d, %3d) is on rank %d",
            comm_coords[0],comm_coords[1],comm_coords[2],comm_coords[3],rank));}
    printLaph("\n\n");}
}

  
int LayoutInfo::get_rank_from_comm_coords(const std::vector<int>& comm_coords)
{
 return comm_map[rank_encode(comm_coords)];
}


int LayoutInfo::rank_encode(const std::vector<int>& comm_coords)
{
 return comm_coords[0]+num_partitions[0]*(comm_coords[1]+num_partitions[1]*(comm_coords[2]
        +num_partitions[2]*comm_coords[3]));
}


    // The input lattice quantity in "src" is assumed to be in lexicographic site order,
    // then "dest" is output in even-odd site order.

void LayoutInfo::lexico_to_evenodd(char* dest, const char* src, int sitebytes)
{
 int locNX = LayoutInfo::getRankLattExtents()[0];
 int locNY = LayoutInfo::getRankLattExtents()[1];
 int locNZ = LayoutInfo::getRankLattExtents()[2];
 int locNT = LayoutInfo::getRankLattExtents()[3];
 const char* lex=src;
 char* even=dest;
 char* odd=dest+((locNX*locNY*locNZ*locNT)/2)*sitebytes;
 for (int t=0;t<locNT;++t)
 for (int z=0;z<locNZ;++z)
 for (int y=0;y<locNY;++y){
    bool oddparity=(((t+z+y+my_start_parity)%2)==1);
    for (int x=0;x<locNX;++x,oddparity=!oddparity){
       if (oddparity){
          std::memcpy(odd,lex,sitebytes); odd+=sitebytes;}
       else{
          std::memcpy(even,lex,sitebytes); even+=sitebytes;}
       lex+=sitebytes;}}
}

    // The input lattice quantity in "src" is assumed to be in even-odd site order,
    // then "dest" is output in  lexicographic site order.

void LayoutInfo::evenodd_to_lexico(char* dest, const char* src, int sitebytes)
{
 int locNX = LayoutInfo::getRankLattExtents()[0];
 int locNY = LayoutInfo::getRankLattExtents()[1];
 int locNZ = LayoutInfo::getRankLattExtents()[2];
 int locNT = LayoutInfo::getRankLattExtents()[3];
 char* lex=dest;
 const char* even=src;
 const char* odd=src+((locNX*locNY*locNZ*locNT)/2)*sitebytes;
 for (int t=0;t<locNT;++t)
 for (int z=0;z<locNZ;++z)
 for (int y=0;y<locNY;++y){
    bool oddparity=(((t+z+y+my_start_parity)%2)==1);
    for (int x=0;x<locNX;++x,oddparity=!oddparity){
       if (oddparity){
          std::memcpy(lex,odd,sitebytes); odd+=sitebytes;}
       else{
          std::memcpy(lex,even,sitebytes); even+=sitebytes;}
       lex+=sitebytes;}}
}


int LayoutInfo::linearSiteIndex_lexico(const vector<int>& loc_coords, const vector<int>& loc_extents)
{
 return loc_coords[0]+loc_extents[0]*(loc_coords[1]+loc_extents[1]*(loc_coords[2]
        +loc_extents[2]*loc_coords[3]));
}

     // Assumes only that number of sites on each MPI process is even

int LayoutInfo::linearSiteIndex_evenodd(const vector<int>& loc_coords, const vector<int>& loc_extents,
                                        int loc_nsites, int start_parity)
{
 return (loc_coords[0]+loc_extents[0]*(loc_coords[1]+loc_extents[1]*(loc_coords[2]
        +loc_extents[2]*loc_coords[3])))/2 + loc_nsites/2*((start_parity+loc_coords[0]
        +loc_coords[1]+loc_coords[2]+loc_coords[3])%2);
}


#if defined(ARCH_PARALLEL)

    // This routine uses MPI_Get_processor_name to determine how the MPI ranks
    // have been distributed among the processor nodes.  A single char is returned:  
    // 'B' for block (consecutive ranks on the same node, until node is filled), 
    // 'C' for cyclic (consecutive ranks on consecutive nodes in a round-robin 
    // fashion), 'O' for other.

char LayoutInfo::getMPIRankDistributionOnNodes()
{
 int num_ranks;
 MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
 vector<std::string> node_names(num_ranks);
 int rank;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 char name[MPI_MAX_PROCESSOR_NAME];
 int resultlength;
 MPI_Get_processor_name(name, &resultlength);  // name does not include '/0' terminating null
 if (rank==0){
    node_names[0]=name;}
 for (int otherrank=1;otherrank<num_ranks;++otherrank){
    if (rank==otherrank){
        MPI_Send(name,resultlength,MPI_CHAR,0,otherrank,MPI_COMM_WORLD);}
    else if (rank==0){
        MPI_Status status;
        MPI_Recv(name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,otherrank,otherrank,MPI_COMM_WORLD,&status);
        int count;
        MPI_Get_count(&status,MPI_CHAR,&count);
        node_names[otherrank]=std::string(name,count);}}
   // now that rank 0 has all of the processor names, it can figure out the
   // distribution, then broadcast it to all ranks
 char result='O';
 if (rank==0){
    set<string> uniqnames(node_names.begin(),node_names.end());
    int num_nodes=uniqnames.size();
    int ranks_per_node=num_ranks/num_nodes;
    if (ranks_per_node*num_nodes==num_ranks){
           // check if a blocking
       if (check_pattern(node_names,num_nodes,ranks_per_node,1,num_nodes)){
          result='B';}
       else if (check_pattern(node_names,ranks_per_node,num_nodes,num_nodes,num_nodes)){
          result='C';}}}
    // broadcast result from rank 0 to all other ranks
 MPI_Bcast(&result,1,MPI_CHAR,0,MPI_COMM_WORLD);
 return result;
}

  //  This routine looks at the strings in the vector "names" which must
  //  have size "nblocks"*"blocksize".  It partitions "names" into "nblock" blocks,
  //  each of size "blocksize" as ordered in the vector.  The first "blockdize"
  //  elements of "names" are put into block 1, the next "blocksize" elements
  //  are put into block 2, and so on.  The routine then checks to
  //  see if each block has "ndistinctinblock" distinct strings.  Lastly,
  //  it checks to see if the total number of distinct strings is
  //  "ndistinctall".  Return true if all these checks pass, false otherwise.

bool LayoutInfo::check_pattern(const std::vector<std::string>& names,
                               int nblocks, int blocksize, int ndistinctinblock,
                               int ndistinctall)
{
 if (int(names.size())!=(nblocks*blocksize)) return false;
 set<string> all;
 int count=0;
 for (int b=0;b<nblocks;++b){
    set<string> block;
    for (int r=0;r<blocksize;++r,++count){
       block.insert(names[count]);}
    if (int(block.size())!=ndistinctinblock){
       return false;}
    all.insert(block.begin(),block.end());}
 return (int(all.size())==ndistinctall);
}


//  The optimal way to distribute the sublattices on the MPI ranks depends on 
//  the choice of MPI rank distribution.  quda_laph will try to distribute the 
//  sublattices such that different sublattices in TIME will be on different 
//  nodes.  In other words, the sublattices assigned to MPI ranks on the same 
//  node will differ only in SPACE, allowing increased efficiencies using 
//  peer-to-peer communication for any calculations that do not depend on
//  inter-time information, such as the LapH eigenvector computations and the
//  hadron sinks. For quark sink computations which are fully four-dimensional, 
//  this layout does not matter. 

//  Typically, mpirun or srun distributes the MPI ranks either in a "cyclic" 
//  or "block" pattern. "block" means consecutive ranks are put on the same node 
//  (but different cores), until the node is filled; whereas "cyclic" means that
//  consecutive ranks are placed on consecutive nodes in a round-robin fashion.
//  For a "block" rank-node pattern, then communication coords (x,y,z,t) should 
//  be ordered in a sequence for which the "x" component varies the fastest.
//  For a "cyclic" rank-node pattern, then communication coords (x,y,z,t) should
//  be ordered in a sequence for which the "t" component varies the fastest.
//  The subroutines below are used for these purposes.  One of these will be
//  passed to quda's initCommsGridQuda.

int LayoutInfo::lex_rank_from_commcoord_t_fastest(const int *comm_coord, void *partitions)
{
 int *dims = reinterpret_cast<int *>(partitions);
 int rank = comm_coord[0];
 for (int i = 1; i < 4; i++) { rank = dims[i] * rank + comm_coord[i]; }
 return rank;
}

int LayoutInfo::lex_rank_from_commcoord_x_fastest(const int *comm_coord, void *partitions)
{
 int *dims = reinterpret_cast<int *>(partitions);
 int rank = comm_coord[3];
 for (int i = 2; i >= 0; i--) { rank = dims[i] * rank + comm_coord[i]; }
 return rank;
}

#endif

// **********************************************************************
}
