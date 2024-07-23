#ifdef TESTING
#include <string>
#include "task_tests.h"
#include "xml_handler.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include "utils.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

using namespace std;
using namespace LaphEnv;

namespace QLTestEnv {

// ************************************************

void get_linear_index_mapping(const vector<int>& local_sizes, int start_parity, 
                              map<vector<int>,int>& index_map)
{
 int evencount=0;
 vector<int> site(4);
 int localvol=local_sizes[0]*local_sizes[1]*local_sizes[2]*local_sizes[3];
 int oddcount=localvol/2;
 for (site[3]=0;site[3]<local_sizes[3];++site[3])
 for (site[2]=0;site[2]<local_sizes[2];++site[2])
 for (site[1]=0;site[1]<local_sizes[1];++site[1])
 for (site[0]=0;site[0]<local_sizes[0];++site[0]){
    int parity=(site[0]+site[1]+site[2]+site[3]+start_parity)%2;
    if (parity==0){
       index_map.insert(make_pair(site,evencount));
       evencount++;}
    else{
       index_map.insert(make_pair(site,oddcount));
       oddcount++;}}
}


void testLayoutInfo(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestLayoutInfo")==0)
 return;
 
 printLaph("Running TestLayoutInfo\n");

 map<vector<int>,int> index_map_even, index_map_odd;
 get_linear_index_mapping(LayoutInfo::getRankLattSizes(),0,index_map_even);
 get_linear_index_mapping(LayoutInfo::getRankLattSizes(),1,index_map_odd);
 
 bool flag=true;

 for (int rank=0;rank<LayoutInfo::getNumRanks();++rank){
    if (LayoutInfo::getMyRank()==rank){
       printf("Hello, this is rank %d\n",rank);
       {const vector<int>& tmp=LayoutInfo::getLattSizes();
       printf("Lattice XYZT sizes: %d %d %d %d\n",tmp[0],tmp[1],tmp[2],tmp[3]);}
       {const vector<int>& tmp=LayoutInfo::getCommNumPartitions();
       printf("XYZTCommNumPartitions: %d %d %d %d\n",tmp[0],tmp[1],tmp[2],tmp[3]);}
       {const vector<int>& tmp=LayoutInfo::getRankLattSizes();
       printf("Local (rank) lattice XYZT sizes: %d %d %d %d\n",tmp[0],tmp[1],tmp[2],tmp[3]);}
       printf("Total number of lattice sites is %d\n",LayoutInfo::getLatticeNumSites());
       printf("Number of lattice sizes on each rank is %d\n",LayoutInfo::getRankLatticeNumSites());
       printf("Number of space-time dimensions is %d\n",LayoutInfo::Ndim);
       {const vector<int>& tmp=LayoutInfo::getMyCommCoords() ;
       printf("My rank coords: %d %d %d %d\n",tmp[0],tmp[1],tmp[2],tmp[3]);}
       printf("My rank is %d\n\n",LayoutInfo::getMyRank());
       for (int k0=0;k0<LayoutInfo::getCommNumPartitions()[0];++k0)
       for (int k1=0;k1<LayoutInfo::getCommNumPartitions()[1];++k1)
       for (int k2=0;k2<LayoutInfo::getCommNumPartitions()[2];++k2)
       for (int k3=0;k3<LayoutInfo::getCommNumPartitions()[3];++k3){
           vector<int> latt_coords(LayoutInfo::Ndim);
           latt_coords[0]=k0*LayoutInfo::getRankLattSizes()[0];
           latt_coords[1]=k1*LayoutInfo::getRankLattSizes()[1];
           latt_coords[2]=k2*LayoutInfo::getRankLattSizes()[2];
           latt_coords[3]=k3*LayoutInfo::getRankLattSizes()[3];
           vector<int> ccheck=LayoutInfo::getCommCoordsFromLatticeCoords(latt_coords);
           if ((ccheck[0]!=k0)||(ccheck[1]!=k1)||(ccheck[2]!=k2)||(ccheck[3]!=k3)){
              printf("check of getCommCoordsFromLatticeCoords failed on this rank");}
           int checkrank=LayoutInfo::getRankFromLatticeCoords(latt_coords);
           printf("This rank thinks rank for lattice coords (%d,%d,%d,%d) is %d\n",
                  latt_coords[0],latt_coords[1],latt_coords[2],latt_coords[3],checkrank);}
       fflush(stdout);}}
#ifdef ARCH_PARALLEL
    comm_barrier();
#endif
 for (int y=0;y<LayoutInfo::getLattSizes()[1];++y)
 for (int x=0;x<LayoutInfo::getLattSizes()[0];++x)
 for (int t=0;t<LayoutInfo::getLattSizes()[3];++t)
 for (int z=0;z<LayoutInfo::getLattSizes()[2];++z){
    vector<int> latt_coords(LayoutInfo::Ndim);
    latt_coords[0]=x; latt_coords[1]=y; latt_coords[2]=z; latt_coords[3]=t;
    int rank, rank_site_linear_index;
    LayoutInfo::getCommInfoFromLatticeCoords(latt_coords,rank,rank_site_linear_index);
    vector<int> comm_coord(LayoutInfo::Ndim);
    vector<int> rank_site(LayoutInfo::Ndim);
    vector<int> start_coord(LayoutInfo::Ndim);
    int start_parity=0;
    for (int dir=0;dir<LayoutInfo::Ndim;++dir){
       comm_coord[dir]=latt_coords[dir]/LayoutInfo::getRankLattSizes()[dir];
       start_coord[dir]=comm_coord[dir]*LayoutInfo::getRankLattSizes()[dir];
       start_parity+=start_coord[dir];
       rank_site[dir]=latt_coords[dir]-start_coord[dir];}
    start_parity=start_parity%2;
    int index=(start_parity)? index_map_odd[rank_site] : index_map_even[rank_site];
    string status;
    if (rank_site_linear_index!=index){ flag=false; status="WRONG";}
    else{ status="correct";}
    printLaph(make_strf("For lattice site (%d,%d,%d,%d)  rank = %d  linear_index = %d should be %d %s",
               latt_coords[0],latt_coords[1],latt_coords[2],latt_coords[3],rank,rank_site_linear_index,index,status));
    }

 uint lexico_lin_index=0;
 for (int t=0;t<LayoutInfo::getLattSizes()[3];++t)
 for (int z=0;z<LayoutInfo::getLattSizes()[2];++z)
 for (int y=0;y<LayoutInfo::getLattSizes()[1];++y)
 for (int x=0;x<LayoutInfo::getLattSizes()[0];++x,++lexico_lin_index){
    int rank=LayoutInfo::getRankFromLexicoLinearIndex(lexico_lin_index);
    vector<int> latt_coords(LayoutInfo::Ndim);
    latt_coords[0]=x; latt_coords[1]=y; latt_coords[2]=z; latt_coords[3]=t;
    if (rank!=LayoutInfo::getRankFromLatticeCoords(latt_coords)) flag=false;}

 bool success=LaphEnv::globalAnd(flag);
 if (success) printLaph("ALL TESTS PASSED!!");
 else         printLaph("Some tests FAILED");

}

// ***********************************************
}
#endif
