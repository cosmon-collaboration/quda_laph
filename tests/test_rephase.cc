#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "quark_smearing_handler.h"
#include "field_ops.h"

//#define PRINT_EVECS
//#define ALL_CONSTANT

using namespace std ;
using namespace LaphEnv ;

void
old_code( vector<LattField> &laph_evecs )
{
 int nev=laph_evecs.size();
 if (nev==0){
    return;}
 for (int v=0;v<nev;++v){
    if (laph_evecs[v].getFieldSiteType()!=FieldSiteType::ColorVector){
       errorLaph("Applying phase convention can only be done to ColorVector fields");}}

 bool dp=(laph_evecs[0].bytesPerWord()==sizeof(std::complex<double>));
 int loc_nsites=LayoutInfo::getRankLatticeNumSites();
 int loc_npsites=loc_nsites/2;
 int start_parity=LayoutInfo::getMyStartParity();
 int nloctime=LayoutInfo::getRankLattExtents()[3];
 int tstride=LayoutInfo::getRankLattExtents()[0]*LayoutInfo::getRankLattExtents()[1]
            *LayoutInfo::getRankLattExtents()[2];
 int incx=FieldNcolor;
 int cbytes=(dp)?sizeof(std::complex<double>):sizeof(std::complex<float>);
 int bps=laph_evecs[0].bytesPerSite();
 complex<double> rephasedp;
 complex<float> rephasesp;

        //  get the rephase factors for each eigvec and local time slice
 vector<char> rephase(nev*nloctime*cbytes);
 if ((LayoutInfo::getMyCommCoords()[0]==0)
   &&(LayoutInfo::getMyCommCoords()[1]==0)
   &&(LayoutInfo::getMyCommCoords()[2]==0)){
    char* rp=rephase.data();
    for (int v=0;v<nev;++v){
       char* fp=reinterpret_cast<char*>(laph_evecs[v].getDataPtr());
       for (int tloc=0;tloc<nloctime;++tloc,rp+=cbytes){
          int parshift=loc_npsites*((start_parity+tloc)%2);
          int start1=((tstride*tloc)/2) + parshift;
          char* x1=fp+bps*start1;         // location of (0,0,0,t)
          if (dp){
             complex<double>* zptr=reinterpret_cast<complex<double>*>(x1);
             complex<double> z(std::conj(*zptr));
             double r=std::abs(z);
             if (r<1e-12){
                errorLaph(make_str("problem applying phase convention: 0-th color component at site",
                          " (0,0,0) has very small magnitude ",r));}  // pray this does not happen!
             rephasedp=z/r;
             std::memcpy(rp,&rephasedp,cbytes);}
          else{
             complex<float>* zptr=reinterpret_cast<complex<float>*>(x1);
             complex<float> z(std::conj(*zptr));
             float r=std::abs(z);
             if (r<1e-8){
                errorLaph(make_str("problem applying phase convention: 0-th color component at site",
                          " (0,0,0) has very small magnitude ",r));}  // pray this does not happen!
             rephasesp=z/r;
             std::memcpy(rp,&rephasesp,cbytes);}}}}

#ifdef ARCH_PARALLEL
    // now broadcast these rephase factors to the ranks that need them
 vector<int> comm_coords(LayoutInfo::Ndim-1);
 comm_coords[LayoutInfo::Ndim-1]=LayoutInfo::getMyCommCoords()[LayoutInfo::Ndim-1];
 comm_coords[0]=0; comm_coords[1]=0; comm_coords[2]=0;
 int orig_sender_rank=LayoutInfo::getRankFromCommCoords(comm_coords);
 int sender_rank=0;
 vector<int> broadcast_ranks;
 int count=0;
 for (comm_coords[0]=0;comm_coords[0]<LayoutInfo::getCommNumPartitions()[0];++comm_coords[0])
 for (comm_coords[1]=0;comm_coords[1]<LayoutInfo::getCommNumPartitions()[1];++comm_coords[1])
 for (comm_coords[2]=0;comm_coords[2]<LayoutInfo::getCommNumPartitions()[2];++comm_coords[2],++count){
    int next_rank=LayoutInfo::getRankFromCommCoords(comm_coords);
    if (next_rank==orig_sender_rank){
       sender_rank=count;}
    broadcast_ranks.push_back(next_rank);}
 if (broadcast_ranks.size()>1){
    if (int(broadcast_ranks.size())==LayoutInfo::getNumRanks()){
       comm_broadcast(rephase.data(),nev*nloctime*cbytes,orig_sender_rank);}
    else{
       MPI_Group world_group;
       MPI_Comm_group(MPI_COMM_WORLD, &world_group);
       MPI_Group time_group;
       MPI_Group_incl(world_group, broadcast_ranks.size(), broadcast_ranks.data(), &time_group);
       MPI_Comm time_comm;
       MPI_Comm_create_group(MPI_COMM_WORLD, time_group, sender_rank, &time_comm);
       if (time_comm != MPI_COMM_NULL){
          int status=MPI_Bcast(rephase.data(), nev*nloctime*cbytes, MPI_BYTE, sender_rank, time_comm);
          if (status!=MPI_SUCCESS){
             errorLaph("Broadcast of re-phase factor in applyLaphPhaseConvention failed");}}
       MPI_Group_free(&time_group);
       MPI_Group_free(&world_group);
       MPI_Comm_free(&time_comm);}}
#endif

    // now apply the rephase factors
 const char* rf=rephase.data();
 for (int v=0;v<nev;++v){
    char* fp=reinterpret_cast<char*>(laph_evecs[v].getDataPtr());
    for (int tloc=0;tloc<nloctime;++tloc,rf+=cbytes){
       int parshift=loc_npsites*((start_parity+tloc)%2);
       int start1=((tstride*tloc)/2) + parshift;
       int stop1=((1+tstride*(tloc+1))/2) + parshift;
       int n1=stop1-start1;
       parshift=loc_npsites*((start_parity+1+tloc)%2);
       int start2=((1+tstride*tloc)/2) + parshift;
       int stop2=((tstride*(tloc+1))/2) + parshift;
       int n2=stop2-start2;
       char* x1=fp+bps*start1;         // location of (0,0,0,t)
       char* x2=fp+bps*start2;         // location of (1,0,0,t)
       for (int c=0;c<FieldNcolor;++c){
          if (dp){
             cblas_zscal(n1,rf,x1,incx);
             cblas_zscal(n2,rf,x2,incx);}
          else{
             cblas_cscal(n1,rf,x1,incx);
             cblas_cscal(n2,rf,x2,incx);}
          x1+=cbytes; x2+=cbytes;}
       }}
#ifdef ARCH_PARALLEL
 comm_barrier();
 MPI_Barrier( MPI_COMM_WORLD ) ;
#endif
}

static void
printevecs( vector<LattField> laphEigvecs , const int rank )
{
#ifdef PRINT_EVECS
  int myrank = rank ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank(  MPI_COMM_WORLD  , &myrank ) ;
#endif
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    if( myrank == rank ) {
      std::cout<<"Ev "<<n<<" | rank "<<rank<<std::endl ;
      complex<double> *pt = (complex<double>*)laphEigvecs[n].getDataPtr() ;
      for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
	for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	  std::cout<<i<<"/"<<c<<" "<<pt[c+FieldNcolor*i]<<std::endl ;
	}
      }
      std::cout<<std::endl ;
      fflush( stdout ) ;
    }
  }
#endif
}

static void
set_constant( vector<LattField> &laphEigvecs )
{
  int myrank = 0 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank ) ;
#endif
  // give them some bullshit values
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    #ifdef ALL_CONSTANT
    complex<double> z( (n+1) , (n+1) ) ;
    setConstantField( laphEigvecs[n], z );
    #else
    complex<double> *ptr = (complex<double>*)laphEigvecs[n].getDataPtr() ;
    const size_t V = (size_t)LayoutInfo::getRankLatticeNumSites() ;
    for( size_t i = 0 ; i < V ; i++ ) {
      for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	const complex<double> z( n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ,
				 n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ) ;
	*ptr = z ;
	ptr++ ;
      }
    }
    #endif
  }
}

int main(int argc, char *argv[]) {
  XMLHandler xml_in;

  if( init_quda_laph(argc, argv, xml_in) != 0 ) {
    exit(1) ;
  }

  int global = 1 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_size( MPI_COMM_WORLD , &global ) ;
#endif
  
  // call rephase here
  const int Nev = 8 ;
  vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);

  cout<<"Constant Eigvecs"<<std::endl ;
  set_constant( laphEigvecs ) ;
  for( int i = 0 ; i < global ; i++ ) {
    printevecs( laphEigvecs , i ) ;
  }
  
  // why is everything a fucking "Handler"
  QuarkSmearingHandler Handle ;
  StopWatch timer ;
  timer.start() ;
  Handle.applyLaphPhaseConvention( laphEigvecs ) ;
  timer.stop() ;
  printLaph(make_strf("\nNew code took = %g seconds\n",
                      timer.getTimeInSeconds()));
  
  
  cout<<"Rephased Eigvecs"<<std::endl ;
  for( int i = 0 ; i < global ; i++ ) {
    printevecs( laphEigvecs , i ) ;
  }

  cout<<"Old code"<<std::endl ;
  vector<LattField> laphEigvecsComp( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecsComp ) ;
  timer.reset() ;
  timer.start() ;
  old_code( laphEigvecsComp ) ;
  timer.stop();
  printLaph(make_strf("\nOld code took = %g seconds\n",
                      timer.getTimeInSeconds()));

  for( int i = 0 ; i < global ; i++ ) {
    printevecs( laphEigvecsComp , i ) ;
  }

  // compute the difference between them
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    double loc_dev = 0 ;
    complex<double> *pt1 = (complex<double>*)laphEigvecs[n].getDataPtr() ;
    complex<double> *pt2 = (complex<double>*)laphEigvecsComp[n].getDataPtr() ;
    for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
      for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	loc_dev += abs( *pt1 - *pt2 ) ;
	pt1++ ; pt2++ ;
      }
    }
    double deviation = 0 ;
    #ifdef ARCH_PARALLEL
    MPI_Allreduce(&loc_dev, &deviation, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
    #else
    deviation = loc_dev ;
    #endif
    int my_rank = 0 ;
    #ifdef ARCH_PARALLEL
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank ) ;
    #endif
    if( my_rank == 0 ) {
      std::cout<<n<<" Deviation New vs Old "<<scientific<<deviation<<std::endl ;
    }
  }
  
  finalize();

  return 0;
}
