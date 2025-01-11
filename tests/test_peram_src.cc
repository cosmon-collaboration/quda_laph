#include "QudaLaphIncludes.h"
#include "init_quda_laph.h"
#include "perambulator_handler.h"
#include "field_ops.h"

#define PRINT_EVECS
//#define ALL_CONSTANT

typedef std::complex<double> dcmplx;
typedef std::complex<float> fcmplx;

using namespace LaphEnv ;

static void
printevecs( std::vector<LattField> laphEigvecs , const int rank )
{
#ifdef PRINT_EVECS
  int myrank = rank ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank(  MPI_COMM_WORLD  , &myrank ) ;
#endif
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    if( myrank == rank ) {
      std::cout<<"Ev "<<n<<" | rank "<<rank<<std::endl ;
      std::complex<double> *pt = (std::complex<double>*)laphEigvecs[n].getDataPtr() ;
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
printferm( LattField &ferm , const int rank )
{
#ifdef PRINT_EVECS
  int myrank = rank ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank(  MPI_COMM_WORLD  , &myrank ) ;
#endif
  if( myrank == rank ) {
    std::complex<double> *pt = reinterpret_cast<std::complex<double>*>(ferm.getDataPtr()) ;
    for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
      for( size_t c = 0 ; c < (size_t)12 ; c++ ) {
	std::cout<<i<<"/"<<c<<" "<<pt[c+12*i]<<std::endl ;
      }
    }
    std::cout<<std::endl ;
    fflush( stdout ) ;
  }
#endif
}

// Creates the source (multiplied by gamma_4) in Dirac-Pauli basis
//   src_spin here is 1,2,3,4

void
old_code(LattField &ferm_src,
	 const void *ev_src_ptr,
	 const int src_time,
	 const int src_spin) {
  printLaph(" Making source for this inversion...");
  ferm_src.reset(FieldSiteType::ColorSpinVector);
  const bool dp = (ferm_src.bytesPerWord() == sizeof(std::complex<double>));
  // initialize source field to zero
  const int loc_nsites = LayoutInfo::getRankLatticeNumSites();
  const int ncmplx_per_site = ferm_src.elemsPerSite();
  const int ncmplx = ncmplx_per_site * loc_nsites;
  int cbytes;
  dcmplx zrhodp;
  fcmplx zrhosp;
  char *zrho;
  if (dp) {
    double *z0 = reinterpret_cast<double *>(ferm_src.getDataPtr());
    std::fill(z0, z0 + 2 * ncmplx, 0.0);
    cbytes = sizeof(std::complex<double>);
    zrho = reinterpret_cast<char *>(&zrhodp);
  } else {
    float *z0 = reinterpret_cast<float *>(ferm_src.getDataPtr());
    std::fill(z0, z0 + 2 * ncmplx, 0.0);
    cbytes = sizeof(std::complex<float>);
    zrho = reinterpret_cast<char *>(&zrhosp);
  }

  const int loc_npsites = loc_nsites / 2;
  const int start_parity = LayoutInfo::getMyStartParity();
  const int mytmin = LayoutInfo::getMyCommCoords()[3] * LayoutInfo::getRankLattExtents()[3];
  const int mytmax = mytmin + LayoutInfo::getRankLattExtents()[3] - 1;
  const int tstride = LayoutInfo::getRankLattExtents()[0] *
                      LayoutInfo::getRankLattExtents()[1] *
                      LayoutInfo::getRankLattExtents()[2];
  const int incx = FieldNcolor;
  const int incy = FieldNcolor * FieldNspin;

  // could be more efficient
  if ((src_time >= mytmin) && (src_time <= mytmax)) {
    int tloc = src_time - mytmin;
    int parshift = loc_npsites * ((start_parity + tloc) % 2);
    int start1 = ((tstride * tloc) / 2) + parshift;
    int stop1 = ((1 + tstride * (tloc + 1)) / 2) + parshift;
    int n1 = stop1 - start1;
    parshift = loc_npsites * ((start_parity + 1 + tloc) % 2);
    int start2 = ((1 + tstride * tloc) / 2) + parshift;
    int stop2 = ((tstride * (tloc + 1)) / 2) + parshift;
    int n2 = stop2 - start2;
    int xstart1 = start1 * incx * cbytes;
    int xstart2 = start2 * incx * cbytes;
    char *ystart1 = ferm_src.getDataPtr() + start1 * incy * cbytes;
    char *ystart2 = ferm_src.getDataPtr() + start2 * incy * cbytes;

    const char *x0 = reinterpret_cast<const char *>(ev_src_ptr);
    zrhodp = dcmplx(1.0, 0.0);
    if (src_spin > 2) {
      zrhodp = -zrhodp;
    } // multiply by gamma_4
    if (!dp) {
      zrhosp = std::complex<float>(real(zrhodp), imag(zrhodp));
    }
    const char *x1 = x0 + xstart1;
    const char *x2 = x0 + xstart2;
    char *y1 = ystart1 + (src_spin - 1) * incx * cbytes;
    char *y2 = ystart2 + (src_spin - 1) * incx * cbytes;
    for (int c = 0; c < FieldNcolor; ++c) {
      if (dp) {
        cblas_zaxpy(n1, (dcmplx *)(zrho), (dcmplx *)(x1), incx, (dcmplx *)(y1),
                    incy);
        cblas_zaxpy(n2, (dcmplx *)(zrho), (dcmplx *)(x2), incx, (dcmplx *)(y2),
                    incy);
      } else {
        cblas_caxpy(n1, (fcmplx *)(zrho), (fcmplx *)(x1), incx, (fcmplx *)(y1),
                    incy);
        cblas_caxpy(n2, (fcmplx *)(zrho), (fcmplx *)(x2), incx, (fcmplx *)(y2),
                    incy);
      }
      x1 += cbytes;
      y1 += cbytes;
      x2 += cbytes;
      y2 += cbytes;
    }
  }
  printLaph("Source for this inversion created");
}


static void
set_constant( std::vector<LattField> &laphEigvecs )
{
  int myrank = 0 ;
#ifdef ARCH_PARALLEL
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank ) ;
#endif
  // give them some bullshit values
  for( size_t n = 0 ; n < laphEigvecs.size() ; n++ ) {
    #ifdef ALL_CONSTANT
    std::complex<double> z( (n+1) , (n+1) ) ;
    setConstantField( laphEigvecs[n], z );
    #else
    std::complex<double> *ptr = (std::complex<double>*)laphEigvecs[n].getDataPtr() ;
    const size_t V = (size_t)LayoutInfo::getRankLatticeNumSites() ;
    for( size_t i = 0 ; i < V ; i++ ) {
      for( size_t c = 0 ; c < (size_t)FieldNcolor ; c++ ) {
	const std::complex<double> z( n+laphEigvecs.size()*(c+FieldNcolor*(i+V*myrank))+1 ,
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
  
  // set eigvecs here
  const int Nev = 2 ;
  std::vector<LattField> laphEigvecs( Nev, FieldSiteType::ColorVector);
  set_constant( laphEigvecs ) ;
  printevecs( laphEigvecs , 0 ) ;
  
  std::cout<<"Evs set"<<std::endl;
  std::vector<void *> evList(Nev);
  for (int n = 0; n < Nev; n++) {
    evList[n] = (void *)(laphEigvecs[n].getDataConstPtr());
  }
  std::cout<<"Evlist made"<<std::endl;
  void **evs_ptr = (void**)evList.data() ;

  for( int mu = 1 ; mu < 5 ; mu++ ) {
    
    LattField ferm_src( FieldSiteType::ColorSpinVector);
    old_code( ferm_src, evs_ptr[0], 0, mu) ;
    printferm( ferm_src , 0 ) ;

    LattField ferm_src_new( FieldSiteType::ColorSpinVector );
    PerambulatorHandler peram = PerambulatorHandler() ;
    peram.make_source( ferm_src_new , evs_ptr[0] , 0, mu) ;
    printferm( ferm_src_new , 0 ) ;
  
    // compute the difference between them
    double loc_dev = 0 ;
    std::complex<double> *pt1 = (std::complex<double>*)ferm_src.getDataPtr() ;
    std::complex<double> *pt2 = (std::complex<double>*)ferm_src_new.getDataPtr() ;
    for( size_t i = 0 ; i < (size_t)LayoutInfo::getRankLatticeNumSites() ; i++ ) {
      for( size_t c = 0 ; c < (size_t)12 ; c++ ) {
	//std::cout<<*pt1<<" "<<*pt2<<" ->"<<*pt1-*pt2<<std::endl ;	
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
      std::cout<<" Deviation New vs Old "<<std::scientific<<deviation<<std::endl ;
    }
  }
  finalize();
  
  return 0;
}
