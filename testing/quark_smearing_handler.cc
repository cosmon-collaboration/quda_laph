#include "quark_smearing_handler.h"
#include "stop_watch.h"
#include "gluon_smearing_handler.h"
#if defined(USE_GSL_CBLAS)
#include "gsl_cblas.h"
#elif defined(USE_OPENBLAS)
#include "openblas.h"
#endif

#include "task_tests.h"
using namespace QLTestEnv;

using namespace std;

namespace LaphEnv {



/*
#ifdef TESTING
#if (QDP_ND == 4)
void printLaphEigenvectors(LaphVectorHandler& laphEigvecs, 
                           const GaugeConfigurationInfo& gaugeinfo,
                           int nEigvecs, TextFileWriter& fout);
#elif (QDP_ND == 3)
void printLaphEigenvectors(LaphVectorHandler& laphEigvecs, int t,
                           const GaugeConfigurationInfo& gaugeinfo,
                           int nEigvecs, TextFileWriter& fout);
#endif
#endif
*/


 // **************************************************************
 // *                                                            *
 // *                                                            *
 // *           QuarkSmearingHandler implementation              *
 // *                                                            *
 // *                                                            *
 // **************************************************************


   // constructors

QuarkSmearingHandler::QuarkSmearingHandler()
     : uPtr(0), gSmearPtr(0), qSmearPtr(0),
       m_read_mode(true), dh_ptr(0)
{
}

QuarkSmearingHandler::QuarkSmearingHandler(const GluonSmearingInfo& gluon_smearing,
                                           const GaugeConfigurationInfo& gauge,
                                           const QuarkSmearingInfo& quark_smearing,
                                           const string& smeared_quark_file_stub,
                                           bool read_mode)
{
 set_info(gluon_smearing,gauge,quark_smearing,smeared_quark_file_stub,read_mode);
}

void QuarkSmearingHandler::setInfo(const GluonSmearingInfo& gluon_smearing,
                                   const GaugeConfigurationInfo& gauge,
                                   const QuarkSmearingInfo& quark_smearing,
                                   const string& smeared_quark_file_stub,
                                   bool read_mode)
{
 clear();
 set_info(gluon_smearing,gauge,quark_smearing,smeared_quark_file_stub,read_mode);
}



void QuarkSmearingHandler::set_info(const GluonSmearingInfo& gluon_smearing,
                                    const GaugeConfigurationInfo& gauge,
                                    const QuarkSmearingInfo& quark_smearing,
                                    const string& smeared_quark_file_stub,
                                    bool read_mode)
{
 smearedQuarkFileStub=tidyString(smeared_quark_file_stub);
 if (smearedQuarkFileStub.empty()){
    errorLaph("empty file name or stub in QuarkSmearingHandler");}

 try{
    uPtr=new GaugeConfigurationInfo(gauge);
    gSmearPtr=new GluonSmearingInfo(gluon_smearing);
    qSmearPtr=new QuarkSmearingInfo(quark_smearing);}
 catch(const std::exception& xp){
    errorLaph("problem in setInfo in QuarkSmearingHandler");}

 m_read_mode=read_mode;

 if (m_read_mode){
/*#if (QDP_ND == 3)
    uint mintime=uPtr->getMinTime();
    uint maxtime=uPtr->getMaxTime();
    FileListInfo files(smearedQuarkFileStub+"_time",
                       mintime,maxtime,false);
    try{ dh_ptr=new DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelKey,
                    LattField>(*this,files,"Laph--SmearedQuarkTimeFile",
                    "LaphEigenvectors");
       bool flag=false;
       for (uint tval=mintime;tval<=maxtime;tval++){
          if (dh_ptr->queryFile(TimeKey(tval))){ flag=true; break;}}
       if (flag)
          QDPIO::cout << "QuarkSmearingHandler file stub appears to be okay"<<endl;
       else{
          QDPIO::cout << "problem with QuarkSmearingHandler file stub: "<<smeared_quark_file_stub<<endl;
          QDP_abort(1);}}
 #elif (QDP_ND == 4)  */

         // read from files (or the NamedObjMap)
    FileListInfo files(smearedQuarkFileStub+"_level",
                       0,qSmearPtr->getNumberOfLaplacianEigenvectors()-1,false);
   // if (!fileExists(files.getFileName(0))){
   //    readTimeSlicesIntoNOM();  // read time slices, reorganize and put into NamedObjMap
   //    FileListInfo nomfiles("NOM_LaphEig_level",
   //                          0,qSmearPtr->getNumberOfLaplacianEigenvectors()-1,false);
   //    files=nomfiles;}
    try{ dh_ptr=new DataGetHandlerMFO<QuarkSmearingHandler,LevelKey,LevelKey,
                        LattField>(*this,files,"Laph--SmearedQuarkLevelFile",
                        "LaphEigenvectors");
       if (dh_ptr->queryFile(LevelKey(0)))
          printLaph("QuarkSmearingHandler file stub appears to be okay");
       else{
          errorLaph(make_strf("problem with QuarkSmearingHandler file stub: %s",
                    smeared_quark_file_stub));}}
    catch(const std::exception& xp){ 
      errorLaph("unable to allocate data handler in QuarkSmearingHandler");}}
 else{
    dh_ptr=0;}
//#if (QDP_ND == 3)
// dph_ptr=0;
//#endif */
}

    // use this to increase the number of Laplacian eigenvectors allowed;
    // need inside QuarkHandler since smearing sub-handler is static.

void QuarkSmearingHandler::updateSmearing(const QuarkSmearingInfo& quark_smearing)
{
 check_info_set("updateQuarkSmearingInfo");
 qSmearPtr->increaseUpdate(quark_smearing);
}


QuarkSmearingHandler::~QuarkSmearingHandler()
{
 clear();
}


void QuarkSmearingHandler::clear()
{
 try{
    delete gSmearPtr;
    delete qSmearPtr;
    delete uPtr;
    delete dh_ptr;} 
 catch(const std::exception& xp) {errorLaph("abort");}
 smearedQuarkFileStub.clear();
 gSmearPtr=0;
 qSmearPtr=0;
 uPtr=0;
// dh_ptr=0;
//#if (QDP_ND == 3)
// Eigenvalues.resize(0);
// delete dph_ptr;
// dph_ptr=0;
//#elif (QDP_ND == 4)
 Eigenvalues.clear();
//#endif
}


           // access to the info

bool QuarkSmearingHandler::isInfoSet() const
{
 return ((gSmearPtr!=0)&&(uPtr!=0)&&(qSmearPtr!=0)
        &&(!smearedQuarkFileStub.empty()));
} 


   // check_mode = 0 means no check, 1 means check for read mode,
   // 2 means check for write mode

void QuarkSmearingHandler::check_info_set(const string& name, int check_mode) const
{
 if (!isInfoSet()){
    errorLaph(make_strf("error in QuarkSmearingHandler: must setInfo before calling %s",name));}
 if (check_mode==1){
    if (!m_read_mode){
       errorLaph(make_strf("error in QuarkSmearingHandler: must be in read mode when calling %s",name));}}
 else if (check_mode==2){
    if (m_read_mode){
       errorLaph(make_strf("error in QuarkSmearingHandler: must not be in read mode when calling %s",name));}}
}

const GluonSmearingInfo& QuarkSmearingHandler::getGluonSmearingInfo() const
{
 check_info_set("getGluonSmearingInfo");
 return *gSmearPtr;
}

const QuarkSmearingInfo& QuarkSmearingHandler::getQuarkSmearingInfo() const
{
 check_info_set("getQuarkSmearingInfo");
 return *qSmearPtr;
}

const GaugeConfigurationInfo& QuarkSmearingHandler::getGaugeConfigurationInfo() const
{
 check_info_set("getGaugeConfigurationInfo");
 return *uPtr;
}

const string& QuarkSmearingHandler::getSmearedQuarkFieldFileStub() const
{
 check_info_set("getSmearedQuarkFieldFileStub");
 return smearedQuarkFileStub;
}
/*
void QuarkSmearingHandler::getFileMap(XMLHandler& xmlout) const
{
 check_info_set("getFileMap",1);
 dh_ptr->getFileMap(xmlout);
}

void QuarkSmearingHandler::outputKeys(XMLHandler& xmlout)
{
 check_info_set("outputKeys",1);
 dh_ptr->outputKeys(xmlout);
}
*/

void QuarkSmearingHandler::failure(const string& message)
{
 errorLaph(make_strf("%s in QuarkSmearingHandler",message));
}


// *******************************************************************************

/*
void QuarkSmearingHandler::LevelDispKey::encode(uint in_val, const DirPath& displace)
{
 if ((in_val>65535)||(displace.Length()>5)){
    QDPIO::cerr << "Unsupported initialization of QuarkSmearingHandler::LevelDispKey"<<std::endl; 
    QDP_abort(1);}
 value=in_val;
 uint path=0;
 for (DirPath::const_iterator seg=displace.begin();seg!=displace.end();++seg){
    path<<=3;
    path|=(*seg)+4;}
 value<<=16;
 value|=path;
}

QuarkSmearingHandler::LevelDispKey::LevelDispKey(uint in_val, const DirPath& displace)
{
 encode(in_val,displace);
}

QuarkSmearingHandler::LevelDispKey::LevelDispKey(XMLHandler& xmlr)
{
 XMLHandler xmlrr(xmlr,"./descendant-or-self::LevelDispKey");
 uint level;
 xmlread(xmlrr,"Level",level,"QuarkSmearingHandler::LevelDispKey");
 DirPath displace(xmlrr);
 encode(level,displace);
}

uint QuarkSmearingHandler::LevelDispKey::getLevel() const
{
 uint level=value>>16;
 return level;
}

DirPath QuarkSmearingHandler::LevelDispKey::getDisplacement() const
{
 uint dispcode=value & 0xFFFFu;
 vector<int> dirs;
 uint dir=dispcode & 0x7u;
 while (dir!=0){
    dirs.push_back(dir-4);
    dispcode>>=3;
    dir=dispcode & 0x7u;}
 DirPath displace;
 for (int k=dirs.size()-1;k>=0;k--){
    displace.addDir(dirs[k]);}
 return displace;
}

void QuarkSmearingHandler::LevelDispKey::output(XMLHandler& xmlw) const 
{
 DirPath displace(getDisplacement());
 uint level=getLevel();
 push(xmlw,"LevelDispKey");
 write(xmlw,"Level",level);
 push(xmlw,"DirPath");
 multi1d<int> dirs(displace.Length());
 int k=0;
 for (DirPath::const_iterator it=displace.begin();it!=displace.end();it++){
     dirs[k++]=*it;}
 write(xmlw,"Directions",dirs);
 pop(xmlw);
 pop(xmlw);
}
*/
 // **********************************************************
/*
#if (QDP_ND == 3)

           // Compute the Laph Eigenvectors and write to file or NamedObjMap.
           // Will not overwrite; user must do a "clear" first.

void QuarkSmearingHandler::computeLaphEigenvectors(
                         const LaphEigenSolverInfo& solver_info,
                         const string& smeared_gauge_file)
{
 check_info_set("computeLaphEigenvectors",2);

 StopWatch rolex,bulova;
 rolex.start();
 double iotime=0.0;

 int nTime = uPtr->getTimeExtent();
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
// int firsttime=uPtr->getMinTime();
// int lasttime=uPtr->getMaxTime();
// if ((mintime>firsttime)&&(mintime<=lasttime)) firsttime=mintime;
// if ((maxtime<lasttime)&&(maxtime>=firsttime)) lasttime=maxtime;

    // fire up the GluonSmearingHandler

 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);

 printLaph(make_str("Computation of Laplacian eigenvectors commencing",
                    " in FieldSmearingHandler"));
 printLaph(make_str("  number of requested eigenvectors = ",nEigvecs));
 printLaph(make_str("            time extent of lattice = ",nTime));

//#ifdef TESTING
// TextFileWriter fout("lapheigvec3d.log");
//#endif

    // fire up the data handler

 FileListInfo writeFiles(smearedQuarkFileStub+"_time",firsttime,lasttime,false);
 DataPutHandlerMFO<QuarkSmearingHandler,TimeKey,LevelKey,
       LattField> DHput(*this,writeFiles,"Laph--SmearedQuarkTimeFile",
                                 "LaphEigenvectors",true,false,
                                 striping_factor,striping_unit);
 
     // do the computation 
  eigensolveQuda(host_evecs_ptr.data(), evals.data(), &eig_param);

void eigensolveQuda(void **h_evecs, double_complex *h_evals, QudaEigParam *param);


 for (int t=firsttime; t<=lasttime; t++){

    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    if (fileExists(oss.str())){
       QDPIO::cout << "file "<<oss.str()<<" already exists"<<endl;
       QDPIO::cout << " will not overwrite"<<endl<<endl;
       continue;}

    QDPIO::cout <<endl<<endl<<"Beginning computation for time = "
                <<t<<" (with "<<lasttime<<" as last)"<<endl<<endl;
    StopWatch iotimer; iotimer.start();
    const multi1d<LatticeColorMatrix>& usmear
                        = gHandler.getSmearedGaugeFieldTimeSlice(t);
    iotimer.stop();
    iotime+=iotimer.getTimeInSeconds();

         // use the Krylov-spectral restarted Lanczos method to compute
         // the eigenvalues/eigenvectors of the Laplacian on this time slice

    LaphVectorHandler laphEigvecs(usmear);
    int chebyshevOrder = solver_info.getChebyshevOrder();
    if (chebyshevOrder>=2){
       laphEigvecs.setChebyshevAccelerationOn(
                      solver_info.getMaximumEigenvalue(),
                      solver_info.getCutoffEigenvalue(),chebyshevOrder);}

    bulova.reset();bulova.start();
    KSRLanczosEigenSolver Comp(laphEigvecs, nEigvecs,
                  solver_info.getResidualTolerance(),
                  solver_info.getKrylovDimension(),
                  solver_info.getMaximumIterations(),
                  laphEigvecs.getSpectrumEnd(),
                  solver_info.getOutputVerbosity());

    if (!Comp.wasSuccessful()){
       QDPIO::cout << endl<<endl<<"Convergence to requested tolerance"
                   << " was NOT achieved...no output"<<endl;
       QDP_abort(1);}
    Eigenvalues.resize(nEigvecs);

     //  multiply each eigenvector by a
     //  phase so that the zero-th color component at site (0,0,0)
     //  is real and positive; this uniquely specifies the overall
     //  phase for each eigenvector, which is important since a
     //  change in these phases can change the effective Laph noise.
     //  Using this convention makes the results robust against
     //  eigenvector deletion and subsequent reconstruction.

    laphEigvecs.applyPhaseConvention();
    bulova.stop();

    QDPIO::cout << endl<<endl<<" SUCCESS: calculation time = "<<bulova.getTimeInSeconds() 
                << " secs" <<endl;
    if (chebyshevOrder>=2)
       QDPIO::cout << "  Chebyshev acceleration was used"<<endl;
    QDPIO::cout << "  Norm of matrix diagonalized was "
                   << Comp.getMatrixNormEstimate() << endl; 
    for (int i=0; i<nEigvecs; i++){  
       QDPIO::cout << "  Raw eigenvalue("<<i<<") = "
                   << Comp.getEigenvalue(i);             
       QDPIO::cout << "  Raw Residual("<<i<<") = "           
                   << Comp.getResidual(i) << endl;
       Eigenvalues[i]=Comp.getEigenvalue(i);}

         // perform a crucial check of the diagonalization
        // if Chebyshev acceleration was used, we need to get the
        // eigenvalues of -Delta

    double offmaxmag=0.0;
    laphEigvecs.resize(nEigvecs+1);
    laphEigvecs.setChebyshevAccelerationOff();
    for (int level=0;level<nEigvecs;level++){
       laphEigvecs.assignMatrixVectorMultiply(nEigvecs,level);
       double rr=sqrt(toDouble(
               laphEigvecs.InnerProductRealPart(nEigvecs,nEigvecs)));
       QDPIO::cout << "magnitude of diagonal element -Laplacian["
                   <<level<<","<<level<<"] = "<<rr<<endl;
       if (chebyshevOrder>=2) Eigenvalues[level]=rr;
       for (int row=0;row<nEigvecs;row++){
          if (row!=level){
             DComplex z=laphEigvecs.InnerProduct(row,nEigvecs);
             double rr=toDouble(sqrt(QDP::real(z)*QDP::real(z)+QDP::imag(z)*QDP::imag(z)));
             if (rr>offmaxmag) offmaxmag=rr;}}}
    QDPIO::cout << "Maximum magnitude of off-diagonal matrix elements"<<endl
                << "       of -Laplacian = " << offmaxmag <<endl<<endl;

#ifdef TESTING
    printLaphEigenvectors(laphEigvecs,t,*uPtr,nEigvecs,fout);
#endif

    iotimer.reset();iotimer.start(); 
    bulova.reset();bulova.start();
    DHput.open(TimeKey(t)); // open file, write header
    for (int k=0;k<nEigvecs;k++)
       DHput.putData(LevelKey(k),laphEigvecs.getVector(k));
    DHput.flush();  // finalize current file
    DHput.close(); 

    gHandler.clearData();  // clear smeared gauge time slices
    iotimer.stop(); bulova.stop();
    QDPIO::cout << "Time to write = "<<bulova.getTimeInSeconds()<<endl;
    iotime+=iotimer.getTimeInSeconds();
    QDPIO::cout << "Computation done for time = "<<t<<endl<<endl;

    }
  
#ifdef TESTING
 fout.close();
#endif

 rolex.stop();
 QDPIO::cout << endl<<"computeLaphEigenvectors: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "Total file I/O time = "<< iotime<<" sec"<<endl;
 QDPIO::cout << "ran successfully" << endl<<endl;

}

           // Compute covariantly displaced Laph Eigenvectors and output to file or NamedObjMap
           // Displacement specified by sequence of directions in "displace",
           // and each segment has length "disp_length"
           
void QuarkSmearingHandler::displaceLaphEigenvectors(
                         const set<DirPath>& displacements, int disp_length,
                         uint time_val, const string& smeared_gauge_file,
                         const std::string& disp_lapheigvec_file_stub)
{
 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);
 displaceLaphEigenvectors(displacements,disp_length,time_val,gHandler,
                          disp_lapheigvec_file_stub);
}


void QuarkSmearingHandler::displaceLaphEigenvectors(
                         const set<DirPath>& displacements, int disp_length,
                         uint time_val, GluonSmearingHandler& gHandler,
                         const std::string& disp_lapheigvec_file_stub)
{
 check_info_set("displaceLaphEigenvectors",1);

 int nTime = uPtr->getTimeExtent();
 uint nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 uint timeval = time_val % nTime;

 START_CODE();
 StopWatch rolex;
 rolex.start();
 double iotime=0.0;
 double dtime=0.0;

 QDPIO::cout << "Displacement of Laplacian eigenvectors commencing"
             << " in FieldSmearingHandler"<<endl;
 QDPIO::cout << "    number of eigenvectors = "<<nEigvecs<<endl;
 QDPIO::cout << "    time extent of lattice = "<<nTime<<endl;
 QDPIO::cout << "     time slice to compute = "<<timeval<<endl;

    //  make set that has nontrivial displacements
 set<DirPath> disp_paths(displacements);
 disp_paths.erase(DirPath());  // remove no-displacement
 StopWatch iotimer; iotimer.start(); 

      // get the undisplaced eigenvectors into memory
 QDPIO::cout << endl<<"Reading LapH eigenvectors into memory"<<endl;
 for (uint eigpair_num=0;eigpair_num<nEigvecs;++eigpair_num){ 
    const LattField& eigvec=getLaphEigenvector(timeval,eigpair_num);}
 iotimer.stop();
 iotime+=iotimer.getTimeInSeconds();

 if (disp_paths.empty()){
    QDPIO::cout << " No displaced eigenvectors required"<<endl;
    QDPIO::cout << "   Total file I/O time = "<< iotime<<" sec"<<endl;
    return;}

 smearedDispQuarkFileStub=tidyString(disp_lapheigvec_file_stub);
 if (smearedDispQuarkFileStub.empty()){
    QDPIO::cerr << "empty file name or stub in QuarkSmearingHandler"<<endl;
    QDP_abort(1);}

    // fire up the data handler

 FileListInfo writeFiles(smearedDispQuarkFileStub+"_time",timeval,timeval,false);
 DataPutHandlerMFO<QuarkSmearingHandler,TimeKey,LevelDispKey,
       LattField> *DHput=0;
 DHput=new DataPutHandlerMFO<QuarkSmearingHandler,TimeKey,LevelDispKey,
       LattField>(*this,writeFiles,"Laph--SmearedQuarkDispTimeFile",
                           "LaphEigenvectors",true,false);

     // do the computation for this time slice; put results into
     // files with the stub name given, concatenated with "_disp" and 
     // suffix  "_time.t" where  t = 0,1,...  

 DHput->open(TimeKey(timeval));  // open file, write header
 QDPIO::cout <<endl<<endl<<"Beginning computation for time = "
             <<timeval<<endl<<endl;

   // loop over the displacements to do  
 for (set<DirPath>::const_iterator dpt=disp_paths.begin();dpt!=disp_paths.end();++dpt){

    iotimer.reset(); iotimer.start();
      // get the stout-smeared gauge field link variables
    const multi1d<LatticeColorMatrix>& usmear
                        = gHandler.getSmearedGaugeFieldTimeSlice(timeval);
    iotimer.stop();
    iotime+=iotimer.getTimeInSeconds();

 //   {XmlBufferWriter xmlo;
 //   dpt->output(xmlo);
 //   QDPIO::cout << "Doing displacement "<<xmlo.str()<<endl;}

       // compute the displacement field for the full path
    LatticeColorMatrix covdisplacer;
    get_path_displacer(usmear,covdisplacer,*dpt,disp_length);
    LattField shifted;
    
       // loop over the eigenvectors

    for (uint eigpair_num=0;eigpair_num<nEigvecs;++eigpair_num){

       iotimer.reset();iotimer.start(); 
       const LattField& eigvec=getLaphEigenvector(timeval,eigpair_num);
       iotimer.stop();
       iotime+=iotimer.getTimeInSeconds();

       StopWatch dtimer; dtimer.start();
       do_path_shift(eigvec,shifted,*dpt,disp_length);
       LattField disp_eigvec=covdisplacer*shifted;
       dtimer.stop();
       dtime+=dtimer.getTimeInSeconds();

       iotimer.reset();iotimer.start();
       DHput->putData(LevelDispKey(eigpair_num,*dpt),disp_eigvec);
       iotimer.stop();
       iotime+=iotimer.getTimeInSeconds();
       } // loop over eigevectors
    }    // loop over displacements

 gHandler.clearData();  // clear smeared gauge time slices
 DHput->flush();
 DHput->close();
 QDPIO::cout << "Computation done for time = "<<timeval<<endl<<endl;

 delete DHput;
 rolex.stop();
 QDPIO::cout << endl<<"displaceLaphEigenvectors: total time = "<< rolex.getTimeInSeconds() << " secs" << endl;
 QDPIO::cout << "   Total file I/O time = "<< iotime<<" sec"<<endl;
 QDPIO::cout << "Total time to displace = "<< dtime<<" sec"<<endl;
 QDPIO::cout << "ran successfully" << endl<<endl;

}

    // Fires up ability to get displaced LapH eigenvectors
    // previously calculations.

void QuarkSmearingHandler::get_displaph_setup()
{
 if (smearedDispQuarkFileStub.empty()){
    QDPIO::cerr << "Empty smeared Displaced QuarkFile stub"<<endl;
    QDP_abort(1);}
 int nTime = uPtr->getTimeExtent();
 delete dph_ptr;
 FileListInfo dispFiles(smearedDispQuarkFileStub+"_time",0,nTime-1,false);
 try{ dph_ptr=new DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelDispKey,
                 LattField>(*this,dispFiles,"Laph--SmearedQuarkDispTimeFile",
                 "LaphEigenvectors");}
 catch(...){ 
   QDPIO::cerr << "unable to allocate data handler in QuarkSmearingHandler"<<endl;
   QDP_abort(1);}
}


  //  shifts a LattField by "dir_path"

void QuarkSmearingHandler::do_path_shift(const LattField& start,
                                         LattField& finish,
                                         const DirPath& dir_path, uint disp_length)
{
 uint nshifts=dir_path.Length()*disp_length;
 if (nshifts==0){
    QDPIO::cerr << "do_path_shift requires a non-trivial shift"<<endl;
    QDP_abort(1);}
 LattField temp;
 const LattField *v1=&start;
 LattField *v2=((nshifts%2)==0)? &temp: &finish;
 LattField *v3=((nshifts%2)==0)? &finish: &temp;
 DirPath::const_iterator seg=dir_path.begin();
 for (int k=0;k<dir_path.Length();++k){
    int dir=*seg;
    for (uint j=0;j<disp_length;++j){
       if (dir>0){
          *v2=shift(*v1,FORWARD,dir-1);}
       else{
          *v2=shift(*v1,BACKWARD,-dir-1);}
       if (v1!=&start) v3=(LattField*)(v1);
       v1=v2; v2=v3;}
    ++seg;}
}

   //  Computes the full covariant displacement field for "dir_path"

void QuarkSmearingHandler::get_path_displacer(const multi1d<LatticeColorMatrix>& usmear,
                                              LatticeColorMatrix& covdisplacer,
                                              const DirPath& dir_path, uint disp_length)
{
 uint nshifts=dir_path.Length()*disp_length;
 if (nshifts==0){
    QDPIO::cerr << "do_path_shift requires a non-trivial shift"<<endl;
    QDP_abort(1);}
 LatticeColorMatrix temp;
 DirPath::const_reverse_iterator seg=dir_path.rbegin();
 bool doprod=false;
 for (int k=0;k<dir_path.Length();++k){
    int dir=*seg;
    for (uint j=0;j<disp_length;++j){
       if (dir>0){
          if (doprod){
             temp=shift(covdisplacer,FORWARD,dir-1);
             covdisplacer=usmear[dir-1]*temp;}
          else{
             covdisplacer=usmear[dir-1];
             doprod=true;}}
       else{
          if (doprod){
             temp=adj(usmear[-dir-1])*covdisplacer;}
          else{
             temp=adj(usmear[-dir-1]);
             doprod=true;}
          covdisplacer=shift(temp,BACKWARD,-dir-1);}}
    ++seg;}
}



// *************************************************************


double QuarkSmearingHandler::estimateLargestLaplacianEigenvalue(
                                   const string& smeared_gauge_file)
{
 check_info_set("estimateLargestLaplacianEigenvalue",2);

 double lambda_max=0.0;
 int nTime = uPtr->getTimeExtent();
 int nEigvecs = 1;
 int mintime = uPtr->getMinTime();
 int maxtime = uPtr->getMaxTime();
 
    // fire up the GluonSmearingHandler

 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);

 QDPIO::cout << "Estimating largest eigenvalue of -Laplacian commencing"
             << " in FieldSmearingHandler"<<endl;
 QDPIO::cout << "            time extent of lattice = "<<nTime<<endl;
 QDPIO::cout << "            minTime = "<<mintime<<", maxTime = " << maxtime<<endl;

 START_CODE();
 StopWatch rolex;
 rolex.start();

     // do the computation for each time slice; put results into
     // files with the stub name given and suffix  ".time_t" where
     // t = 0,1,...  

 for (int t=mintime; t<=maxtime; t++){

    QDPIO::cout <<endl<<endl<<"Beginning computation for time = "<<t<<endl<<endl;
    const multi1d<LatticeColorMatrix>& usmear
                        = gHandler.getSmearedGaugeFieldTimeSlice(t);

         // use the Krylov-spectral restarted Lanczos method to compute
         // the eigenvalues/eigenvectors of the Laplacian on this time slice

    LaphVectorHandler laphEigvecs(usmear);
    laphEigvecs.setInitialToRandom();

    KSRLanczosEigenSolver Comp(laphEigvecs, nEigvecs, 1e-3, 24, 24, 'U');

    if (!Comp.wasSuccessful()){
       QDPIO::cout << endl<<endl<<"Convergence to requested tolerance"
                   << " was NOT achieved...no output"<<endl;}

    double this_max = Comp.getEigenvalue(0);
    QDPIO::cout << "largest eigenvalue on this time slice = "
                << this_max<<endl;

    lambda_max=max(lambda_max,this_max);
    
    }
  

 rolex.stop();
 QDPIO::cout << endl<<"computeLaphEigenvectors: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

 return lambda_max;
}

     // Provides access to the quark smearing eigenvectors
     // while in read mode.

const LattField& QuarkSmearingHandler::getLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 check_info_set("getLaphEigenvector",1);
 return dh_ptr->getData(TimeKey(timeslice),LevelKey(eigpair_num));
}

bool QuarkSmearingHandler::queryLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 check_info_set("queryLaphEigenvector",1);
 return dh_ptr->queryData(TimeKey(timeslice),LevelKey(eigpair_num));
}

void QuarkSmearingHandler::removeLaphEigenvector(
                                   int timeslice, int eigpair_num)
{
 if (!m_read_mode) return;
 return dh_ptr->removeData(TimeKey(timeslice),LevelKey(eigpair_num));
}

void QuarkSmearingHandler::clearLaphEigenvectors()
{
 if (!m_read_mode) return;
 dh_ptr->clearData();
 Eigenvalues.resize(0);
}



const LattField& QuarkSmearingHandler::getDisplacedLaphEigenvector(
                                   int timeslice, int eigpair_num, const DirPath& dirpath)
{
 check_info_set("getLaphEigenvector",1);
 if (dph_ptr==0) get_displaph_setup();
 return dph_ptr->getData(TimeKey(timeslice),LevelDispKey(eigpair_num,dirpath));
}

bool QuarkSmearingHandler::queryDisplacedLaphEigenvector(
                                   int timeslice, int eigpair_num, const DirPath& dirpath)
{
 check_info_set("queryLaphEigenvector",1);
 if (dph_ptr==0) get_displaph_setup();
 return dph_ptr->queryData(TimeKey(timeslice),LevelDispKey(eigpair_num,dirpath));
}

void QuarkSmearingHandler::removeDisplacedLaphEigenvector(
                                   int timeslice, int eigpair_num, const DirPath& dirpath)
{
 if (dph_ptr==0) return;
 return dph_ptr->removeDataMEM(TimeKey(timeslice),LevelDispKey(eigpair_num,dirpath));
}

void QuarkSmearingHandler::clearDisplacedLaphEigenvectors()
{
 if (dph_ptr) dph_ptr->clearDataMEM();
}


#elif (QDP_ND == 4)

           // The 3-d computation of the Laph eigenvectors produces
           // one file for each time slice, and each file contains all
           // eigen-levels up to the requested number of eigenvectors.
           // This routine reads all of these time files and combines
           // all time slices of EACH level into a 4D field in TheNamedObjMap.

void QuarkSmearingHandler::readTimeSlicesIntoNOM()
{
 check_info_set("readTimeSlices",1);
 double rtimer=0.0,wtimer=0.0;
 START_CODE();
 StopWatch rolex,seiko;
 rolex.start();
 seiko.start();
 
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 int mintime = uPtr->getMinTime();
 int maxtime = uPtr->getMaxTime();
 int nTime = maxtime - mintime + 1;
 int levelstart=0;

    // first, check to see if already in TheNamedObjMap
 FileListInfo checkFiles("LaphEig_level",levelstart,nEigvecs-1);
 bool flag=true;
 for (int level=levelstart;level<nEigvecs;++level){
    if (!(TheNamedObjMap::Instance().check(checkFiles.getFileName(level)))){
       flag=false; break;}}
 if (flag){ return;}

    // check that all "time" files are available
 for (int t=mintime;t<=maxtime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    if (!fileExists(oss.str())){
         QDPIO::cerr << "file "<<oss.str()<<" does not exist"<<endl;
         QDPIO::cerr << " readTimeSlices fails...."<<endl;
         QDP_abort(1);}}

   // create the data handlers (the "get" is only used for checking
   // the headers and extracting the eigenvalues)

 FileListInfo readFiles(smearedQuarkFileStub+"_time",mintime,maxtime);
 bool overwrite=true;
 FileListInfo writeFiles("NOM_LaphEig_level",levelstart,nEigvecs-1,overwrite);
 Eigenvalues.resize(nEigvecs); // storage for eigenvalues
 for (int iLev=0; iLev < nEigvecs; ++iLev)
   Eigenvalues[iLev].resize(nTime);

 m_read_mode=false;
 {DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelKey,
       LattField> DHget(*this,readFiles,"Laph--SmearedQuarkTimeFile",
                                 "LaphEigenvectors");}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Input check done; time="<<rtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();
 DataPutHandlerMFO<QuarkSmearingHandler,LevelKey,LevelKey,
       LattField> DHput(*this,writeFiles,"Laph--SmearedQuarkLevelFile",
                                 "LaphEigenvectors",true,false);
 seiko.stop();
 wtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Put handler set up; time="<<wtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();


    // use a "Handle" in case of exception
 multi1d< Handle<IOMap<LevelKey,
          TimeSliceOf<LattField> > > > rdm_ptr(nTime);
 for (int t=mintime;t<=maxtime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    rdm_ptr[t-mintime]=new IOMap<LevelKey,TimeSliceOf<LattField> >;
    try { rdm_ptr[t-mintime]->openReadOnly(oss.str(),"Laph--SmearedQuarkTimeFile");}
    catch(...) { QDPIO::cerr << "read failure in readTimeSlices"; 
                 QDP_abort(1);}}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Ready for read; time="<<rtimer<<" seconds"<<endl;

   //  now collect time slices of each level and output to file

 LattField laph_eigvecs = zero;
 TimeSliceOf<LattField> buffer(laph_eigvecs);

   // create each level file

 for (int level=levelstart;level<nEigvecs;level++){
    double rtime=0.0, wtime=0.0;
    seiko.reset();seiko.start();
    QDPIO::cout << "creating file for level "<<level<<endl;
    LevelKey key(level);
    DHput.open(key);  // creates files, writes header
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    for (int t=mintime;t<=maxtime;t++){
       buffer.setCurrentTime(t);
       try{ rdm_ptr[t-mintime]->get(key,buffer); }
       catch(...){ QDPIO::cerr << "Lookup error in readTimeSlices"<<endl;
                   QDP_abort(1);}}
    seiko.stop();
    rtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    DHput.putData(key,laph_eigvecs);
    DHput.flush();  // finalize current file
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    QDPIO::cout << " level "<<level<<"  total read time = "<< rtime  << " secs" << endl;
    QDPIO::cout << " level "<<level<<" total write time = "<< wtime  << " secs" << endl;
    rtimer+=rtime;
    wtimer+=wtime;
    }

 m_read_mode=true;
 smearedQuarkFileStub="NOM_LaphEig";
 rolex.stop();
 QDPIO::cout << endl<<"readTimeSlices: total time = "<< rolex.getTimeInSeconds()  << " secs" << endl;
 QDPIO::cout <<"  total read time = "<< rtimer  << " secs" << endl;
 QDPIO::cout <<" total write time = "<< wtimer  << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

}

           // The 3-d computation of the Laph eigenvectors produces
           // one file for each time slice, and each file contains all
           // eigen-levels up to the requested number of eigenvectors.
           // This routine re-organizes these files, producing one file
           // for each level, but each file contains all time slices.

void QuarkSmearingHandler::combineTimeSlices(int striping_factor, int striping_unit)
{
 check_info_set("combineTimeSlices",2);

 double rtimer=0.0,wtimer=0.0;
 START_CODE();
 StopWatch rolex,seiko;
 rolex.start();
 seiko.start();
 
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 int mintime = uPtr->getMinTime();
 int maxtime = uPtr->getMaxTime();
 int nTime = maxtime - mintime + 1;

    // check that all "time" files are available
 for (int t=mintime;t<=maxtime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    if (!fileExists(oss.str())){
         QDPIO::cerr << "file "<<oss.str()<<" does not exist"<<endl;
         QDPIO::cerr << " combineTimeSlices fails...."<<endl;
         QDP_abort(1);}}

    // determine the smallest level to start with (based on existing files)
 int levelstart;
 for (levelstart=nEigvecs-1;levelstart>=0;levelstart--){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << levelstart;
    if (fileExists(oss.str())) break;}
 levelstart++;
 if (levelstart==nEigvecs) return;

   // create the data handlers (the "get" is only used for checking
   // the headers and extracting the eigenvalues)

 FileListInfo readFiles(smearedQuarkFileStub+"_time",mintime,maxtime);
 FileListInfo writeFiles(smearedQuarkFileStub+"_level",levelstart,nEigvecs-1);
 Eigenvalues.resize(nEigvecs); // storage for eigenvalues
 for (int iLev=0; iLev < nEigvecs; ++iLev)
   Eigenvalues[iLev].resize(nTime);

 {DataGetHandlerMFO<QuarkSmearingHandler,TimeKey,LevelKey,
       LattField> DHget(*this,readFiles,"Laph--SmearedQuarkTimeFile",
                                 "LaphEigenvectors");}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Input check done; time="<<rtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();
 DataPutHandlerMFO<QuarkSmearingHandler,LevelKey,LevelKey,
       LattField> DHput(*this,writeFiles,"Laph--SmearedQuarkLevelFile",
                                 "LaphEigenvectors",true,false,
                                 striping_factor,striping_unit);
 seiko.stop();
 wtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Put handler set up; time="<<wtimer<<" seconds"<<endl;
 seiko.reset();seiko.start();


    // use a "Handle" in case of exception
 multi1d< Handle<IOMap<LevelKey,
          TimeSliceOf<LattField> > > > rdm_ptr(nTime);
 for (int t=mintime;t<=maxtime;t++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_time." << t;
    rdm_ptr[t-mintime]=new IOMap<LevelKey,TimeSliceOf<LattField> >;
    try { rdm_ptr[t-mintime]->openReadOnly(oss.str(),"Laph--SmearedQuarkTimeFile");}
    catch(...) { QDPIO::cerr << "read failure in combineTimeSlices"; 
                 QDP_abort(1);}}
 seiko.stop();
 rtimer+=seiko.getTimeInSeconds();
 QDPIO::cout << "Ready for read; time="<<rtimer<<" seconds"<<endl;

   //  now collect time slices of each level and output to file

 LattField laph_eigvecs = zero;
 TimeSliceOf<LattField> buffer(laph_eigvecs);

   // create each level file

 for (int level=levelstart;level<nEigvecs;level++){
    double rtime=0.0, wtime=0.0;
    seiko.reset();seiko.start();
    QDPIO::cout << "creating file for level "<<level<<endl;
    LevelKey key(level);
    DHput.open(key);  // creates files, writes header
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    for (int t=mintime;t<=maxtime;t++){
       buffer.setCurrentTime(t);
       try{ rdm_ptr[t-mintime]->get(key,buffer); }
       catch(...){ QDPIO::cerr << "Lookup error in combineTimeSlices"<<endl;
                   QDP_abort(1);}}
    seiko.stop();
    rtime+=seiko.getTimeInSeconds();
    seiko.reset();seiko.start();
    DHput.putData(key,laph_eigvecs);
    DHput.flush();  // finalize current file
    seiko.stop();
    wtime+=seiko.getTimeInSeconds();
    QDPIO::cout << " level "<<level<<"  total read time = "<< rtime  << " secs" << endl;
    QDPIO::cout << " level "<<level<<" total write time = "<< wtime  << " secs" << endl;
    rtimer+=rtime;
    wtimer+=wtime;
    }

 rolex.stop();
 QDPIO::cout << endl<<"combineTimeSlices: total time = "<< rolex.getTimeInSeconds()  << " secs" << endl;
 QDPIO::cout <<"  total read time = "<< rtimer  << " secs" << endl;
 QDPIO::cout <<" total write time = "<< wtimer  << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

}

*/
           // Compute the Laph Eigenvectors simultaneously on all
           // time slices (from mintime to maxtime in the *u_ptr).
           // Put results into TheNamedObjMap or to level files.

void QuarkSmearingHandler::computeLaphEigenvectors(
                         const LaphEigenSolverInfo& solver_info,
                         const string& smeared_gauge_file)
//                         int striping_factor, int striping_unit,
//                         int mintime, int maxtime)
{
 check_info_set("computeLaphEigenvectors",2);

 StopWatch rolex,bulova;
 rolex.start();
 double iotime=0.0;

 int nTime = uPtr->getTimeExtent();
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
// int firsttime=uPtr->getMinTime();
// int lasttime=uPtr->getMaxTime();
// if ((mintime>firsttime)&&(mintime<=lasttime)) firsttime=mintime;
// if ((maxtime<lasttime)&&(maxtime>=firsttime)) lasttime=maxtime;

    // check to see if output files exist already
 for (int v=0;v<nEigvecs;++v){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << v;
    if (fileExists(oss.str())){
       printLaph(make_str("file ",oss.str()," already exists:\n",
                 " ... will not overwrite so computation aborted\n\n"));
       return;}}

    // fire up the GluonSmearingHandler

 GluonSmearingHandler gHandler(*gSmearPtr,*uPtr,smeared_gauge_file);

 printLaph(make_str("Computation of Laplacian eigenvectors commencing",
                    " in FieldSmearingHandler"));
 printLaph(make_str("  number of requested eigenvectors = ",nEigvecs));
 printLaph(make_str("            time extent of lattice = ",nTime));
// QDPIO::cout << "     first time to compute = "<<firsttime<<endl;
// QDPIO::cout << "      last time to compute = "<<lasttime<<endl;

    // fire up the data handler

 FileListInfo writeFiles(smearedQuarkFileStub+"_level",0,nEigvecs-1);
 DataPutHandlerMFO<QuarkSmearingHandler,LevelKey,LevelKey,
       LattField> DHput(*this,writeFiles,"Laph--SmearedQuarkLevelFile",
                                 "LaphEigenvectors",true,false);
//                                 striping_factor,striping_unit);
      // read smeared gauge field 
 StopWatch iotimer; iotimer.start();
// const vector<LattField>& usmear = gHandler.getSmearedGaugeField();
 gHandler.copyDataToDevice();
 printLaph(make_str("\n\nSmeared gauge field obtained from ",smeared_gauge_file,"\n"));
 
 iotimer.stop();
 iotime+=iotimer.getTimeInSeconds();

         // use the Krylov-spectral restarted Lanczos method (blocked version) to compute
         // the eigenvalues/eigenvectors of minus the 3d Laplacian on all time slices
 bulova.start();
 QudaInvertParam eig_inv_param;
 QudaEigParam eig_param;
 solver_info.setQudaParam(eig_inv_param,eig_param,*qSmearPtr);

   // allocate space for eigenvectors and eigenvalues on the host
 vector<LattField> laphEigvecs(nEigvecs,FieldSiteType::ColorVector);
 vector<complex<double>> laphEigvals(nEigvecs*nTime);
 __complex__ double* h_evals=(__complex__ double*)(laphEigvals.data());
 vector<void*> h_evecs(nEigvecs);
 for (int v=0;v<nEigvecs;++v){
    h_evecs[v]=laphEigvecs[v].getDataPtr();}

 {complex<double> z(1.0/double(FieldNcolor*LayoutInfo::getLatticeNumSites()),0.0);
  setConstantField(laphEigvecs[0],z);}

   // now determine the eigenvectors!!
 eigensolveQuda(h_evecs.data(), h_evals, &eig_param);

     ///   CHECK CONVERGENCE!!!
/*
 if (!Comp.wasSuccessful()){
    QDPIO::cout << endl<<endl<<"Convergence to requested tolerance"
                << " was NOT achieved...no output"<<endl;
    QDP_abort(1);}
 Eigenvalues.resize(nEigvecs);
*/
 bulova.stop();
 printLaph(make_strf("\n\n SUCCESS: calculation time = %g seconds\n",
           bulova.getTimeInSeconds()));
 if (solver_info.getChebyshevOrder()>=2){
    printLaph(make_str("  Chebyshev acceleration was used"));}

     //  multiply each eigenvector by a
     //  phase so that the zero-th color component at site (0,0,0)
     //  is real and positive; this uniquely specifies the overall
     //  phase for each eigenvector, which is important since a
     //  change in these phases can change the effective Laph noise.
     //  Using this convention makes the results robust against
     //  eigenvector deletion and subsequent reconstruction.

 bulova.reset(); bulova.start();
 applyLaphPhaseConvention(laphEigvecs);
 bulova.stop();
 printLaph(make_strf(" Phase shift convention applied in %g seconds\n",
           bulova.getTimeInSeconds()));

 if (solver_info.getOutputVerbosity()>0){
    printLaph("Eigenvalues of -Laplacian:");
    for (int t=0;t<nTime;++t)
    for (int v=0;v<nEigvecs;++v){
       printLaph(make_strf(" eigenvalue[%3d](t=%3d) = %20.15f",v,t,
                 real(laphEigvals[nEigvecs*t+v])));}}    // NOT CORRECT YET

//  Do some checks of the results

 checkLaphEigvecComputation(laphEigvecs,gHandler.getSmearedGaugeField());

 iotimer.reset();iotimer.start(); 
// bulava.reset();bulava.start();  */
 for (int k=0;k<nEigvecs;k++){
    DHput.open(LevelKey(k)); // open file, write header
    DHput.putData(LevelKey(k),laphEigvecs[k]);
    DHput.close();     // finalize current file
    laphEigvecs[k].clear();} 

 gHandler.clearData();  // clear smeared gauge time slices
 iotimer.stop();
 printLaph(make_strf("Time to write laph eigenvectors = %g seconds",iotimer.getTimeInSeconds()));
 iotime+=iotimer.getTimeInSeconds();
 printLaph("Computation done \n");
  
 rolex.stop();
 printLaph(make_strf("computeLaphEigenvectors: total time = %g seconds",
           rolex.getTimeInSeconds()));
 printLaph(make_strf("Total file I/O time = %g seconds",iotime));
 printLaph("ran successfully\n");

}



// *************************************************************

    //  This is temporary; should be done on the device

void QuarkSmearingHandler::applyLaphPhaseConvention(vector<LattField>& laph_evecs)
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
 int nloctime=LayoutInfo::getRankLattSizes()[3];
 int tstride=LayoutInfo::getRankLattSizes()[0]*LayoutInfo::getRankLattSizes()[1]
            *LayoutInfo::getRankLattSizes()[2];
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
#endif
}

// **********************************************************************************

void QuarkSmearingHandler::checkLaphEigvecComputation(const vector<LattField>& laphEigvecs,
                                               const vector<LattField>& smeared_gauge_field)
{
 int nTime = uPtr->getTimeExtent();
 int nEigvecs = qSmearPtr->getNumberOfLaplacianEigenvectors();
 Array<double> block_eigvals(nEigvecs,nTime); 
 vector<double> offmaxmag(nTime,0.0);
 LattField minus_delta_on_ev(FieldSiteType::ColorVector);
 for (int v=0;v<nEigvecs;++v){
       // act with the spatial Laplacian on this eigenvector  M=-Delta  M x = lambda x
    applyMinusSpatialLaplacian(minus_delta_on_ev,laphEigvecs[v],smeared_gauge_field);
    vector<complex<double>> eigvals(getTimeSlicedInnerProducts(laphEigvecs[v],minus_delta_on_ev));
    vector<complex<double>> norm(getTimeSlicedInnerProducts(laphEigvecs[v],laphEigvecs[v]));
    for (int t=0;t<nTime;++t){
       printLaph(make_strf(" eigval[%d] (t=%d) is (%20.15f, %20.15f)",v,t,real(eigvals[t]),imag(eigvals[t])));
       printLaph(make_strf(" norm[%d] (t=%d) is (%20.15f, %20.15f)",v,t,real(norm[t]),imag(norm[t])));
       block_eigvals(v,t)=real(eigvals[t])/real(norm[t]);}
    for (int vv=0;vv<nEigvecs;++vv){
       if (vv!=v){
          vector<complex<double>> offdiag(getTimeSlicedInnerProducts(laphEigvecs[vv],minus_delta_on_ev));
          for (int t=0;t<nTime;++t){
             double offdiagmag=std::abs(offdiag[t]);
             if (offdiagmag>offmaxmag[t]){
                offmaxmag[t]=offdiagmag;}}}}}
 for (int t=0;t<nTime;++t){
    printLaph(make_str("\nCheck of Laph eigenvectors for time = ",t));
    printLaph(make_strf("  Maximum magnitude of offdiagonal <eigvec[v] | -Delta | eigvec[vv]> = %g",
              offmaxmag[t]));
    for (int v=0;v<nEigvecs;++v){
       printLaph(make_strf("   Eigenvalue[%3d] = %20.15f",v,block_eigvals(v,t)));}}
}
 //   LattField res(FieldSiteType::ColorVector);
   // lattice_addto(delta_on_ev, laphEigvecs[v],
   //                const complex<double>& zcoef=complex<double>(1.0,0.0))


/* int nblocks=lasttime-firsttime+1;
 const multi1d<double>& matnorms=Comp.getMatrixNormEstimate();
 for (int block=0;block<nblocks;block++){
    QDPIO::cout << "  Norm of matrix block "<<block<<" diagonalized was "
                   << matnorms[block] << endl; 
    for (int i=0; i<nEigvecs; i++){  
       QDPIO::cout << "  Raw eigenvalue("<<i<<") = "
                   << Comp.getEigenvalue(i)[block];             
       QDPIO::cout << "  Raw Residual("<<i<<") = "           
                   << Comp.getResidual(i)[block] << endl;}}
 for (int i=0; i<nEigvecs; i++){  
    Eigenvalues[i]=Comp.getEigenvalue(i);}
       // perform a crucial check of the diagonalization
       // if Chebyshev acceleration was used, we need to get the
       // eigenvalues of -Delta
 double offmaxmag=0.0;
 laphEigvecs.resize(nEigvecs+1);
 laphEigvecs.setChebyshevAccelerationOff();
 auto complGatherer = laphEigvecs.getComplexGathererToRoot();
 auto doublGatherer = laphEigvecs.getDoubleGathererToRoot();
 for (int level=0;level<nEigvecs;level++){
    laphEigvecs.assignMatrixVectorMultiply(nEigvecs,level);
    multi1d<double> rb=doublGatherer(laphEigvecs.InnerProductRealPart(nEigvecs,nEigvecs));
    if (Layout::primaryNode()) {
      for (int block=0;block<nblocks;block++){
	rb[block]=sqrt(rb[block]);
	QDPIO::cout << "magnitude of block "<<block<<" diagonal element -Laplacian["
	  <<level<<","<<level<<"] = "<<rb[block]<<endl;}
      if (chebyshevOrder>=2) Eigenvalues[level]=rb;
    }
    for (int row=0;row<nEigvecs;row++){
       if (row!=level){
          multi1d<DComplex> z=complGatherer(laphEigvecs.InnerProduct(row,nEigvecs));
	  if (Layout::primaryNode()) {
	    for (int b=0;b<nblocks;b++){
	      double rr=toDouble(sqrt(QDP::real(z[b])*QDP::real(z[b])+QDP::imag(z[b])*QDP::imag(z[b])));
	      if (rr>offmaxmag) offmaxmag=rr;}}}}}
 QDPIO::cout << "Maximum magnitude of off-diagonal matrix elements"<<endl
   << "       of -Laplacian = " << offmaxmag <<endl<<endl;

#ifdef TESTING
 TextFileWriter fout("lapheigvec4d.log");
 printLaphEigenvectors(laphEigvecs,*uPtr,nEigvecs,fout);
 fout.close();
#endif
*/


// *************************************************************

     // Provides access to the quark smearing eigenvectors
     // while in read mode.

const LattField& QuarkSmearingHandler::getLaphEigenvector(int eigpair_num)
{
 check_info_set("getLaphEigenvector",1);
 LevelKey key(eigpair_num);
 return dh_ptr->getData(key,key);
}

bool QuarkSmearingHandler::queryLaphEigenvector(int eigpair_num)
{
 check_info_set("queryLaphEigenvector",1);
 LevelKey key(eigpair_num);
 return dh_ptr->queryData(key,key);
}

void QuarkSmearingHandler::removeLaphEigenvector(int eigpair_num)
{
 if (!m_read_mode) return;
 LevelKey key(eigpair_num);
 return dh_ptr->removeData(key,key);
}

void QuarkSmearingHandler::clearLaphEigenvectors()
{
 if (!m_read_mode) return;
 dh_ptr->clearData();
 Eigenvalues.resize(0);
}

void QuarkSmearingHandler::closeLaphLevelFiles()
{
 if (!m_read_mode) return;
 dh_ptr->closeAll();
}

bool QuarkSmearingHandler::checkAllLevelFilesExist()
{
 check_info_set("checkAllLevelFilesExist");

 for (int level=0;level< qSmearPtr->getNumberOfLaplacianEigenvectors();level++){
    ostringstream oss;
    oss << smearedQuarkFileStub << "_level." << level;
    if (!fileExists(oss.str())){
       printLaph(make_strf("needed file %s does not exists",oss.str()));
       return false; }}
 return true;
}
/*
#endif


#if (QDP_ND == 3)

void QuarkSmearingHandler::writeHeader(XMLHandler& xmlw, const TimeKey& fkey, int suffix)
{
 if (suffix!=fkey.value){
    QDPIO::cerr << "time suffix does not match file key"<<endl;
    QDP_abort(1);}
 push(xmlw, "LaphEigenvectors");
 gSmearPtr->output(xmlw);
 uPtr->output(xmlw);
 qSmearPtr->output(xmlw);
 fkey.output(xmlw);
 write(xmlw,"Eigenvalues",Eigenvalues);
 pop(xmlw);
}

bool QuarkSmearingHandler::checkHeader(XMLHandler& xml_in, int suffix)
{
 if (xml_tag_count(xml_in,"LaphEigenvectors")!=1) return false;
 XMLHandler xmlr(xml_in,"./descendant-or-self::LaphEigenvectors");
 try {
    GaugeConfigurationInfo gauge_check(xmlr);
    uPtr->checkEqual(gauge_check);
    GluonSmearingInfo gsmear_check(xmlr);
    gSmearPtr->checkEqual(gsmear_check);
    QuarkSmearingInfo qsmear_check(xmlr);
    if (m_read_mode){
       qSmearPtr->checkOK(qsmear_check);}
    else{
       qSmearPtr->checkEqual(qsmear_check);}
    TimeKey time_check(xmlr);
    if (time_check.value!=suffix){throw(std::invalid_argument("error"));}
    if (!m_read_mode) read(xmlr,"Eigenvalues",Eigenvalues);}
 catch(...) { return false;}
 return true;
}

#elif (QDP_ND == 4)
*/
void QuarkSmearingHandler::writeHeader(XMLHandler& xmlw, const LevelKey& fkey, 
                                       int suffix)
{
 if (suffix!=fkey.value){
    errorLaph("level suffix does not match file key");}
 xmlw.set_root("LaphEigenvectors");
 XMLHandler xmltmp;
 gSmearPtr->output(xmltmp); xmlw.put_child(xmltmp);
 uPtr->output(xmltmp); xmlw.put_child(xmltmp);
 qSmearPtr->output(xmltmp); xmlw.put_child(xmltmp);
 fkey.output(xmltmp); xmlw.put_child(xmltmp);
// int nblocks=Eigenvalues[suffix].size();
// multi1d<double> evals(nblocks);
// for (int b=0;b<nblocks;b++) evals[b]=Eigenvalues[suffix][b];
// write(xmlw,"Eigenvalues",evals);
}

bool QuarkSmearingHandler::checkHeader(XMLHandler& xml_in, int suffix)
{
 if (xml_tag_count(xml_in,"LaphEigenvectors")!=1) return false;
 XMLHandler xmlr(xml_in,"LaphEigenvectors");
 try {
    //GaugeConfigurationInfo gauge_check(xmlr);   //  TO MATCH USQCD CHANGE LATER!!
    //uPtr->checkEqual(gauge_check);
    GluonSmearingInfo gsmear_check(xmlr);
    gSmearPtr->checkEqual(gsmear_check);
    QuarkSmearingInfo qsmear_check(xmlr);
    if (m_read_mode){
       qSmearPtr->checkEqual(qsmear_check);
       LevelKey level_check(xmlr);
       if (level_check.value!=suffix){throw(std::invalid_argument("error"));}}
  //  else{
  //     qSmearPtr->checkEqual(qsmear_check);
  //     TimeKey time_check(xmlr);
  //     if (time_check.value!=suffix){throw(std::invalid_argument("error"));}
  //     multi1d<double> eigvals;
  //     read(xmlr,"Eigenvalues",eigvals);
  //     for (int iL=0; iL<eigvals.size(); ++iL)
  //        Eigenvalues[iL][suffix-uPtr->getMinTime()]=eigvals[iL];}
    }
 catch(const std::exception& xp) { return false;}
 return true;
}
/*
#endif
*/

// ******************************************************************
}
