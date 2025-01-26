#include "tasks.h"
#include "perambulator_handler.h"
#include "stop_watch.h"

using namespace std;

namespace LaphEnv {


// ***********************************************************************
// *                                                                     *
// *  Task to compute the temporal correlations for a few low-lying      *
// *  hadrons using the full distillation method. NOTE: the stochastic   *
// *  LapH method is NOT used in this class. Only local hadron operators *
// *  are used.  The correlators are evaluated on a single configuration *
// *  for several time slices and **for several values of Nev**, the     *
// *  number of LapH eigenvectors. These correlator values can then be   *
// *  used later in effective masses to tune the LapH smearing.  In      *
// *  other words, to determine which value of Nev produces the least    *
// *  amount of excited-state contamination and the least amount of      *
// *  statistical error.  The pion and nucleon correlators are evaluated *
// *  when requested. Assumes smeared gauge field and the Laph           *
// *  eigenvectors are available, as well as the quark perambulators.    *
// *  XML input must have the form:                                      *
// *                                                                     *
// *    <Task>                                                           *
// *     <Name> LAPH_SMEAR_TUNE </Name>                                  *
// *     <GaugeConfigurationInfo> ... </GaugeConfigurationInfo>          *
// *     <GluonStoutSmearingInfo> ... </GluonStoutSmearingInfo>          *
// *     <QuarkActionInfo> ... </QuarkActionInfo>                        *
// *     <QuarkLaphSmearingInfo> ... </QuarkLaphSmearingInfo>            *
// *     <SmearedGaugeFileName> ... </SmearedGaugeFileName>              *
// *     <SmearedQuarkFileStub> ... </SmearedQuarkFileStub>              *
// *     <QuarkPerambulators>                                            *
// *         <FileListInfo>  ... </FileListInfo>   (for quarks)          *
// *         <UpperSpinComponentsOnly/> (optional for nucleons)          *
// *     </QuarkPerambulators>                                           *
// *     <ZMomentumQuantum> 3 </ZMomentumQuantum>                        *
// *     <TimeSeparations> 3 4 5 </TimeSeparations>                      *
// *     <SmearingNeigvalues>16 20 24 ...</SmearingNeigvalues>           *
// *     <OutputData> out.dat </OutputData>   (name of output file)      *
// *                                                                     *
// *     <DoPionA/>   (Evaluate all Nev all at the same time)            *
// *     <DoPionS/>   (Evaluate Nev one at a time in sequence)           *
// *     <DoNucleonP/> (Triplets on all nodes, correlator on primary)    *
// *     <DoNucleonA/> (Triplets and correlator on all nodes)            *
// *     <DoNucleonS/> (Triplets and correlator on all nodes but         *
// *                     in sequence for Nev so less memory used)        *
// *                                                                     *
// *     <DoPion/> is same as <DoPionA/>                                 *
// *     <DoNucleon/> is same as <DoNucleonA/>                           *
// *                                                                     *
// *    </Task>                                                          *
// *                                                                     *
// *  To include the pion correlator, include the tag <DoPion/>.  To     *
// *  include the nucleon correlator, include either the tag             *
// *  <DoNucleonA/> or <DoNucleonP/> or <DoNucleonS/>.  If you have a    *
// *  small lattice and only a few Nev to compute, choose <DoNucleonP/>  *
// *  for which the baryon triplets are computed using all MPI processes,*
// *  but the correlators are evaluated only on the primary process.  If *
// *  you have a large lattice and a large number of Nev to compute,     *
// *  choose <DoNucleonA/> which computes the baryon triplets and the    *
// *  correlators distributing the work load and memory usage over all   *
// *  MPI processes. <DoNucleonA/> and <DoNucleonS/> distribute across   *
// *  all nodes, but differ in how the results for different Nev are     *
// *  done: <DoNucleonA/> computes all Nev results simultaneously and    *
// *  uses more memory, whereas <DoNucleonS/> computes the Nev in        *
// *  sequence, using less memory.                                       *
// *                                                                     *
// *  Only simple single-site operators are used:                        *
// *    Pion:    d-bar u   (1 3 + 2 4 - 3 1 - 4 2)   chi - psi format    *
// *    Nucleon:   u u d   ( 1 2 1 - 1 1 2 )           G1g_1 SS 0        *
// *    Nucleon:   u u d   ( 2 2 1 - 1 2 2 )           G1g_2 SS 0        *
// *                                                                     *
// *  The pion and nucleon are taken to have momentum 2*Pi*Zmom/Lz       *
// *  where Zmom is specified in the <ZMomentumQuantum> tag.             *
// *                                                                     *
// ***********************************************************************


void doLaphSmearTune(XMLHandler& xmltask)
{ /*
 XMLHandler xmlr(xmltask);
 string gauge_xml;
 GaugeConfigurationInfo gaugeinfo(xmlr,gauge_xml);
// gaugeinfo.setRelaxedChecking();
 GluonSmearingInfo gSmear(xmlr);
 QuarkActionInfo quark(xmlr,gaugeinfo);
 QuarkSmearingInfo qSmear(xmlr);

 string smeared_gauge_file;
 xmlread(xmlr,"SmearedGaugeFileName",smeared_gauge_file,"SmearTune");
 smeared_gauge_file=tidyString(smeared_gauge_file);

 string smeared_quark_filestub;
 xmlread(xmlr,"SmearedQuarkFileStub",smeared_quark_filestub,"SmearTune");
 smeared_quark_filestub=tidyString(smeared_quark_filestub);
 FileListInfo qfiles(xmlr);
 
 vector<uint> neig_smear;
 xmlread(xmlr,"SmearingNeigvalues",neig_smear,"SmearTune");
    // convert "neig_smear" to a std::vector and sort
// vector<uint> smear_neigs(neig_smear.size());
// for (uint k=0;k<neig_smear.size();++k){
//    smear_neigs[k]=neig_smear[k];}
 std::sort(smear_neigs.begin(), smear_neigs.end());
 std::unique(smear_neigs.begin(), smear_neigs.end());
 if (smear_neigs.size()==0){
    printLaph("Nothing to compute in SmearTune");
    return;}
 if (smear_neigs[0]==0){
    printLaph("All Neigs>0 required in SmearTune");
    return;}

 string outdata;
 xmlread(xmlr,"OutputData",outdata,"SmearTune");
 outdata=tidyString(outdata);
 if (outdata.length()==0){
    QDPIO::cerr << "Empty string for output data file; aborting"<<endl;
    QDP_abort(1);}

 int zmom;
 xmlread(xmlr,"ZMomentumQuantum",zmom,"SmearTune");
 int Textent=gaugeinfo.getTimeExtent();
 multi1d<int> tseps;
 xmlread(xmlr,"TimeSeparations",tseps,"SmearTune");
 set<int> tsepvals;
 for (uint k=0;k<tseps.size();++k){
    if ((tseps[k]>=0)&&(tseps[k]<Textent)){
       tsepvals.insert(tseps[k]);}}
 if (tsepvals.size()==0){
    QDPIO::cout << "No valid time separations"<<endl;
    return;}

 bool upper_spin_only=false;
 if (xml_tag_count(xmlr,"UpperSpinComponentsOnly")>=1){
    upper_spin_only=true;}

 QDPIO::cout << endl << endl;
 QDPIO::cout << " ***********************************************************"<<endl;
 QDPIO::cout << " *                                                         *"<<endl;
 QDPIO::cout << " *      Task: Laph Smear Tune                              *"<<endl;
 QDPIO::cout << " *                                                         *"<<endl;
 QDPIO::cout << " ***********************************************************"<<endl;
 QDPIO::cout << endl;
 QDPIO::cout <<endl<<gaugeinfo.output()<<endl;
 QDPIO::cout << "XML header in the gauge configuration:"<<endl;
 QDPIO::cout << gauge_xml<<endl<<endl;
 QDPIO::cout <<endl<<endl<<"Gluon Smearing:"<<endl<<gSmear.output()<<endl<<endl;
 QDPIO::cout <<endl<<endl<<"Quark Smearing:"<<endl<<qSmear.output()<<endl<<endl;
 QDPIO::cout <<"SmearedQuarkFileStub: "<<smeared_quark_filestub<<endl;
 QDPIO::cout << endl<<"QuarkAction:"<<endl<< quark.output()<<endl<<endl;
 QDPIO::cout << "ZMomentumQuantum: "<<zmom<<endl;
 QDPIO::cout << "Time Separations: ";
 for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt) QDPIO::cout << " "<<*tt;
 QDPIO::cout <<endl<< "NeigenvectorSmear: ";
 for (int k=0;k<smear_neigs.size();++k) QDPIO::cout << " "<<smear_neigs[k];
 QDPIO::cout <<endl<<endl<< "Output file: "<<outdata<<endl<<endl;

 XmlBufferWriter xml_out;
 push(xml_out,"LAPH_SMEAR_TUNE");
 gaugeinfo.output(xml_out);
 gSmear.output(xml_out);
 qSmear.output(xml_out);
 quark.output(xml_out);
 pop(xml_out);
 xmlout << xml_out.str();

 bool dopionA=false, dopionS=false, donucleonP=false, donucleonA=false, donucleonS=false;
 if (xml_tag_count(xmlr,"DoPion")>0) dopionA=true;
 if (xml_tag_count(xmlr,"DoPionA")>0) dopionA=true;
 if (xml_tag_count(xmlr,"DoPionS")>0) dopionS=true;
 if (xml_tag_count(xmlr,"DoNucleon")>0) donucleonA=true;
 if (xml_tag_count(xmlr,"DoNucleonA")>0) donucleonA=true;
 if (xml_tag_count(xmlr,"DoNucleonS")>0) donucleonS=true;
 if (xml_tag_count(xmlr,"DoNucleonP")>0) donucleonP=true;

 if (upper_spin_only){
    if ((dopionA)||(dopionS)){
        QDPIO::cout << "To compute pions, must have all four Dirac spin components"<<endl;
        return;}}

 PerambulatorHandler Q(gaugeinfo,gSmear,qSmear,quark,qfiles,
                       smeared_quark_filestub,upper_spin_only);
 set<int> source_times;
 Q.getSourceTimes(source_times);

 if (dopionA){
    calcPionCorrelatorA(Q,zmom,smear_neigs,source_times,outdata+"_pion",tsepvals);}
 else if (dopionS){
    calcPionCorrelatorS(Q,zmom,smear_neigs,source_times,outdata+"_pion2",tsepvals);}
 QMP_barrier();
 if (donucleonA){
    calcNucleonCorrelatorBlocked(Q,zmom,smear_neigs,source_times,outdata+"_nucleon",tsepvals);}
 else if (donucleonS){
    calcNucleonCorrelatorBlocked2(Q,zmom,smear_neigs,source_times,outdata+"_nucleon",tsepvals);}
 else if (donucleonP){
    calcNucleonCorrelator(Q,zmom,smear_neigs,source_times,outdata+"_nucleon",tsepvals);} */
}


// ********************************************************************
   
/*
void LaphSmearTuneInlineMeas::getLaphEigenvectors(PerambulatorHandler& Q, int timeval,
                                                  std::vector<const LatticeColorVector*>& laph_eigvecs)
{
 StopWatch bulova;
 bulova.start();
 clearLaphEigenvectors(Q,laph_eigvecs);
// uint nEigvecs = Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 uint nEigvecs = m_nLaphEvs;
 laph_eigvecs.resize(nEigvecs);
 QuarkSmearingHandler& qsHandler=Q.getQuarkSmearingHandler();
      // get the eigenvectors into memory
 QDPIO::cout << endl<<"Reading LapH eigenvectors into memory for time "<<timeval<<endl;
 for (uint eigpair_num=0;eigpair_num<nEigvecs;++eigpair_num){ 
    const LatticeColorVector& eigvec=qsHandler.getLaphEigenvector(timeval,eigpair_num);
    laph_eigvecs[eigpair_num]=&eigvec;}
 bulova.stop();
 QDPIO::cout << "LapH eigenvectors read in time "<<bulova.getTimeInSeconds()<<" seconds"<<endl;
}


void LaphSmearTuneInlineMeas::clearLaphEigenvectors(PerambulatorHandler& Q, 
                                                    std::vector<const LatticeColorVector*>& laph_eigvecs)
{
 laph_eigvecs.clear();
 Q.getQuarkSmearingHandler().clearLaphEigenvectors();
}


// **************************************************************************************

   // Meson doublets are small enough to easily fit into the memory of one node,
   // so we do not spread meson_doublets over several nodes.  meson_doublets will reside
   // only on the primary node.
   //
   //          M[m,k](t) = sum_x v[a,m](x,t)^* v[a,k](x,t)
   
   
void LaphSmearTuneInlineMeas::evaluateMesonDoublets(PerambulatorHandler& Q, int pz,
                                        std::vector<const LatticeColorVector*>& laph_eigvecs,
                                        multi2d<Complex>& meson_doublets, uint indmin, uint indlim)
{
 QMP_barrier();
 StopWatch bulova;
 bulova.start();
// uint nEv=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 LatticeComplex* momphase=0;
 LatticeComplex* singlet=new LatticeComplex;
 LatticeComplex* singlet0=singlet;

   // set up the exp(-I*pz*z) phases for momentum in z-direction
 if (pz!=0){
    momphase=new LatticeComplex;
    singlet0=new LatticeComplex;
    double twopi=6.2831853071795864770;
    int zdir=2;
    Double momquantum_z=twopi*pz / toDouble(Layout::lattSize()[zdir]);
    LatticeReal p_dot_r = Layout::latticeCoordinate(zdir)*momquantum_z; 
    *momphase = cmplx(cos(p_dot_r),-sin(p_dot_r));}

 uint nEv=m_nLaphEvs;
 int nstore=int(indlim)-int(indmin);
 if (nstore>0){
    meson_doublets.resize(nstore,nEv);} 
 const int nsites = Layout::sitesOnNode();
 for (uint aEv=0; aEv<nEv; ++aEv) {
   bool save=((nstore>0)&&(aEv>=indmin)&&(aEv<indlim));
   const LatticeColorVector& aLHS = *(laph_eigvecs[aEv]);
   for (uint bEv=0; bEv<nEv; ++bEv) {
     const LatticeColorVector& aRHS = *(laph_eigvecs[bEv]);
     *singlet0 = localInnerProduct(aLHS, aRHS);
     if (pz!=0){
        *singlet=(*momphase)*(*singlet0);}
     Complex sum=zero;
     for (int i=0;i<nsites;i++){
        sum.elem()+=singlet->elem(i);}
     QDPInternal::globalSum(sum);
     if (save){
        meson_doublets(aEv-indmin,bEv)=sum;}}}
 bulova.stop();
 delete momphase;
 delete singlet;
 if (pz!=0){
    delete singlet0;}
 QDPIO::cout << "Meson doublets computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
 QMP_barrier();
}


   // Baryon triplets are large enough that they might need to spread over several nodes.
   // Instead of using a multi3d for baryon_triplets, we use a multi1d<multi2d> in case
   // we need to split it among different nodes.  
   
   //       B[k,l,m](t) = sum_x epsilson(a,b,c)  v[a,k](x,t) v[b,l](x,t) v[c,m](x,t)
   
   // Triplets are antisymmetric in k,l,m.
   
         //  This version is mainly for debugging purposes.


void LaphSmearTuneInlineMeas::evaluateBaryonTriplets(PerambulatorHandler& Q, int pz,
                                                     std::vector<const LatticeColorVector*>& laph_eigvecs,
                                                     multi1d<multi2d<Complex>>& baryon_triplets)
{
 QMP_barrier();
 StopWatch bulova;
 bulova.start();
// uint nEv=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 LatticeComplex* momphase=0;
 LatticeComplex* singlet=new LatticeComplex;
 LatticeComplex* singlet0=singlet;

   // set up the exp(-I*pz*z) phases for momentum in z-direction
 if (pz!=0){
    momphase=new LatticeComplex;
    singlet0=new LatticeComplex;
    double twopi=6.2831853071795864770;
    int zdir=2;
    Double momquantum_z=twopi*pz / toDouble(Layout::lattSize()[zdir]);
    LatticeReal p_dot_r = Layout::latticeCoordinate(zdir)*momquantum_z; 
    *momphase = cmplx(cos(p_dot_r),-sin(p_dot_r));}

 uint nEv=m_nLaphEvs;
 if (Layout::primaryNode()){
     baryon_triplets.resize(nEv);
     for (int k=0;k<nEv;++k){
        baryon_triplets[k].resize(nEv,nEv);}}
 const int nsites = Layout::sitesOnNode();
 for (uint aEv=0; aEv<nEv; ++aEv) {
   const LatticeColorVector& v1 = *(laph_eigvecs[aEv]);
   for (uint bEv=aEv+1; bEv<nEv; ++bEv) {
     const LatticeColorVector& v2 = *(laph_eigvecs[bEv]);
     LatticeColorVector vdiq = colorCrossProduct(v1, v2);
     for (uint cEv=bEv+1; cEv<nEv; ++cEv) {
        const LatticeColorVector& v3 = *(laph_eigvecs[cEv]);
        *singlet0 = colorVectorContract(vdiq, v3);
        if (pz!=0){
           *singlet=(*momphase)*(*singlet0);}
        Complex sum=zero;
        for (int i=0;i<nsites;i++){
           sum.elem()+=singlet->elem(i);}
        QDPInternal::globalSum(sum);
        if (Layout::primaryNode()){
           baryon_triplets[aEv](bEv,cEv)=sum;}}}}
   // now fill in other elements using antisymmetry of the triplets
 if (Layout::primaryNode()){
    for (uint aEv=0; aEv<nEv; ++aEv)
    for (uint bEv=0; bEv<nEv; ++bEv){
       baryon_triplets[aEv](aEv,bEv)=zero;
       baryon_triplets[aEv](bEv,aEv)=zero;
       baryon_triplets[aEv](bEv,bEv)=zero;}
    for (uint aEv=0; aEv<nEv; ++aEv) {
      for (uint bEv=0;bEv<aEv; ++bEv) {
        for (uint cEv=0; cEv<bEv; ++cEv) {
           baryon_triplets[aEv](bEv,cEv)=-baryon_triplets[cEv](bEv,aEv);
           baryon_triplets[aEv](cEv,bEv)=baryon_triplets[cEv](bEv,aEv);
           baryon_triplets[bEv](aEv,cEv)=baryon_triplets[cEv](bEv,aEv);
           baryon_triplets[cEv](aEv,bEv)=-baryon_triplets[cEv](bEv,aEv);
           baryon_triplets[bEv](cEv,aEv)=-baryon_triplets[cEv](bEv,aEv);}}}}
 bulova.stop();
 delete momphase;
 delete singlet;
 if (pz!=0){
    delete singlet0;}
 QDPIO::cout << "Baryon triplets computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
 QMP_barrier();
}


   //  This version evaluates the baryon triplets, but the results are spread throughout
   //  the nodes based on the first index.

void LaphSmearTuneInlineMeas::evaluateBaryonTripletsBlocked(PerambulatorHandler& Q, int pz, uint ind1min, uint ind1lim,
                                                     std::vector<const LatticeColorVector*>& laph_eigvecs,
                                                     multi1d<multi2d<Complex>>& baryon_triplets)
{
 QMP_barrier();
 StopWatch bulova;
 bulova.start();
// uint nEv_orig=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 LatticeComplex* momphase=0;
 LatticeComplex* singlet=new LatticeComplex;
 LatticeComplex* singlet0=singlet;

   // set up the exp(-I*pz*z) phases for momentum in z-direction
 if (pz!=0){
    momphase=new LatticeComplex;
    singlet0=new LatticeComplex;
    double twopi=6.2831853071795864770;
    int zdir=2;
    Double momquantum_z=twopi*pz / toDouble(Layout::lattSize()[zdir]);
    LatticeReal p_dot_r = Layout::latticeCoordinate(zdir)*momquantum_z; 
    *momphase = cmplx(cos(p_dot_r),-sin(p_dot_r));}

 uint nEv=m_nLaphEvs;
 uint nMPI=Layout::numNodes();
 uint ind1size=ind1lim-ind1min;
 baryon_triplets.resize(ind1size);
 for (int k=0;k<ind1size;++k){
     baryon_triplets[k].resize(nEv,nEv);}
 const int nsites = Layout::sitesOnNode();
 for (uint aEv=0; aEv<nEv; ++aEv) {
   const LatticeColorVector& v1 = *(laph_eigvecs[aEv]);
   for (uint bEv=aEv+1; bEv<nEv; ++bEv) {
     const LatticeColorVector& v2 = *(laph_eigvecs[bEv]);
     LatticeColorVector vdiq = colorCrossProduct(v1, v2);
     for (uint cEv=bEv+1; cEv<nEv; ++cEv) {
        const LatticeColorVector& v3 = *(laph_eigvecs[cEv]);
        *singlet0 = colorVectorContract(vdiq, v3);
        if (pz!=0){
           *singlet=(*momphase)*(*singlet0);}
        Complex sum=zero;
        for (int i=0;i<nsites;i++){
           sum.elem()+=singlet->elem(i);}
        QDPInternal::globalSum(sum);
        if ((aEv>=ind1min)&&(aEv<ind1lim)){
           baryon_triplets[aEv-ind1min](bEv,cEv)=sum;
           baryon_triplets[aEv-ind1min](cEv,bEv)=-sum;}
        if ((bEv>=ind1min)&&(bEv<ind1lim)){
           baryon_triplets[bEv-ind1min](cEv,aEv)=sum;
           baryon_triplets[bEv-ind1min](aEv,cEv)=-sum;}
        if ((cEv>=ind1min)&&(cEv<ind1lim)){
           baryon_triplets[cEv-ind1min](aEv,bEv)=sum;
           baryon_triplets[cEv-ind1min](bEv,aEv)=-sum;}}}}
 QMP_barrier();
   // now fill in other elements using antisymmetry of the triplets
 for (uint aEv=ind1min; aEv<ind1lim; ++aEv)
 for (uint bEv=0; bEv<nEv; ++bEv){
    baryon_triplets[aEv-ind1min](aEv,bEv)=zero;
    baryon_triplets[aEv-ind1min](bEv,aEv)=zero;
    baryon_triplets[aEv-ind1min](bEv,bEv)=zero;}
 bulova.stop();
 delete momphase;
 delete singlet;
 if (pz!=0){
    delete singlet0;}
 QDPIO::cout << "Baryon blocked triplets computed: time = "
             <<bulova.getTimeInSeconds()<<" seconds"<<endl; 
 QMP_barrier();
}


   //  This routine cycles the baryon_triplet blocks.  The data on node "n" goes to node "n+1",
   //  while the data on the last rank goes to the first rank.  Assume m_ntripactive is even.

void LaphSmearTuneInlineMeas::cycleBaryonTripletsBlocked(multi1d<multi2d<Complex>>& baryon_triplets, 
                                                         int cyclenum)
{
 QMP_barrier();
 int thisnode=Layout::nodeNumber();
 int outnode=0,innode=0,datanode=0;
 uint inind1min=0,inind1lim=0;
 uint blocksize=0;
 uint nEv=0;
 multi1d<multi2d<Complex>> tmp;

 if (thisnode<m_ntripactive){
    nEv=baryon_triplets[0].size1();
    outnode=thisnode+1;
    if (outnode==m_ntripactive) outnode=0;
    innode=thisnode-1;
    if (innode==-1) innode=m_ntripactive-1;
    datanode=(thisnode+m_ntripactive-cyclenum)%m_ntripactive;
    getIndex1Limits(datanode,nEv,inind1min,inind1lim);
    blocksize=inind1lim-inind1min;
    tmp.resize(blocksize);
    for (uint b=0;b<blocksize;b++){
       tmp[b].resize(nEv,nEv);}

    //  first set of transfers:  even nodes send, odd receive
    if ((thisnode%2)==0){
       sendBaryonTriplets(baryon_triplets,outnode);}
    else{
       recvBaryonTriplets(tmp,innode);}}
 QMP_barrier();

    //  second set of transfers:  odd nodes send, even receive
 if (thisnode<m_ntripactive){
    if ((thisnode%2)==1){
       sendBaryonTriplets(baryon_triplets,outnode);}
    else{
       recvBaryonTriplets(tmp,innode);}}
 QMP_barrier();
   // all transfer done; now copy tmp to baryon_triplets
 if (thisnode<m_ntripactive){
    baryon_triplets=tmp;}
}


   //  Sends the baryon triplets from this process to another process.  Waits until the send is completed.
   //  Sends in 16MB chunks.

void LaphSmearTuneInlineMeas::sendBaryonTriplets(const multi1d<multi2d<Complex>>& baryon_triplets, 
                                                 int dest_process)
{
 uint blocksize=baryon_triplets.size();
 uint nEv=baryon_triplets[0].size1();
 uint Nsend=2*blocksize*nEv*nEv;
 double *sendbuf=new double[Nsend];
 uint count=0;
 for (uint b=0;b<blocksize;++b)
 for (uint i=0;i<nEv;++i)
 for (uint j=0;j<nEv;++j){
    sendbuf[count++]=baryon_triplets[b](i,j).elem().elem().elem().real();
    sendbuf[count++]=baryon_triplets[b](i,j).elem().elem().elem().imag();}
 int chunksize=16777216; // 16MB
 int remaindersize=Nsend*sizeof(double); 
 char *ptr=(char*)sendbuf;
 if (remaindersize<chunksize) chunksize=remaindersize;
 while (ptr!=0){
    QDPInternal::sendToWait(ptr,dest_process,chunksize);
    ptr+=chunksize; remaindersize-=chunksize;
    if (chunksize>remaindersize) chunksize=remaindersize;
    if (chunksize==0) ptr=0;}
 delete sendbuf;
}

   //  Receives the baryon triplets at this process from another process.  Waits until the receive is completed.
   //  baryon_triplet should already be appropriately sized.

void LaphSmearTuneInlineMeas::recvBaryonTriplets(multi1d<multi2d<Complex>>& baryon_triplets, 
                                                 int src_process)
{
 uint blocksize=baryon_triplets.size();
 uint nEv=baryon_triplets[0].size1();
 uint Nsend=2*blocksize*nEv*nEv;
 double *recvbuf=new double[Nsend];
 int chunksize=16777216; // 16MB
 int remaindersize=Nsend*sizeof(double); 
 char *ptr=(char*)recvbuf;
 if (remaindersize<chunksize) chunksize=remaindersize;
 while (ptr!=0){
    QDPInternal::recvFromWait(ptr,src_process,chunksize);
    ptr+=chunksize; remaindersize-=chunksize;
    if (chunksize>remaindersize) chunksize=remaindersize;
    if (chunksize==0) ptr=0;}
 uint count=0;
 for (uint b=0;b<blocksize;++b)
 for (uint i=0;i<nEv;++i)
 for (uint j=0;j<nEv;++j){
    Double re=recvbuf[count++];
    Double im=recvbuf[count++];
    baryon_triplets[b](i,j)=cmplx(re,im);}
 delete recvbuf;
}



    // Compute the min and lim for index 1 for this MPI process.  Return the active number
    // of processes needed to handle the baryon triplets.
 
uint LaphSmearTuneInlineMeas::getIndex1Limits(uint thisnode, uint nEv, uint& ind1min, uint& ind1lim)
{
 int nMPI=Layout::numNodes();
 uint b1=nEv/nMPI;
 uint b2=nEv%nMPI;
 ind1min=(thisnode<=b2)?(thisnode*(b1+1)):(b2*(b1+1)+(thisnode-b2)*b1);
 ind1lim=(thisnode<b2)?(ind1min+b1+1):(ind1min+b1);
 return (b1==0)?b2:nMPI;
}


  // Computing the mesons involves a summation over n1=0..nEv-1 and n2=0..nEv-1,
  // among other summations.  This routine takes "nMPI", "nEv", and determines
  // the ranges of n1 and n2 to be handled by the current node "thisnode".
  // This node will handle n1=n1min..n1lim-1 and n2=n2min..n2lim-1

void LaphSmearTuneInlineMeas::getMesonEvIndexLimits(uint nEv, uint thisnode, 
                                       uint& n1min, uint& n1lim, uint& n2min, uint& n2lim)
{
 int nMPI=Layout::numNodes();
 uint nrows=int(floor(0.01+sqrt(double(nMPI))));
 uint ncols=nMPI/nrows;
 uint diff=nMPI%(nrows*ncols);
 while ((nrows>=1)&&(diff>0)){
    --nrows; ncols=nMPI/nrows;
    diff=nMPI%(nrows*ncols);}
 uint nodecol=thisnode%ncols;
 uint noderow=thisnode/ncols;
 uint b1=nEv/nrows;
 uint b2=nEv%nrows;
 n1min=(noderow<=b2)?(noderow*(b1+1)):(b2*(b1+1)+(noderow-b2)*b1);
 n1lim=(noderow<b2)?(n1min+b1+1):(n1min+b1);
 b1=nEv/ncols;
 b2=nEv%ncols;
 n2min=(nodecol<=b2)?(nodecol*(b1+1)):(b2*(b1+1)+(nodecol-b2)*b1);
 n2lim=(nodecol<b2)?(n2min+b1+1):(n2min+b1);
}



// ********************************************************************************************

   //  NOTE: below is not very optimized code since these routines are not meant
   //  to be run many many times.  Should be run just a few times when determining
   //  the smearing.
   
   //  For the pion correlator, the pion operator is taken to have the form
   //
   //       d-bar u  gamma_5 pion   (1 3 + 2 4 - 3 1 - 4 2)   chi - psi format
   //
   //  which gives us
   //
   //    C(t-t0) =   M[m,k](t)  PQ(k,alpha,t | k',beta,t0)  M[k',m'](t0)  PQ(m,alpha,t | m' beta,t0)^*
   //
   //  "neig_smear" contains the Nev values to get results for, and should be sorted in
   //  ascending order.

      // This version computes the different Nev values in sequence one at a time.
      // Should be faster for just using only a few values of Nev.

void LaphSmearTuneInlineMeas::calcPionCorrelatorS(PerambulatorHandler& Q, int pz,
                                 const vector<uint>& neig_smear, const set<int>& source_times,
                                 const string& outdata, const set<int>& tsepvals)
{
 QMP_barrier();
 StopWatch bulova;
 bulova.start();
 int Textent=Q.getGaugeConfigurationInfo().getTimeExtent();
 uint nEv_orig=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 m_nLaphEvs=neig_smear.back();
 if (m_nLaphEvs>nEv_orig){
    QDPIO::cout << "Invalid number of LapH eigenvectors in SmearTune"<<endl;
    return;}
 uint nEigvecs=m_nLaphEvs;
 uint Nspin=QDP::Ns;
 TextFileWriter dataout(outdata);
 multi2d<Complex> pion(neig_smear.size(),tsepvals.size());

 for (set<int>::iterator it=source_times.begin();it!=source_times.end();++it){
    int src_time=*it;
    QDPIO::cout <<endl<<endl<< "Calculating pion correlator for source time = "<<src_time<<endl;

      // load laph eigenvectors into memory for source time
    vector<const LatticeColorVector*> laph_eigvecs;
    getLaphEigenvectors(Q,src_time,laph_eigvecs);

      // evaluate meson doublets for source time
    multi2d<Complex> meson_src_doublets;
    uint stlim=(Layout::primaryNode())?m_nLaphEvs:0;
    evaluateMesonDoublets(Q,-pz,laph_eigvecs,meson_src_doublets,0,stlim);
    clearLaphEigenvectors(Q,laph_eigvecs);

    uint tsepind=0;
    for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
       uint tsep=*tt;

       QDPIO::cout << endl<<"Beginning meson correlator calculations for time separation "<<tsep<<endl<<endl;
       int snk_time=(src_time+tsep)%Textent;

          // load laph eigenvectors into memory for sink time
       getLaphEigenvectors(Q,snk_time,laph_eigvecs);
          // evaluate meson doublets for sink time
       multi2d<Complex> meson_snk_doublets;
       evaluateMesonDoublets(Q,pz,laph_eigvecs,meson_snk_doublets,0,stlim);
       clearLaphEigenvectors(Q,laph_eigvecs);

          // load perambulators into memory (only on primary node)
       QDPIO::cout << "Loading perambulators into memory"<<endl;
       StopWatch casio; casio.start();
       multi2d<multi2d<Complex>> PQ(Nspin,Nspin);
       for (int snk_spin=1;snk_spin<=Nspin;snk_spin++)   
       for (int src_spin=1;src_spin<=Nspin;src_spin++){
          multi2d<Complex> buff(Q.getData(snk_time,snk_spin,src_time,src_spin,nEigvecs));
          if (Layout::primaryNode()){
             PQ(snk_spin-1,src_spin-1)=buff;}}
       casio.stop();
       QDPIO::cout << "Perambulators loaded in time "<<casio.getTimeInSeconds()<<" seconds"<<endl;
       QDPIO::cout << "Beginning computations..."<<endl;

          // evaluate the correlators (only on primary node)
          
       if (Layout::primaryNode()){
          for (int j=0;j<neig_smear.size();++j){
             uint Nev=neig_smear[j];
             QDPIO::cout << "Evaluating correlators for Nev = "<<Nev<<endl;
             if ((Nev>0)&&(Nev<=nEigvecs)){
                Complex corrt=zero;
                for (int snk_spin=0;snk_spin<Nspin;snk_spin++)   
                for (int src_spin=0;src_spin<Nspin;src_spin++)
                for (int n1=0;n1<Nev;++n1)
                for (int n2=0;n2<Nev;++n2){
                   Complex resA=zero, resB=zero;
                   for (int nn=0;nn<Nev;++nn){
                      resA+=meson_snk_doublets(n1,nn)*PQ(snk_spin,src_spin)(nn,n2);
                      resB+=meson_src_doublets(n2,nn)*conj(PQ(snk_spin,src_spin)(n1,nn));}
                   corrt+=resA*resB;}
                 QDPIO::cout << "result = "<<corrt<<endl;
                 pion(j,tsepind)=corrt;
             }}}
       QMP_barrier();
       }

    if (Layout::primaryNode()){
       for (int j=0;j<neig_smear.size();++j){
          dataout << "Pion correlator for Z-momentum "<<pz<<" and for Nev = "<<neig_smear[j]
                  <<" and source time "<<src_time<<"\n";
          uint tsepind=0;
          for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
             uint tsep=*tt;
             ostringstream oss;
             oss << std::setprecision(12)<<pion(j,tsepind).elem().elem().elem().real();
             dataout << "pion("<<tsep<<") = "<<oss.str()<<"\n";}}}
    }
 dataout.close();
 bulova.stop();
 QDPIO::cout << "Pion correlators computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
}


      // This version computes the different Nev values all at once.
      // Should be faster when using a large number of Nev values.
      // It distributes the final correlator computation among the nodes as well.

void LaphSmearTuneInlineMeas::calcPionCorrelatorA(PerambulatorHandler& Q, int pz,
                                 const vector<uint>& neig_smear, const set<int>& source_times,
                                 const string& outdata, const set<int>& tsepvals)
{
 StopWatch bulova;
 bulova.start();
 int Textent=Q.getGaugeConfigurationInfo().getTimeExtent();
 uint nEv_orig=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 m_nLaphEvs=neig_smear.back();
 if (m_nLaphEvs>nEv_orig){
    QDPIO::cout << "Invalid number of LapH eigenvectors in SmearTune"<<endl;
    return;}
 uint nEigvecs=m_nLaphEvs;
 uint Nspin=QDP::Ns;
 TextFileWriter dataout(outdata);
 multi2d<Complex> pion(neig_smear.size(),tsepvals.size());
 uint nsmear=neig_smear.size();

     // get the n1,n2 ranges that this node will handle
 int thisnode=Layout::nodeNumber();
 uint n1min,n1lim,n2min,n2lim;
 getMesonEvIndexLimits(nEigvecs,thisnode,n1min,n1lim,n2min,n2lim);
 
 for (set<int>::iterator it=source_times.begin();it!=source_times.end();++it){
    int src_time=*it;
    QDPIO::cout <<endl<<endl<< "Calculating pion correlator for source time = "<<src_time<<endl;

      // load laph eigenvectors into memory for source time
    vector<const LatticeColorVector*> laph_eigvecs;
    getLaphEigenvectors(Q,src_time,laph_eigvecs);

      // evaluate meson doublets for source time
    multi2d<Complex> meson_src_doublets;
    evaluateMesonDoublets(Q,-pz,laph_eigvecs,meson_src_doublets,n2min,n2lim);
    clearLaphEigenvectors(Q,laph_eigvecs);

    uint tsepind=0;
    for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
       uint tsep=*tt;

       QDPIO::cout << endl<<"Beginning meson correlator calculations for time separation "<<tsep<<endl<<endl;
       int snk_time=(src_time+tsep)%Textent;

          // load laph eigenvectors into memory for sink time
       getLaphEigenvectors(Q,snk_time,laph_eigvecs);
          // evaluate meson doublets for sink time
       multi2d<Complex> meson_snk_doublets;
       evaluateMesonDoublets(Q,pz,laph_eigvecs,meson_snk_doublets,n1min,n1lim);
       clearLaphEigenvectors(Q,laph_eigvecs);

          // load perambulators into memory (on all nodes!)
       QDPIO::cout << "Loading perambulators into memory"<<endl;
       StopWatch casio; casio.start();
       multi2d<multi2d<Complex>> PQ(Nspin,Nspin);
       for (int snk_spin=1;snk_spin<=Nspin;snk_spin++)   
       for (int src_spin=1;src_spin<=Nspin;src_spin++){
          multi2d<Complex> buff(Q.getData(snk_time,snk_spin,src_time,src_spin,nEigvecs));
          PQ(snk_spin-1,src_spin-1)=buff;}
       casio.stop();
       QDPIO::cout << "Perambulators loaded in time "<<casio.getTimeInSeconds()<<" seconds"<<endl;
       QDPIO::cout << "Beginning computations..."<<endl;
       QMP_barrier();
 
          // evaluate the correlators
       multi1d<Complex> corrt(nsmear); corrt=zero;
       for (triLooper tl(neig_smear,n1min,n1lim,n2min,n2lim);tl.notdone;++tl){
          tl.doContractTerm(corrt,meson_snk_doublets,meson_src_doublets,PQ);}
       QMP_barrier();
       QDPInternal::globalSumArray(corrt);
       for (uint s=0;s<nsmear;++s){
          QDPIO::cout << "result["<<s<<"] = "<<corrt[s]<<endl;
          pion(s,tsepind)=corrt[s];}
       }

    if (Layout::primaryNode()){
       for (int j=0;j<nsmear;++j){
          dataout << "Pion correlator for Z-momentum "<<pz<<" and for Nev = "<<neig_smear[j]
                  <<" and source time "<<src_time<<"\n";
          uint tsepind=0;
          for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
             uint tsep=*tt;
             ostringstream oss;
             oss << std::setprecision(12)<<pion(j,tsepind).elem().elem().elem().real();
             dataout << "pion("<<tsep<<") = "<<oss.str()<<"\n";}}}
    }
 dataout.close();
 bulova.stop();
 QDPIO::cout << "Pion correlators computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
}

  // ************************************************************


   //  NOTE: below is not very optimized code since these routines are not meant
   //  to be run many many times.  Should be run just a few times when determining
   //  the smearing.
   
   //  For the nucleon correlator, the nucleon operator is taken to have the form
   //
   //       Phi[uud]_{121} - Phi[uud]_{112}
   //
   //  The nucleon correlator is given in terms of the baryon triplet B(k,l,m)(t).
   //  PQ(s1,s2)(k,k')(t|t0) are the quark perambulators.  The times below are
   //  suppressed, and s1,s2 are zero-offset spin values.
   //
   //    C(t-t0) =   BB(k',l,m)  PQ(0,0)(k',k) B(k,l,m)(t0)^* 
   //
   //  where  BB(k,l,m) =  B(k,l',m')(t)  [  PQ(0,0)(l',l) PQ(1,1)(m',m)
   //
   //             - PQ(0,1)(l',l) PQ(1,0)(m',m)  ]
   //
   //  "neig_smear" contains the Nev values to get results for, and should be sorted in
   //  ascending order.
   
         //  This version is for debugging purposes only.  Baryon triplets
         //  are computed using all MPI processes, but then the final
         //  correlators in terms of the triplets are done only by the primary process.

void LaphSmearTuneInlineMeas::calcNucleonCorrelator(PerambulatorHandler& Q, int pz,
                                 const vector<uint>& neig_smear, const set<int>& source_times,
                                 const string& outdata, const set<int>& tsepvals)
{
 StopWatch bulova;
 bulova.start();
 int Textent=Q.getGaugeConfigurationInfo().getTimeExtent();
 uint nEv_orig=Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 m_nLaphEvs=neig_smear.back();
 if (m_nLaphEvs>nEv_orig){
    QDPIO::cout << "Invalid number of LapH eigenvectors in SmearTune"<<endl;
    return;}
 uint nEigvecs=m_nLaphEvs;
 uint Nspin=QDP::Ns/2;
 TextFileWriter dataout(outdata);
 multi2d<Complex> nucleon(neig_smear.size(),tsepvals.size());

 for (set<int>::iterator it=source_times.begin();it!=source_times.end();++it){
    int src_time=*it;
    QDPIO::cout <<endl<<endl<< "Calculating nucleon correlator for source time = "<<src_time<<endl;

      // load laph eigenvectors into memory for source time
    vector<const LatticeColorVector*> laph_eigvecs;
    getLaphEigenvectors(Q,src_time,laph_eigvecs);

      // evaluate baryon triplets for source time
    multi1d<multi2d<Complex>> baryon_src_triplets;
    evaluateBaryonTriplets(Q,pz,laph_eigvecs,baryon_src_triplets);
    clearLaphEigenvectors(Q,laph_eigvecs);

    uint tsepind=0;
    for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
       uint tsep=*tt;

       QDPIO::cout << endl<<"Beginning baryon correlator calculations for time separation "<<tsep<<endl<<endl;
       int snk_time=(src_time+tsep)%Textent;

          // load laph eigenvectors into memory for sink time
       getLaphEigenvectors(Q,snk_time,laph_eigvecs);
          // evaluate baryon triplets for sink time
       multi1d<multi2d<Complex>> baryon_snk_triplets;
       evaluateBaryonTriplets(Q,pz,laph_eigvecs,baryon_snk_triplets);
       clearLaphEigenvectors(Q,laph_eigvecs);

          // load perambulators into memory (only on primary node)
       QDPIO::cout << "Loading perambulators into memory"<<endl;
       StopWatch rolex; rolex.start();
       multi2d<multi2d<Complex>> PQ(Nspin,Nspin);
       for (int snk_spin=1;snk_spin<=Nspin;snk_spin++)   
       for (int src_spin=1;src_spin<=Nspin;src_spin++){
          multi2d<Complex> buff(Q.getData(snk_time,snk_spin,src_time,src_spin,nEigvecs));
          if (Layout::primaryNode()){
             PQ(snk_spin-1,src_spin-1)=buff;}}
       rolex.stop();
       QDPIO::cout << "Perambulators read in time "<<rolex.getTimeInSeconds()<<" seconds"<<endl;
       QDPIO::cout << "Beginning computations..."<<endl;
       
          // evaluate the correlators (only on primary node)
          
       if (Layout::primaryNode()){
          for (int j=0;j<neig_smear.size();++j){
             uint Nev=neig_smear[j];
             QDPIO::cout << "Evaluating correlators for Nev = "<<Nev<<endl;
             if ((Nev>0)&&(Nev<=nEigvecs)){

                multi3d<Complex> B1(Nev,Nev,Nev);
                for (int k=0;k<Nev;++k)
                for (int l=0;l<Nev;++l)
                for (int m=0;m<Nev;++m){
                   Complex res=zero;
                   for (int q=0;q<Nev;++q){
                      res+=baryon_snk_triplets[k](l,q)*PQ(1,1)(q,m);}
                   B1(k,l,m)=res;}
                multi3d<Complex> BB(Nev,Nev,Nev);
                for (int k=0;k<Nev;++k)
                for (int l=0;l<Nev;++l)
                for (int m=0;m<Nev;++m){
                   Complex res=zero;
                   for (int q=0;q<Nev;++q){
                      res+=B1(k,q,m)*PQ(0,0)(q,l);}
                   BB(k,l,m)=res;}
                for (int k=0;k<Nev;++k)
                for (int l=0;l<Nev;++l)
                for (int m=0;m<Nev;++m){
                   Complex res=zero;
                   for (int q=0;q<Nev;++q){
                      res+=baryon_snk_triplets[k](l,q)*PQ(1,0)(q,m);}
                   B1(k,l,m)=res;}
                for (int k=0;k<Nev;++k)
                for (int l=0;l<Nev;++l)
                for (int m=0;m<Nev;++m){
                   Complex res=zero;
                   for (int q=0;q<Nev;++q){
                      res+=B1(k,q,m)*PQ(0,1)(q,l);}
                   BB(k,l,m)-=res;}
                Complex corrt=zero;
                for (int k=0;k<Nev;++k)
                for (int l=0;l<Nev;++l)
                for (int m=0;m<Nev;++m){
                   Complex res=zero;
                   for (int q=0;q<Nev;++q){
                      res+=BB(q,l,m)*(PQ(0,0)(q,k)+PQ(1,1)(q,k));}  // average over row 1 and row 2
                   corrt+=res*conj(baryon_src_triplets[k](l,m));}
                QDPIO::cout << "result = "<<1.5*corrt<<endl;   // same scale as inline_scale
                nucleon(j,tsepind)=1.5*corrt;
             }}}
       QMP_barrier();
       }
    if (Layout::primaryNode()){
       for (int j=0;j<neig_smear.size();++j){
          dataout << "Nucleon correlator for Z-momentum "<<pz<<" and for Nev = "<<neig_smear[j]
                  <<" and source time "<<src_time<<"\n";
          uint tsepind=0;
          for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
             uint tsep=*tt;
             ostringstream oss;
             oss << std::setprecision(12)<<nucleon(j,tsepind).elem().elem().elem().real();
             dataout << "nucleon("<<tsep<<") = "<<oss.str()<<"\n";}}}
    }
 dataout.close();
 bulova.stop();
 QDPIO::cout << "Nucleon correlators computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
}



  // **************************************************************************************
  

LaphSmearTuneInlineMeas::triLooper::triLooper(const std::vector<uint>& limits, uint ind1min, uint ind1lim, 
                                              uint ind2min, uint ind2lim) : m_limits(limits)
{
 Nspin=QDP::Ns;
 notdone=false;
 if (m_limits.empty()) return;
 lim1=ind1lim; lim2=ind2lim;
 min1=ind1min; min2=ind2min;
 lim3=m_limits.back();
 if (lim3<lim1) lim1=lim3;
 if (lim3<lim2) lim2=lim3;
 if ((lim1<=min1)||(lim2<=min2)||(lim3==0)) return;
 i1=min1;
 i2=min2;
 i3=0;
 uint imax=max(max(i1,i2),i3);
 indxmin=0;
 while (m_limits[indxmin]<=imax) indxmin++;
 nsmear=m_limits.size();
 sums1.resize(nsmear,Nspin*Nspin);
 sums2.resize(nsmear,Nspin*Nspin);
 sums1=zero;
 sums2=zero;
 notdone=true;
}

LaphSmearTuneInlineMeas::triLooper& LaphSmearTuneInlineMeas::triLooper::operator++()
{
 i3++;
 if (i3==lim3){
    i3=0;
    i2++;
    if (i2==lim2){
       i2=min2;
       i1++;
       if (i1==lim1){
          i1=min1;
          notdone=false;}}}
 uint imax=max(max(i1,i2),i3);
 indxmin=0;
 while (m_limits[indxmin]<=imax) indxmin++;
 return *this;
}


LaphSmearTuneInlineMeas::triLooper& LaphSmearTuneInlineMeas::triLooper::operator++(int i)
{
 return operator++();
}
 

LaphSmearTuneInlineMeas::triLooper& LaphSmearTuneInlineMeas::triLooper::doContractTerm(
                               multi1d<Complex>& corrt, const multi2d<Complex>& meson_snk_doublets,  
                               const multi2d<Complex>& meson_src_doublets, 
                               const multi2d<multi2d<Complex>>& PQ)
{
 uint count=0;
 for (int snk_spin=0;snk_spin<Nspin;snk_spin++)   
 for (int src_spin=0;src_spin<Nspin;src_spin++){
    Complex res1=meson_snk_doublets(i1-min1,i3)*PQ(snk_spin,src_spin)(i3,i2);
    Complex res2=meson_src_doublets(i2-min2,i3)*conj(PQ(snk_spin,src_spin)(i1,i3));
    for (uint s=indxmin;s<nsmear;++s){
       sums1(s,count)+=res1;
       sums2(s,count)+=res2;}
    count++;}
 if (i3==(lim3-1)){
    uint imax=max(i1,i2);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<nsmear;++s){
       Complex res=zero;
       for (uint count=0;count<Nspin*Nspin;++count){
          res+=sums1(s,count)*sums2(s,count);}
       corrt[s]+=res;}
    sums1=zero;
    sums2=zero;}
 return *this;
}


// *****************************************************************************


LaphSmearTuneInlineMeas::quadLooper::quadLooper(const std::vector<uint>& limits, uint ind1min, uint ind1lim)
                                                  : m_limits(limits)
{
 initialize(ind1min,ind1lim,0,limits.back());
}


LaphSmearTuneInlineMeas::quadLooper::quadLooper(const std::vector<uint>& limits, uint ind1min, uint ind1lim, 
                                                uint ind4min, uint ind4lim)    : m_limits(limits)
{
 initialize(ind1min,ind1lim,ind4min,ind4lim);
}


void LaphSmearTuneInlineMeas::quadLooper::initialize(uint ind1min, uint ind1lim, uint ind4min, uint ind4lim)
{
 notdone=false;
 if (m_limits.empty()) return;
 lim23=m_limits.back();
 if (lim23==0) return;
 if ((ind1min>=lim23)||(ind4min>=lim23)) return;
 min1=ind1min;
 lim1=ind1lim;
 if (lim1<=min1) return;
 min4=ind4min;
 lim4=ind4lim;
 if (lim4<=min4) return;
 i1=min1; o1=0;
 i2=0; 
 i3=0; 
 i4=min4; o4=0;
 uint imax=max(max(i1,i2),max(i3,i4));
 indxmin=0;
 while (m_limits[indxmin]<=imax) indxmin++;
 nsmear=m_limits.size();
 sums.resize(nsmear);
 sums=zero;
 notdone=true;
}

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::operator++()
{
 i4++; o4++;
 if (i4==lim4){
    i4=min4; o4=0;
    i3++;
    if (i3==lim23){
       i3=0;
       i2++;
       if (i2==lim23){
          i2=0;
          i1++; o1++;
          if (i1==lim1){
             i1=min1; o1=0;
             notdone=false;}}}}
 uint imax=max(max(i1,i2),max(i3,i4));
 indxmin=0;
 while (m_limits[indxmin]<=imax) indxmin++;
 return *this;
}


LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::operator++(int i)
{
 return operator++();
}
 

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doTerm1(multi1d<multi3d<Complex>>& B,        
                                          const multi1d<multi2d<Complex>>& baryon_snk_triplets,  
                                          const multi2d<multi2d<Complex>>& PQ)
{
 for (uint s=indxmin;s<sums.size();++s){
    sums[s]+=baryon_snk_triplets[o1](i2,i4)*PQ(1,1)(i4,i3);}
 if (i4==(lim4-1)){
    uint imax=max(max(i1,i2),i3);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<sums.size();++s){
       B[s](o1,i2,i3)=sums[s];}
    sums=zero;}
 return *this;
}

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doTerm2(multi1d<multi3d<Complex>>& B, 
                                          const multi1d<multi3d<Complex>>& B1, const multi2d<multi2d<Complex>>& PQ)
{
 for (uint s=indxmin;s<sums.size();++s){
    sums[s]+=B1[s](o1,i4,i3)*PQ(0,0)(i4,i2);}
 if (i4==(lim4-1)){
    uint imax=max(max(i1,i2),i3);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<sums.size();++s){
       B[s](o1,i2,i3)=sums[s];}
    sums=zero;}
 return *this;
}

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doTerm3(multi1d<multi3d<Complex>>& B,        
                                          const multi1d<multi2d<Complex>>& baryon_snk_triplets,  
                                          const multi2d<multi2d<Complex>>& PQ)
{
 for (uint s=indxmin;s<sums.size();++s){
    sums[s]+=baryon_snk_triplets[o1](i2,i4)*PQ(1,0)(i4,i3);}
 if (i4==(lim4-1)){
    uint imax=max(max(i1,i2),i3);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<sums.size();++s){
       B[s](o1,i2,i3)=sums[s];}
    sums=zero;}
 return *this;
}

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doTerm4(multi1d<multi3d<Complex>>& B, 
                                          const multi1d<multi3d<Complex>>& B1, const multi2d<multi2d<Complex>>& PQ)
{
 for (uint s=indxmin;s<sums.size();++s){
    sums[s]+=B1[s](o1,i4,i3)*PQ(0,1)(i4,i2);}
 if (i4==(lim4-1)){
    uint imax=max(max(i1,i2),i3);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<sums.size();++s){
       B[s](o1,i2,i3)-=sums[s];}
    sums=zero;}
 return *this;
}

LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doSrcTerm(multi1d<multi3d<Complex>>& BS,        
                                          const multi1d<multi2d<Complex>>& baryon_src_triplets,  
                                          const multi2d<multi2d<Complex>>& PQ)
{
 for (uint s=indxmin;s<sums.size();++s){
    sums[s]+=(PQ(0,0)(i1,i4)+PQ(1,1)(i1,i4))*conj(baryon_src_triplets[o4](i2,i3));}  // average over row 1 and row 2
 if (i4==(lim4-1)){
    uint imax=max(max(i1,i2),i3);
    uint smin=0;
    while (m_limits[smin]<=imax) smin++;
    for (uint s=smin;s<sums.size();++s){
       BS[s](o1,i2,i3)+=sums[s];}
    sums=zero;}
 return *this;
}


LaphSmearTuneInlineMeas::quadLooper& LaphSmearTuneInlineMeas::quadLooper::doContractTerm(multi1d<Complex>& corrt,        
                                            const multi1d<multi3d<Complex>>& BB, const multi1d<multi3d<Complex>>& BS)
{
 for (uint s=indxmin;s<sums.size();++s){
    corrt[s]+=BB[s](o1,i2,i3)*BS[s](o1,i2,i3);}
 return *this;
}

 
 // **********************************************************************************************

   //  This is the important version to use.  The nucleon correlator is given by
   //
   //    C(t-t0) =   BB(k',l,m)  PQ(0,0)(k',k) B(k,l,m)(t0)^* 
   //
   //  where  BB(k,l,m) =  B(k,l',m')(t)  [  PQ(0,0)(l',l) PQ(1,1)(m',m)
   //
   //             - PQ(0,1)(l',l) PQ(1,0)(m',m)  ]
   //
   //  Baryon triplets are computed on all processes are usual, splitting up
   //  the 3D lattice across the MPI processes.  Each process maintains
   //  only one part of the baryon triplets:  ind1min..ind1lim-1 of the first index.


void LaphSmearTuneInlineMeas::calcNucleonCorrelatorBlocked(PerambulatorHandler& Q, int pz,
                                 const vector<uint>& smear_neigs, const set<int>& source_times,
                                 const string& outdata, const set<int>& tsepvals)
{
 StopWatch bulova;
 bulova.start();
 uint nEigvecs_orig = Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 uint Nspin=QDP::Ns/2;
 if (smear_neigs.back()>nEigvecs_orig){
    QDPIO::cout << "Too large an Neigs requested in calcNucleon Correlator"<<endl;
    return;}
 m_nLaphEvs=smear_neigs.back();
 uint nEv=m_nLaphEvs;
 uint nsmear=smear_neigs.size();
     //  Determine ind1min and ind1lim for this MPI rank
 int thisnode=Layout::nodeNumber();
 uint ind1min,ind1lim;
 m_ntripactive=getIndex1Limits(thisnode,nEv,ind1min,ind1lim);
 if ((m_ntripactive%2)==1){
    QDPIO::cout << "Nucleon correlator calculations requires an EVEN number of active triplets processes"<<endl;
    return;}
 multi2d<Complex> nucleon(nsmear,tsepvals.size());
 TextFileWriter dataout(outdata);
 int Textent=Q.getGaugeConfigurationInfo().getTimeExtent();
 uint smearindmin=0;
 while (smear_neigs[smearindmin]<=ind1min) smearindmin++;


 for (set<int>::iterator it=source_times.begin();it!=source_times.end();++it){
    int src_time=*it;
    QDPIO::cout <<endl<<endl<< "Calculating nucleon correlator for source time = "<<src_time<<endl;

      // load laph eigenvectors into memory for source time
    vector<const LatticeColorVector*> laph_eigvecs;
    getLaphEigenvectors(Q,src_time,laph_eigvecs);

      // evaluate baryon triplets for source time
    multi1d<multi2d<Complex>> baryon_src_triplets;
    evaluateBaryonTripletsBlocked(Q,pz,ind1min,ind1lim,laph_eigvecs,baryon_src_triplets);
    clearLaphEigenvectors(Q,laph_eigvecs);

  //  for (int tsep=tsepmin;tsep<=tsepmax;++tsep){
    uint tsepind=0;
    for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
       uint tsep=*tt;

       QDPIO::cout << endl<<"Beginning baryon correlator calculations for time separation "<<tsep<<endl<<endl;
       int snk_time=(src_time+tsep)%Textent;
       multi1d<Complex> corrt(nsmear); corrt=zero;

          // load laph eigenvectors into memory for sink time
       getLaphEigenvectors(Q,snk_time,laph_eigvecs);
          // evaluate baryon triplets for sink time
       multi1d<multi2d<Complex>> baryon_snk_triplets;
       evaluateBaryonTripletsBlocked(Q,pz,ind1min,ind1lim,laph_eigvecs,baryon_snk_triplets);
       clearLaphEigenvectors(Q,laph_eigvecs);

          // load perambulators into memory
       QDPIO::cout << "Loading perambulators into memory"<<endl;
       StopWatch rolex; rolex.start();
       multi2d<multi2d<Complex>> PQ(Nspin,Nspin);
       for (int snk_spin=1;snk_spin<=Nspin;snk_spin++)   
       for (int src_spin=1;src_spin<=Nspin;src_spin++){
          multi2d<Complex> buff(Q.getData(snk_time,snk_spin,src_time,src_spin,m_nLaphEvs));
             PQ(snk_spin-1,src_spin-1)=buff;}
       rolex.stop();
       QDPIO::cout << "Perambulators read in time "<<rolex.getTimeInSeconds()<<" seconds"<<endl;
       QDPIO::cout << "Beginning computations..."<<endl;
       rolex.reset(); rolex.start();
          //  calculate the BB matrices
       uint ind1size=ind1lim-ind1min;
       multi1d<multi3d<Complex>> BB(nsmear);
       multi1d<multi3d<Complex>> BB1(nsmear);
       if (ind1lim>ind1min){
          for (uint s=smearindmin;s<nsmear;++s){
             BB[s].resize(ind1size,smear_neigs[s],smear_neigs[s]);
             BB1[s].resize(ind1size,smear_neigs[s],smear_neigs[s]);
             BB[s]=zero;}
          for (quadLooper ql(smear_neigs,ind1min,ind1lim);ql.notdone;++ql){
             ql.doTerm1(BB1,baryon_snk_triplets,PQ);} 
          for (quadLooper ql(smear_neigs,ind1min,ind1lim);ql.notdone;++ql){
             ql.doTerm2(BB,BB1,PQ);}
          for (quadLooper ql(smear_neigs,ind1min,ind1lim);ql.notdone;++ql){
             ql.doTerm3(BB1,baryon_snk_triplets,PQ);}
          for (quadLooper ql(smear_neigs,ind1min,ind1lim);ql.notdone;++ql){
             ql.doTerm4(BB,BB1,PQ);}
          }
       QMP_barrier();

       if (ind1lim>ind1min){
          for (uint s=smearindmin;s<nsmear;++s){
             BB1[s]=zero;}}
          // evaluate the src side
          // cycle the src baryon triplets   
       for (int src_blk=0;src_blk<m_ntripactive;++src_blk){
          //StopWatch casio; casio.start();
          //QDPIO::cout << "Block "<<src_blk;
          if (src_blk>0){
             cycleBaryonTripletsBlocked(baryon_src_triplets,src_blk);}
          int srcnode=(thisnode+m_ntripactive-src_blk)%m_ntripactive;
          uint srcind1min,srcind1lim;
          getIndex1Limits(srcnode,nEv,srcind1min,srcind1lim);
          for (quadLooper ql(smear_neigs,ind1min,ind1lim,srcind1min,srcind1lim);ql.notdone;++ql){
             ql.doSrcTerm(BB1,baryon_src_triplets,PQ);}
          //casio.stop();
          / *QDPIO::cout << " completed in time "<<casio.getTimeInSeconds()<<" seconds "<<endl;* /}
       cycleBaryonTripletsBlocked(baryon_src_triplets,m_ntripactive);  // one last cycle to put back into original locations
           
          // now do the final contraction
       if (ind1lim>ind1min){
          for (quadLooper ql(smear_neigs,ind1min,ind1lim,0,1);ql.notdone;++ql){
             ql.doContractTerm(corrt,BB,BB1);}}
       QMP_barrier();
       QDPInternal::globalSumArray(corrt);
       rolex.stop();
       QDPIO::cout << "done: time "<<rolex.getTimeInSeconds()<<" seconds"<<endl;

       for (uint s=0;s<nsmear;++s){
          QDPIO::cout << "result["<<s<<"] = "<<1.5*corrt[s]<<endl;   // same scale as inline_scale
          nucleon(s,tsepind)=1.5*corrt[s];}}

    if (Layout::primaryNode()){
       for (int j=0;j<nsmear;++j){
          dataout << "Nucleon correlator for Z-momentum "<<pz<<" and for Nev = "<<smear_neigs[j]
                  <<" and source time "<<src_time<<"\n";
          uint tsepind=0;
          for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
             uint tsep=*tt;
             ostringstream oss;
             oss << std::setprecision(12)<<nucleon(j,tsepind).elem().elem().elem().real();
             dataout << "nucleon("<<tsep<<") = "<<oss.str()<<"\n";}}}
    }
 dataout.close();
 bulova.stop();
 QDPIO::cout << "Nucleon correlators computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
}



  //  This version uses less memory.  Does the different Nev values in sequence.

void LaphSmearTuneInlineMeas::calcNucleonCorrelatorBlocked2(PerambulatorHandler& Q, int pz,
                                 const vector<uint>& smear_neigs, const set<int>& source_times,
                                 const string& outdata, const set<int>& tsepvals)
{
 StopWatch bulova;
 bulova.start();
 uint nEigvecs_orig = Q.getQuarkSmearingInfo().getNumberOfLaplacianEigenvectors();
 uint Nspin=QDP::Ns/2;
 if (smear_neigs.back()>nEigvecs_orig){
    QDPIO::cout << "Too large an Neigs requested in calcNucleon Correlator"<<endl;
    return;}
 m_nLaphEvs=smear_neigs.back();
 uint nEv=m_nLaphEvs;
 uint nsmear=smear_neigs.size();
     //  Determine ind1min and ind1lim for this MPI rank
 int thisnode=Layout::nodeNumber();
 uint ind1min,ind1lim;
 m_ntripactive=getIndex1Limits(thisnode,nEv,ind1min,ind1lim);
 if ((m_ntripactive%2)==1){
    QDPIO::cout << "Nucleon correlator calculations requires an EVEN number of active triplets processes"<<endl;
    return;}
 multi2d<Complex> nucleon(nsmear,tsepvals.size());
 TextFileWriter dataout(outdata);
 int Textent=Q.getGaugeConfigurationInfo().getTimeExtent();


 for (set<int>::iterator it=source_times.begin();it!=source_times.end();++it){
    int src_time=*it;
    QDPIO::cout <<endl<<endl<< "Calculating nucleon correlator for source time = "<<src_time<<endl;

      // load laph eigenvectors into memory for source time
    vector<const LatticeColorVector*> laph_eigvecs;
    getLaphEigenvectors(Q,src_time,laph_eigvecs);

      // evaluate baryon triplets for source time
    multi1d<multi2d<Complex>> baryon_src_triplets;
    evaluateBaryonTripletsBlocked(Q,pz,ind1min,ind1lim,laph_eigvecs,baryon_src_triplets);
    clearLaphEigenvectors(Q,laph_eigvecs);

    uint tsepind=0;
    for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
       uint tsep=*tt;

       QDPIO::cout << endl<<"Beginning baryon correlator calculations for time separation "<<tsep<<endl<<endl;
       int snk_time=(src_time+tsep)%Textent;

          // load laph eigenvectors into memory for sink time
       getLaphEigenvectors(Q,snk_time,laph_eigvecs);
          // evaluate baryon triplets for sink time
       multi1d<multi2d<Complex>> baryon_snk_triplets;
       evaluateBaryonTripletsBlocked(Q,pz,ind1min,ind1lim,laph_eigvecs,baryon_snk_triplets);
       clearLaphEigenvectors(Q,laph_eigvecs);

          // load perambulators into memory
       QDPIO::cout << "Loading perambulators into memory"<<endl;
       StopWatch rolex; rolex.start();
       multi2d<multi2d<Complex>> PQ(Nspin,Nspin);
       for (int snk_spin=1;snk_spin<=Nspin;snk_spin++)   
       for (int src_spin=1;src_spin<=Nspin;src_spin++){
          multi2d<Complex> buff(Q.getData(snk_time,snk_spin,src_time,src_spin,m_nLaphEvs));
             PQ(snk_spin-1,src_spin-1)=buff;}
       rolex.stop();
       QDPIO::cout << "Perambulators read in time "<<rolex.getTimeInSeconds()<<" seconds"<<endl;
       QDPIO::cout << "Beginning computations..."<<endl;
       rolex.reset(); 

          // loop sequentially over the Nev values
       for (uint s=0;s<nsmear;++s){
       
          uint Nev=smear_neigs[s];
          QDPIO::cout << "Starting computation for Nev = "<<Nev<<"  s = "<<s<<endl;
          int ind1size=int(min(ind1lim,Nev))-int(ind1min);
          multi3d<Complex> BB,BB1;
          rolex.reset(); rolex.start();
          //  calculate the BB matrices
          if (ind1size>0){
             BB.resize(ind1size,Nev,Nev);
             BB1.resize(ind1size,Nev,Nev);
             BB=zero;
             for (int k=0;k<ind1size;++k)
             for (int l=0;l<Nev;++l)
             for (int m=0;m<Nev;++m){
                Complex res=zero;
                for (int q=0;q<Nev;++q){
                   res+=baryon_snk_triplets[k](l,q)*PQ(1,1)(q,m);}
                BB1(k,l,m)=res;}
             for (int k=0;k<ind1size;++k)
             for (int l=0;l<Nev;++l)
             for (int m=0;m<Nev;++m){
                Complex res=zero;
                for (int q=0;q<Nev;++q){
                   res+=BB1(k,q,m)*PQ(0,0)(q,l);}
                BB(k,l,m)=res;}
             for (int k=0;k<ind1size;++k)
             for (int l=0;l<Nev;++l)
             for (int m=0;m<Nev;++m){
                Complex res=zero;
                for (int q=0;q<Nev;++q){
                   res+=baryon_snk_triplets[k](l,q)*PQ(1,0)(q,m);}
                BB1(k,l,m)=res;}
             for (int k=0;k<ind1size;++k)
             for (int l=0;l<Nev;++l)
             for (int m=0;m<Nev;++m){
                Complex res=zero;
                for (int q=0;q<Nev;++q){
                   res+=BB1(k,q,m)*PQ(0,1)(q,l);}
                BB(k,l,m)-=res;}
             }
          QMP_barrier();

          if (ind1size>0){
                BB1=zero;}
          // evaluate the src side
          // cycle the src baryon triplets   
          for (int src_blk=0;src_blk<m_ntripactive;++src_blk){
             //StopWatch casio; casio.start();
             //QDPIO::cout << "Block "<<src_blk;
             if (src_blk>0){
                cycleBaryonTripletsBlocked(baryon_src_triplets,src_blk);}
             int srcnode=(thisnode+m_ntripactive-src_blk)%m_ntripactive;
             uint srcind1min,srcind1lim;
             getIndex1Limits(srcnode,nEv,srcind1min,srcind1lim);
             int srcsize=int(min(srcind1lim,Nev))-int(srcind1min);
             if (srcsize>0){
                for (int o1=0;o1<ind1size;++o1){
                   int i1=o1+ind1min;
                   for (int i2=0;i2<Nev;++i2)
                   for (int i3=0;i3<Nev;++i3){
                      Complex sum=zero;
                      for (int o4=0;o4<srcsize;++o4){
                        int i4=o4+srcind1min;
                        sum+=(PQ(0,0)(i1,i4)+PQ(1,1)(i1,i4))*conj(baryon_src_triplets[o4](i2,i3));}  // average over row 1 and row 2
                      BB1(o1,i2,i3)+=sum;}}}
             //casio.stop();
             / *QDPIO::cout << " completed in time "<<casio.getTimeInSeconds()<<" seconds "<<endl;* /}
          cycleBaryonTripletsBlocked(baryon_src_triplets,m_ntripactive);  // one last cycle to put back into original locations
           
          // now do the final contraction
          Complex corrt=zero;
          if (ind1size>0){
             for (int k=0;k<ind1size;++k)
             for (int l=0;l<Nev;++l)
             for (int m=0;m<Nev;++m){
                corrt+=BB(k,l,m)*BB1(k,l,m);}}
          QMP_barrier();
          QDPInternal::globalSum(corrt);
          rolex.stop();
          QDPIO::cout << "done: time "<<rolex.getTimeInSeconds()<<" seconds"<<endl;

          QDPIO::cout << "result["<<s<<"] = "<<1.5*corrt<<endl;   // same scale as inline_scale
          nucleon(s,tsepind)=1.5*corrt;}}

    if (Layout::primaryNode()){
       for (int j=0;j<nsmear;++j){
          dataout << "Nucleon correlator for Z-momentum "<<pz<<" and for Nev = "<<smear_neigs[j]
                  <<" and source time "<<src_time<<"\n";
          uint tsepind=0;
          for (set<int>::const_iterator tt=tsepvals.begin();tt!=tsepvals.end();++tt,++tsepind){
             uint tsep=*tt;
             ostringstream oss;
             oss << std::setprecision(12)<<nucleon(j,tsepind).elem().elem().elem().real();
             dataout << "nucleon("<<tsep<<") = "<<oss.str()<<"\n";}}}
    }
 dataout.close();
 bulova.stop();
 QDPIO::cout << "Nucleon correlators computed: time = "<<bulova.getTimeInSeconds()<<" seconds"<<endl; 
}

*/

// ******************************************************************
}
