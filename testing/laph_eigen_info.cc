#include "laph_eigen_info.h"
#include "layout_info.h"
#include "quda_info.h"

using namespace std;

namespace LaphEnv {


// *************************************************************

   // XmlReader constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(const XMLHandler& xmlin)
{
 XMLHandler xml_in(xmlin);
 xml_tag_assert(xml_in,"LaphEigenSolverInfo");
 XMLHandler xmlr(xml_in, "LaphEigenSolverInfo");
 extract_info_from_reader(xmlr);
}


void LaphEigenSolverInfo::extract_info_from_reader(XMLHandler& xmlr)
{
 xmlread(xmlr,"ResidualTolerance", tolerance, "LaphEigenSolverInfo");
 xmlread(xmlr,"MaxIterations", maxIterations, "LaphEigenSolverInfo");
 xmlread(xmlr,"KrylovDimension", dimKrylov, "LaphEigenSolverInfo");

 chebyshevOrder=1;
 if (xml_tag_count(xmlr,"ChebyshevOrder")==1)
    xmlread(xmlr,"ChebyshevOrder", chebyshevOrder, "LaphEigenSolverInfo");
 if (chebyshevOrder<1) chebyshevOrder=1;

 maxEigenvalue=15.0;
 if (xml_tag_count(xmlr,"MaxEigenvalue")==1){
    xmlread(xmlr,"MaxEigenvalue",maxEigenvalue,"LaphEigenSolverInfo");
    if (maxEigenvalue<12.0){
       xmlreadfail(xmlr,"LaphEigenSolverInfo",
           "invalid MaxEigenvalue in LaphEigenSolverInfo: should exceed 12");}}
 else{
    if (chebyshevOrder>1){
       xmlreadfail(xmlr,"LaphEigenSolverInfo",
          "when Chebyshev acceleration used, must supply a maximum eigenvalue");}}

 cutoffEigenvalue=0.5;
 if (xml_tag_count(xmlr,"CutoffEigenvalue")==1){
    xmlread(xmlr,"CutoffEigenvalue",cutoffEigenvalue,"LaphEigenSolverInfo");
    if ((cutoffEigenvalue<0.03)||(cutoffEigenvalue>(maxEigenvalue-0.9))){
       xmlreadfail(xmlr,"LaphEigenSolverInfo",
           "invalid CutoffEigenvalue in LaphEigenSolverInfo");}}
 else{
    if (chebyshevOrder>1){
       xmlreadfail(xmlr,"LaphEigenSolverInfo",
           "when Chebyshev acceleration used, must supply a cutoff eigenvalue");}}

 startVector="equal_components";
 if (xml_tag_count(xmlr,"StartingVectorType")==1){
    string svread;
    xmlread(xmlr,"StartingVectorType",svread,"LaphEigenSolverInfo");
    svread=tidyString(svread);
    if (svread=="random") startVector="random";
    else if (svread=="equal_components") startVector="equal_components";
    else{
       xmlreadfail(xmlr,"LaphEigenSolverInfo",
           "Invalid <StartingVectorType> tag in LaphEigenSolverInfo");}}

 if (xml_tag_count(xmlr,"OutputVerbosity")==1){
    xmlread(xmlr,"OutputVerbosity",outputVerbosity,"LaphEigenSolverInfo");
    if (outputVerbosity<0) outputVerbosity=0;
    if (outputVerbosity>2) outputVerbosity=2;}
 else
    outputVerbosity=0;

 if ((maxIterations<4)||(dimKrylov<6)||(tolerance<=0)){
    xmlreadfail(xmlr,"LaphEigenSolverInfo",
        "invalid input parameters in LaphEigenSolverInfo");}
}


 // *************************************************************

    // copy constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(const LaphEigenSolverInfo& in) 
            :  maxIterations(in.maxIterations),
               dimKrylov(in.dimKrylov),
               tolerance(in.tolerance),
               chebyshevOrder(in.chebyshevOrder),
               maxEigenvalue(in.maxEigenvalue),
               cutoffEigenvalue(in.cutoffEigenvalue),
               startVector(in.startVector),
               outputVerbosity(in.outputVerbosity) {}

LaphEigenSolverInfo& LaphEigenSolverInfo::operator=(
               const LaphEigenSolverInfo& in)
{
 maxIterations=in.maxIterations;
 startVector=in.startVector;
 dimKrylov=in.dimKrylov;
 tolerance=in.tolerance;
 chebyshevOrder=in.chebyshevOrder;
 maxEigenvalue=in.maxEigenvalue;
 cutoffEigenvalue=in.cutoffEigenvalue;
 outputVerbosity=in.outputVerbosity;
 return *this;
}


string LaphEigenSolverInfo::output(int indent) const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.output(indent);
}

void LaphEigenSolverInfo::output(XMLHandler& xmlout) const
{
 xmlout.set_root("LaphEigenSolverInfo");
 xmlout.put_child("ResidualTolerance", make_string(tolerance));
 xmlout.put_child("MaxIterations", make_string(maxIterations));
 xmlout.put_child("KrylovDimension", make_string(dimKrylov));
 if (chebyshevOrder>1){
    xmlout.put_child("ChebyshevOrder", make_string(chebyshevOrder));
    xmlout.put_child("MaxEigenvalue",make_string(maxEigenvalue));
    xmlout.put_child("CutoffEigenvalue",make_string(cutoffEigenvalue));}
 xmlout.put_child("StartingVectorType", startVector);
 xmlout.put_child("OutputVerbosity", make_string(outputVerbosity));
}


void LaphEigenSolverInfo::setQudaParam(QudaInvertParam& eig_inv_param, 
                                       QudaEigParam& eig_param,
                                       const QuarkSmearingInfo& qsmear) const
{
    // build eig_inv_param first
 eig_inv_param = newQudaInvertParam();
 eig_inv_param.dslash_type = QUDA_LAPLACE_DSLASH;
 eig_inv_param.inv_type = QUDA_CG_INVERTER;
 eig_inv_param.tol = tolerance;
 eig_inv_param.maxiter = maxIterations;
 eig_inv_param.reliable_delta = 0.1;
 eig_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
 eig_inv_param.mass = -2.0;                  // for massless Laplacian
 eig_inv_param.kappa = 1.0 / (8.0 + eig_inv_param.mass);
 eig_inv_param.dagger = QUDA_DAG_NO;
 eig_inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
 eig_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
 eig_inv_param.dirac_order = QUDA_DIRAC_ORDER;

 eig_inv_param.laplace3D = LayoutInfo::Ndim-1;  // spatial laplacian
 eig_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
 eig_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
 eig_inv_param.cpu_prec = QudaInfo::get_cpu_prec();
 eig_inv_param.cuda_prec = QudaInfo::get_cuda_prec();
 if (getOutputVerbosity()==0){
    eig_inv_param.verbosity = QUDA_SILENT;}
 else if (getOutputVerbosity()==1){
    eig_inv_param.verbosity = QUDA_SUMMARIZE;}
 else{
    eig_inv_param.verbosity = QUDA_VERBOSE;}

// eig_inv_param.extlib_type = solver_ext_lib;
// eig_inv_param.native_blas_lapack = (native_blas_lapack ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE);
 eig_inv_param.struct_size = sizeof(eig_inv_param);

    // build eig_param next
 eig_param = newQudaEigParam();
 eig_param.invert_param = &eig_inv_param;
 eig_param.eig_type = QUDA_EIG_TR_LANCZOS_3D;
 eig_param.spectrum = QUDA_SPECTRUM_SR_EIG;   // -Laplacian eigenvalues all positive
 eig_param.n_conv = qsmear.getNumberOfLaplacianEigenvectors();
 eig_param.max_restarts = maxIterations;
 eig_param.max_ortho_attempts = 6;
 eig_param.tol = tolerance;
 eig_param.n_kr = dimKrylov;
 eig_param.n_ev = eig_param.n_conv;
 eig_param.use_poly_acc = (chebyshevOrder>1) ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
 eig_param.poly_deg = chebyshevOrder;
 eig_param.a_min = cutoffEigenvalue/6.0;   // due to quda operator definition
 eig_param.a_max = maxEigenvalue/6.0;      // due to quda operator definition
 eig_param.n_ev_deflate = 0;


// eig_param.ortho_block_size = eig_ortho_block_size;
// eig_param.block_size
//   = (eig_param.eig_type == QUDA_EIG_TR_LANCZOS || eig_param.eig_type == QUDA_EIG_IR_ARNOLDI) ? 1 : eig_block_size;
// eig_param.qr_tol = eig_qr_tol;
// eig_param.batched_rotate = eig_batched_rotate;
 eig_param.require_convergence = QUDA_BOOLEAN_TRUE;
// eig_param.check_interval = eig_check_interval;

 eig_param.use_norm_op = QUDA_BOOLEAN_FALSE;
 eig_param.use_dagger = QUDA_BOOLEAN_FALSE;
 eig_param.use_pc = QUDA_BOOLEAN_FALSE;
 eig_param.compute_gamma5 = QUDA_BOOLEAN_FALSE;
 eig_param.compute_svd = QUDA_BOOLEAN_FALSE;
 eig_param.use_eigen_qr = QUDA_BOOLEAN_FALSE;
 eig_param.ortho_dim = LayoutInfo::Ndim-1;
 eig_param.ortho_dim_size_local=LayoutInfo::getRankLattSizes()[LayoutInfo::Ndim-1];
 eig_param.struct_size = sizeof(eig_param);
}


/*

  // Parameter set for solving eigenvalue problems.
  typedef struct QudaEigParam_s {
    / ** Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size * /
    size_t struct_size;

    // EIGENSOLVER PARAMS
    //-------------------------------------------------
    / ** Used to store information pertinent to the operator ** /
    QudaInvertParam *invert_param;

    / ** Type of eigensolver algorithm to employ ** /
    QudaEigType eig_type;

    / ** Use Polynomial Acceleration ** /
    QudaBoolean use_poly_acc;

    / ** Degree of the Chebysev polynomial ** /
    int poly_deg;

    / ** Range used in polynomial acceleration ** /
    double a_min;
    double a_max;

    / ** Whether to preserve the deflation space between solves.  If
        true, the space will be stored in an instance of the
        deflation_space struct, pointed to by preserve_deflation_space * /
    QudaBoolean preserve_deflation;

    / ** This is where we store the deflation space.  This will point
        to an instance of deflation_space. When a deflated solver is enabled, the deflation space will be obtained from this.  * /
    void *preserve_deflation_space;

    / ** If we restore the deflation space, this boolean indicates
        whether we are also preserving the evalues or recomputing
        them.  For example if a different mass shift is being used
        than the one used to generate the space, then this should be
        false, but preserve_deflation would be true * /
    QudaBoolean preserve_evals;

    / ** What type of Dirac operator we are using ** /
    / ** If !(use_norm_op) && !(use_dagger) use M. ** /
    / ** If use_dagger, use Mdag ** /
    / ** If use_norm_op, use MdagM ** /
    / ** If use_norm_op && use_dagger use MMdag. ** /
    / ** If use_pc for any, then use the even-odd pc version ** /
    QudaBoolean use_dagger;
    QudaBoolean use_norm_op;
    QudaBoolean use_pc;

    / ** Use Eigen routines to eigensolve the upper Hessenberg via QR ** /
    QudaBoolean use_eigen_qr;

    / ** Performs an MdagM solve, then constructs the left and right SVD. ** /
    QudaBoolean compute_svd;

    / ** Performs the \gamma_5 OP solve by Post multipling the eignvectors with
        \gamma_5 before computing the eigenvalues * /
    QudaBoolean compute_gamma5;

    / ** If true, the solver will error out if the convergence criteria are not met ** /
    QudaBoolean require_convergence;

    / ** Which part of the spectrum to solve ** /
    QudaEigSpectrumType spectrum;

    / ** Size of the eigenvector search space ** /
    int n_ev;
    / ** Total size of Krylov space ** /
    int n_kr;
    / ** Max number of locked eigenpairs (deduced at runtime) ** /
    int nLockedMax;
    / ** Number of requested converged eigenvectors ** /
    int n_conv;
    / ** Number of requested converged eigenvectors to use in deflation ** /
    int n_ev_deflate;
    / ** Tolerance on the least well known eigenvalue's residual ** /
    double tol;
    / ** Tolerance on the QR iteration ** /
    double qr_tol;
    / ** For IRLM/IRAM, check every nth restart ** /
    int check_interval;
    / ** For IRLM/IRAM, quit after n restarts ** /
    int max_restarts;
    / ** For the Ritz rotation, the maximal number of extra vectors the solver may allocate ** /
    int batched_rotate;
    / ** For block method solvers, the block size ** /
    int block_size;
    / ** For block method solvers, quit after n attempts at block orthonormalisation ** /
    int max_ortho_attempts;
    / ** For hybrid modifeld Gram-Schmidt orthonormalisations ** /
    int ortho_block_size;

    / ** In the test function, cross check the device result against ARPACK ** /
    QudaBoolean arpack_check;
    / ** For Arpack cross check, name of the Arpack logfile ** /
    char arpack_logfile[512];

    / ** Name of the QUDA logfile (residua, upper Hessenberg/tridiag matrix updates) ** /
    char QUDA_logfile[512];

    / ** The orthogonal direction in the 3D eigensolver ** /
    int ortho_dim;

    / ** The size of the orthogonal direction in the 3D eigensolver, local ** /
    int ortho_dim_size_local;

    //-------------------------------------------------

    // EIG-CG PARAMS
    //-------------------------------------------------
    int nk;
    int np;

    / ** Whether to load eigenvectors * /
    QudaBoolean import_vectors;

    / ** The precision of the Ritz vectors * /
    QudaPrecision cuda_prec_ritz;

    / ** The memory type used to keep the Ritz vectors * /
    QudaMemoryType mem_type_ritz;

    / ** Location where deflation should be done * /
    QudaFieldLocation location;

    / ** Whether to run the verification checks once set up is complete * /
    QudaBoolean run_verify;

    / ** Filename prefix where to load the null-space vectors * /
    char vec_infile[256];

    / ** Filename prefix for where to save the null-space vectors * /
    char vec_outfile[256];

    / ** The precision with which to save the vectors * /
    QudaPrecision save_prec;

    / ** Whether to inflate single-parity eigen-vector I/O to a full
        field (e.g., enabling this is required for compatability with
        MILC I/O) * /
    QudaBoolean io_parity_inflate;

    / ** Whether to save eigenvectors in QIO singlefile or partfile format * /
    QudaBoolean partfile;

    / ** The Gflops rate of the eigensolver setup * /
    double gflops;

    / **< The time taken by the eigensolver setup * /
    double secs;

    / ** Which external library to use in the deflation operations (Eigen) * /
    QudaExtLibType extlib_type;
    //-------------------------------------------------
  } QudaEigParam;*/

// *************************************************************
}
