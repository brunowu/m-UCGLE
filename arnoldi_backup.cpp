#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <complex>
#include <unistd.h>


//Block KrylovSchur METHODs to approximate the eigevalues

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp" //Anasazi interface to Tpetra
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

using Tpetra::CrsMatrix;
using Tpetra::Map;
using std::vector;


int main( int argc, char *argv[] ){
  //MPI_Init( &argc, &argv );
  Teuchos::GlobalMPISession mpisess (&argc, &argv, &std::cout);

  int arank, asize;

  int exit_type = 0;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &asize );
  MPI_Comm_rank( MPI_COMM_WORLD, &arank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(arank == 0){
    printf("Info ]> The Comm world size of ERAM is %d \n", asize);
  }

  #ifndef HAVE_TPETRA_COMPLEX_DOUBLE
  #  error "Anasazi: This test requires Scalar = std::complex<double> to be enabled in Tpetra."
  #else

  #endif

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  typedef std::complex<double>                ST;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Tpetra::MultiVector<ST>             MV;
  typedef MV::global_ordinal_type             GO; //global indices for matrices
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;


  bool success = false;

  const ST ONE = SCT::one ();

  int info = 0;

  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  const int MyPID = comm->getRank ();

  /*stantard MPI functionality*/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "Arnoldi in Trilinos said: my rank = " << rank << "\n";
  ///////////////////////////////

  bool verbose = false;
  bool debug = false;
  bool insitu = false;
  bool herm = false;
  std::string which("LM");
  std::string filename;
  int nev = 5;
  int blockSize = 1;
  MT tol = 1.0e-6;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("herm","nonherm",&herm,"Solve Hermitian or non-Hermitian problem.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix (assumes non-Hermitian unless specified otherwise).");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;
  if (filename == "") {
    // get default based on herm
    if (herm) {
      filename = "mhd1280b.cua";
    }
    else {
      filename = "mhd1280a.cua";
    }
  }

  if (MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  // Get the data from the HB file
  int dim,dim2,nnz;
  int rnnzmax;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  if (MyPID == 0) {
    info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
    // find maximum NNZ over all rows
    cout << "Matrix Info: Row Num = " << dim << ", Col Num = " << dim2 << ", nnz = " << nnz << "\n";
    vector<int> rnnz(dim,0);
    for (int *ri=rowind; ri<rowind+nnz; ++ri) {
      ++rnnz[*ri-1];
    }
    rnnzmax = *std::max_element(rnnz.begin(),rnnz.end());
  }
  else {
    // address uninitialized data warnings
    dvals = NULL;
    colptr = NULL;
    rowind = NULL;
  }

  Teuchos::broadcast(*comm,0,&info);
  Teuchos::broadcast(*comm,0,&nnz);
  Teuchos::broadcast(*comm,0,&dim);
  Teuchos::broadcast(*comm,0,&rnnzmax);
  if (info == 0 || nnz < 0) {
    if (MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }


 // create map
  RCP<const Map<> > map = rcp (new Map<> (dim, 0, comm));
  RCP<CrsMatrix<ST> > K = rcp (new CrsMatrix<ST> (map, rnnzmax));
  if (MyPID == 0) {
    // Convert interleaved doubles to complex values
    // HB format is compressed column. CrsMatrix is compressed row.
    const double *dptr = dvals;
    const int *rptr = rowind;
    for (int c=0; c<dim; ++c) {
      for (int colnnz=0; colnnz < colptr[c+1]-colptr[c]; ++colnnz) {
        K->insertGlobalValues (static_cast<GO> (*rptr++ - 1), tuple<GO> (c), tuple (ST (dptr[0], dptr[1])));
        dptr += 2;
      }
    }
  }
  if (MyPID == 0) {
    // Clean up.
    free( dvals );
    free( colptr );
    free( rowind );
  }
  K->fillComplete();
  //cout << "Matrix read by Arnoldi: \n" << *K << endl;



// Create initial vectors
  RCP<MV> ivec = rcp( new MV(map,blockSize) );
  ivec->randomize ();

  // Create eigenproblem
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,ivec) );
  //
  // Inform the eigenproblem that the operator K is symmetric
  problem->setHermitian(herm);
  //
  // Set the number of eigenvalues requested
  problem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  bool boolret = problem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
           << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }


  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  /*
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
  if (verbose) {
    verbosity += Anasazi::IterationDetails;
  }
  */
  if (debug) {
    verbosity += Anasazi::Debug;
  }



  // Eigensolver parameters
  int numBlocks = 8;
  int maxRestarts = 10;
  //
  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "In Situ Restarting", insitu );

//
  // Create the solver manager
  Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> MySolverMgr(problem, MyPL);

  int length = 5;
  int i, end = 0;

  while(!end){

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    success = (returnCode == Anasazi::Converged);

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    RCP<MV> evecs = sol.Evecs;
    int numev = sol.numVecs;
    ST *Evalues  = new ST [numev];

    if (numev > 0) {
      std::ostringstream os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os.precision(6);

      // Compute the direct residual
      std::vector<MT> normV( numev );
      Teuchos::SerialDenseMatrix<int,ST> T (numev, numev);
      for (int i=0; i<numev; i++) {
        T(i,i) = ST(sol.Evals[i].realpart,sol.Evals[i].imagpart);
        Evalues[i] = ST(sol.Evals[i].realpart,sol.Evals[i].imagpart);
      }

      RCP<MV> Kvecs = MVT::Clone( *evecs, numev );

      OPT::Apply( *K, *evecs, *Kvecs );

      MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
      MVT::MvNorm( *Kvecs, normV );

      os << "Direct residual norms computed in BlockKrylovSchurComplex_test.exe" << endl
         << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
         << "----------------------------------------" << endl;
      for (int i=0; i<numev; i++) {
        if ( SCT::magnitude(T(i,i)) != SCT::zero() ) {
          normV[i] = SCT::magnitude(normV[i]/T(i,i));
        }
        os << std::setw(20) << T(i,i) << std::setw(20) << normV[i] << endl;
        success = (normV[i] < tol);
      }
      
      if (MyPID==0) {
        cout << endl << os.str() << endl;
      }
    }

    if (MyPID==0) {
      if (success)
        cout << "End Result: TEST PASSED" << endl;
      else
        cout << "End Result: TEST FAILED" << endl;
    }

    mpi_lsa_com_cplx_array_send(&COMM_FATHER, &numev, Evalues);
    printf("Arnoldi send eigenvalues to FATHER\n");

    //check if any type to receive
    if(!mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type)){
      if(arank == 0){
        printf("Info ]> ERAM Receive signal information from Father\n");
      }
      //exit if receive the exit signal
      if(exit_type == 666){
        if(arank == 0){
          printf("Info ]> ERAM exit\n");
        }
        end = 1;
        break;
      }
    }
  }


  int exit_signal = 777;
  mpi_lsa_com_type_send(&COMM_FATHER, &exit_signal);

  MPI_Comm_free(&COMM_FATHER);

  //return 0;
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

}
