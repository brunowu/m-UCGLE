#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <complex>
#include <unistd.h>


//Block KrylovSchur METHODs to approximate the eigevalues

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziTpetraAdapter.hpp" //Anasazi interface to Tpetra
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

// I/O for Matrix-Market files
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

using Tpetra::CrsMatrix;
using Tpetra::Map;
using std::vector;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;
using std::cout;
using std::endl;
using Tpetra::global_size_t;
using Tpetra::Import;


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

  typedef double                              ST;
  typedef std::complex<double>                CST;
  typedef Tpetra::Map<>::global_ordinal_type    GO;
  typedef Tpetra::Map<>::local_ordinal_type     LO;
  typedef Tpetra::MultiVector<ST,LO,GO>     MV;
  typedef Tpetra::MultiVector<CST,LO,GO>     CMV;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  typedef Tpetra::CrsMatrix<ST,LO,GO>     MAT;



  bool success = false;

  const ST ONE = SCT::one ();

  int info = 0;

  //
  // Get the default communicator
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int MyPID = comm->getRank ();

  /*stantard MPI functionality*/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "Arnoldi in Trilinos said: my rank = " << rank << "\n";
  ///////////////////////////////

  Teuchos::oblackholestream blackhole;

  bool verbose = false;
  bool debug = false;
  bool insitu = false;
  bool herm = false;
  std::string which("LM");
  std::string filename;
  int nev = 10;
  int blockSize = 1;
  MT tol = 1.0e-1;
  bool printMatrix   = false;
  bool allprint      = false;
  int numBlocks = 100;
  int maxRestarts = 100;

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
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  cmdp.setOption("numBlocks",&numBlocks,"Number of blocks in Krylov basis.");
  cmdp.setOption("maxRestarts",&maxRestarts,"Number of restarts allowed.");

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
      filename = "utm300.mtx";
    }
  }else{
    printf("Parser the MATRIX file NAME = %s ! \n", filename.c_str());
  }

  std::ostream& out = ( (allprint || (MyPID == 0)) ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  if (MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }


  if (MyPID == 0){
    std::cout << "Anasazi Scalar type == " << Teuchos::typeName(ONE) << std::endl;
  }
  


  if (MyPID == 0){
    std::cout << "ERAM ]> Building problem... " << Teuchos::typeName(ONE) << std::endl;
  }
  
  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename, comm);

  if (MyPID == 0){
    std::cout << "ERAM ]> Matrix Loaded ... " << Teuchos::typeName(ONE) << std::endl;
  }

  //cout << "Matrix read by Arnoldi: \n" << *K << endl;
  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  } else if( verbose ){
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  //A->describe(*fos, Teuchos::VERB_EXTREME);

  // get the maps
  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  GO nrows = dmnmap->getGlobalNumElements();
  RCP<Map<LO,GO> > root_map
    = rcp( new Map<LO,GO>(nrows,MyPID == 0 ? nrows : 0,0,comm) );
  RCP<MV> Xhat = rcp( new MV(root_map,blockSize) );
  RCP<Import<LO,GO> > importer = rcp( new Import<LO,GO>(dmnmap,root_map) );


// Create initial vectors
  RCP<MV> ivec = rcp( new MV(dmnmap,blockSize) );
  ivec->randomize ();

  // Create eigenproblem
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(A,ivec) );
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
  //int verbosity = Anasazi::Errors + Anasazi::Warnings;
  
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
  if (verbose) {
    verbosity += Anasazi::IterationDetails;
  }
  
  if (debug) {
    verbosity += Anasazi::Debug;
  }



  // Eigensolver parameters

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
  int numev;

  while(!end){

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    success = (returnCode == Anasazi::Converged);

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    numev = sol.numVecs;
    printf("NUMVECS = %d\n", numev);

    CST *Evalues  = new CST [numev];

    if (numev > 0) {
      std::ostringstream os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os.precision(6);

      // Compute the direct residual
      std::vector<MT> normV( numev );
      for (int i=0; i<numev; i++) {
        Evalues[i] = CST(sol.Evals[i].realpart,sol.Evals[i].imagpart);
        printf("Evals[%d] = %f + %fi\n", i+1,sol.Evals[i].realpart, sol.Evals[i].imagpart);
      }
   
      if (MyPID==0) {
        cout << endl << os.str() << endl;
      }
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
