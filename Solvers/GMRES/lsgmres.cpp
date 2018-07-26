#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresLsSolMgr.hpp"

// I/O for Matrix-Market files
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

#include "LSResUpdate.hpp"

int main(int argc, char *argv[]){
//Tpetra::ScopeGuard tpetraScope(&argc,&argv);
	Teuchos::GlobalMPISession mpisess (&argc, &argv, &std::cout);

  	typedef double 									Scalar;
  	typedef Teuchos::ScalarTraits<Scalar>         	SCT;
  	typedef SCT::magnitudeType                  	MT;

  	typedef Tpetra::Map<>::local_ordinal_type 		LO;
  	typedef Tpetra::Map<>::global_ordinal_type 		GO;
	typedef Tpetra::Operator<Scalar,int>         	OP;
  	typedef Tpetra::CrsMatrix<Scalar,LO,GO> 		MAT;
  	typedef Tpetra::MultiVector<Scalar,LO,GO> 		MV;
  	typedef Belos::MultiVecTraits<Scalar,MV>    	MVT;
  	typedef Belos::OperatorTraits<Scalar,MV,OP> 	OPT;

  	const Scalar ONE  = SCT::one();

  	using Tpetra::global_size_t;
  	using Tpetra::Map;
  	using Tpetra::Import;
  	using Teuchos::RCP;
  	using Teuchos::rcp;

  	//
  	// Get the default communicator
  	//
  	Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  	int myRank = comm->getRank();

  	Teuchos::oblackholestream blackhole;

  	bool printMatrix   = false;
  	bool allprint      = false;
  	bool verbose 	   = (myRank==0);
  	bool debug 		   = false;
  	std::string filename("utm300.mtx");
  	int frequency 	   = -1;
  	int numVectors 	   = 2;
	int blocksize 	   = 100;
  	int numblocks 	   = 5;
  	double tol 		   = 1.0e-5;
  	bool precond 	   = false;
  	bool dumpdata 	   = false;
  	bool reduce_tol;
  	Teuchos::ParameterList mptestpl;
  	Belos::ReturnType ret;
  	double t1, t2, t3, t4, t5, t6;

  	Teuchos::CommandLineProcessor cmdp(false,true);
  	cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  	cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  	cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
 	cmdp.setOption("tol",&tol,"Relative residual tolerance used by solver.");
 	cmdp.setOption("num-rhs",&numVectors,"Number of right-hand sides to be solved for.");
 	cmdp.setOption("block-size",&blocksize,"Block size to be used by the solver.");
 	cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
	cmdp.setOption("num-blocks",&numblocks,"Number of blocks in the Krylov basis.");
	cmdp.setOption("reduce-tol","fixed-tol",&reduce_tol,"Require increased accuracy from higher precision scalar types.");
  	cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  	cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  	cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  	cmdp.setOption("dump-data","no-dump-data",&dumpdata,"Dump raw data to data.dat.");

  	if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    	return -1;
  	}

  	if (debug) {
    	verbose = true;
  	}

  	if (!verbose) {
    	frequency = -1;
  	}

	mptestpl.set( "Block Size", blocksize );
 	mptestpl.set( "Num Blocks", numblocks);

  	mptestpl.set<MT>( "Convergence Tolerance", tol );
  	mptestpl.set( "Timer Label",  Teuchos::typeName(ONE) );


	int verbLevel = Belos::Errors + Belos::Warnings;
/*
	if (debug) {
    	verbLevel |= Debug;
  	}
*/
  	if (verbose) {
    	verbLevel |= Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
  	}
  	mptestpl.set( "Verbosity", verbLevel );
  	if (verbose) {
    	if (frequency > 0) {
      		mptestpl.set( "Output Frequency", frequency );
    	}
  	}

  	std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
  	RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

	if (myRank==0){
		std::cout << "Testing Scalar == " << Teuchos::typeName(ONE) << std::endl;
	}

  	/*build the problem by loading matrix market from local file systems*/
  	if (myRank==0){
  		std::cout << "GMRES ]> Building problem..." << std::endl;
  	}
  	
  	
	//	Teuchos::TimeMonitor localtimer(btimer); 
	t1 = MPI_Wtime();

  	RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename, comm);
  	if(myRank == 0){
  		std::cout << "GMRES ]> Load matrix: " << filename << " from local filesystem" << std::endl;
  	}

  	if( printMatrix ){
    	A->describe(*fos, Teuchos::VERB_EXTREME);
  	}

  	else if( verbose ){
    	std::cout << std::endl << A->description() << std::endl << std::endl;
  	}

  	// get the maps
  	RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  	RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  	GO nrows = dmnmap->getGlobalNumElements();
  	RCP<Map<LO,GO> > root_map
    	= rcp( new Map<LO,GO>(nrows,myRank == 0 ? nrows : 0,0,comm) );
  	RCP<MV> Xhat = rcp( new MV(root_map,numVectors) );
  	RCP<Import<LO,GO> > importer = rcp( new Import<LO,GO>(dmnmap,root_map) );

  	// generate randomly the final solution of systems as X
  	RCP<MV> X = rcp(new MV(dmnmap,numVectors));
  	//X->randomize();
  	MVT::MvRandom(*X);
  	//X->describe(*fos, Teuchos::VERB_EXTREME);

  	/*generate the right hand sides B by given X and A */
  	RCP<MV> B = rcp(new MV(dmnmap,numVectors));
  	OPT::Apply( *A, *X, *B );
  	MVT::MvInit( *X, 0.0 );

	RCP<Belos::LinearProblem<Scalar,MV,OP> > problem = rcp( new Belos::LinearProblem<Scalar,MV,OP>(A,X,B) );

	problem->setLabel(Teuchos::typeName(SCT::one()));

	TEUCHOS_TEST_FOR_EXCEPT(problem->setProblem() == false);

//  LSResUpdate(problem);


  t2 = MPI_Wtime();

	if(myRank == 0) printf("Building TIME = %f\n", t2 - t1 ); 
 	
	/////////////////////////////////////////////////
	/*construct GMRES solver*/

	RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
	if (myRank==0){
  		std::cout << "GMRES ]> Construct solver ..." << std::endl;
  	}
  	
  	t3 = MPI_Wtime();
  	//Teuchos::TimeMonitor localtimer(ctimer); 
  	solver = rcp(new Belos::BlockGmresSolMgr<Scalar,MV,OP>( problem, rcp(&mptestpl,false) ));
  	t4 = MPI_Wtime();
  	
  	if(myRank == 0) printf("Construct TIME = %f\n", t4 - t3 ); 

  	if(myRank == 0) printf("GMRES ]> Solving problem ...\n" ); 

	int diter;

	t5 = MPI_Wtime();
	try {
      ret = solver->solve();
    }

    catch (std::exception &e) {
      std::cout << "Caught exception: " << std::endl << e.what() << std::endl;
      ret = Belos::Unconverged;
    }

    t6 = MPI_Wtime();

	if(myRank == 0) printf("Solving TIME = %f\n", t6 - t5 ); 

	diter = solver->getNumIters();

	if (ret == Belos::Converged){
		if(myRank == 0){
			printf("GMRES ]> This GMRES Component is converged\n");
		}
	}
	else{
		if(myRank == 0){
			printf("GMRES ]> This GMRES Component cannot be converged with given parameters\n");
		}
	}


	return 0;


}