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

// I/O for Matrix-Market files
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

int main(int argc, char *argv[]) {
  //Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  Teuchos::GlobalMPISession mpisess (&argc, &argv, &std::cout);

  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>         SCT;
  typedef SCT::magnitudeType                  MT;

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
  typedef Belos::MultiVecTraits<Scalar,MV>    MVT;

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
  bool verbose = (myRank==0);
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  const size_t numVectors = 1;

  RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename, comm);
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

  // Create random X
  RCP<MV> X = rcp(new MV(dmnmap,numVectors));
  X->randomize();
  RCP<MV> vec_sol = MVT::Clone(*X,numVectors);
  vec_sol->putScalar(0.0);

  RCP<MV> x_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> r0_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> w_1_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> w1_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> w0_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> r1_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> sol_tmp = MVT::Clone(*vec_sol,numVectors);
  RCP<MV> vec_tmp = MVT::Clone(*vec_sol,numVectors);

  x_tmp->putScalar(0.0);
  r0_tmp->putScalar(0.0);
  w_1_tmp->putScalar(0.0);
  w1_tmp->putScalar(0.0);
  w0_tmp->putScalar(0.0);
  r1_tmp->putScalar(0.0);
  sol_tmp->putScalar(0.0); //later should be replaced by the current solution in GMRES
  vec_tmp->putScalar(0.0);

  RCP<MV> vec_rhs = rcp(new MV(rngmap,numVectors));
  vec_rhs->putScalar(1); //should be replaced by getRHS

  int ls_power = 5;

  int i, j;

  std::vector<MT> normV( numVectors );

  //simulation of the received data from LS
  double data_tmp [32] = {10.000000,-0.720766,-1.395545,-0.660261,-0.527755,-0.457798,-0.368613,-0.319541,-0.259549,-0.189444,-0.168309
                        -0.119417,0.326302,0.576631,0.609613,0.611939,0.612093,0.612104,0.612104,0.612104,0.612104,0.612104,0.000000,
                       0.075973,0.042991,0.040666,0.040511,0.040501,0.040500,0.040500,0.040500,0.040500};



  int size = (int) data_tmp[0];
  double alpha;
  double *eta = new double [size];
  double *beta = new double [size];
  double *delta = new double [size];

  std::memcpy(&alpha,&data_tmp[1], 1*sizeof(double));
  std::memcpy(eta,&data_tmp[2], size*sizeof(double));
  std::memcpy(beta,&data_tmp[2 + size], size*sizeof(double));
  std::memcpy(delta,&data_tmp[2 + 2 * size], size*sizeof(double));


  //A->apply(*x_tmp,*r0_tmp);
  //r0_tmp->describe(*fos, Teuchos::VERB_EXTREME);

  for(j = 0; j < ls_power; j++){
    /* r0 = b-Ax*/
    /* put A*x into vec_tmp */
    A->apply(*sol_tmp,*vec_tmp);
    MVT::MvAddMv( -ONE, *vec_tmp, ONE, *vec_rhs, *r0_tmp );
    w0_tmp = MVT::CloneCopy(*r0_tmp);
    x_tmp = MVT::CloneCopy(*w0_tmp);
    MVT::MvScale(*x_tmp, eta[0]);
    //vec_tmp->describe(*fos, Teuchos::VERB_EXTREME);
    w_1_tmp->putScalar(0.0);
    /* depending on the ls polynom size (int)data_tmp[0] = size */
    for(i = 0; i < size - 1; i++){
      /* w1=-alpha*w0 - delta[i]*w_1 ((  y = alpha x + delta y. ))*/
      w1_tmp = MVT::CloneCopy(*w_1_tmp);
      /*AXPBY*/
      MVT::MvAddMv( -alpha, *w1_tmp, delta[i], *w0_tmp, *w1_tmp );
      /* w1 = w1 - A*w0 */
      A->apply(*w0_tmp,*vec_tmp);
      /* y = alpha x + y.*/
      MVT::MvAddMv( ONE, *w1_tmp, -ONE, *vec_tmp, *w1_tmp );
      /* w1 = w1/beta[i] */
      MVT::MvScale(*w1_tmp, 1/beta[i]);
      /* w_1 = w0 */
      w_1_tmp = MVT::CloneCopy(*w0_tmp);
      /* w0 = w1*/
      w0_tmp = MVT::CloneCopy(*w1_tmp);
      /* x = x + (w0 * eta[i] ) */
      MVT::MvAddMv( ONE, *x_tmp, eta[i + 1], *w0_tmp, *x_tmp );      
    }
    /* update solution, x1= x1+x*/
    MVT::MvAddMv( ONE, *sol_tmp, ONE, *w0_tmp, *sol_tmp );  
    /* put A*x into VEC_TEMP */
    A->apply(*sol_tmp,*vec_tmp);   
    /* now put residual (-A*x + f) into vec_vv(0) */
    MVT::MvAddMv( -ONE, *vec_tmp, ONE, *vec_rhs, *r1_tmp );  
    /* compute norm and see if it's below epsilon */
    MVT::MvNorm( *r1_tmp, normV);

    std::vector<Scalar>::iterator iter_begin = normV.begin();
    std::vector<Scalar>::iterator iter_end   = normV.end();
    std::vector<Scalar>::iterator iter;

    for(iter = iter_begin; iter != iter_end; ++iter){
      std::cout << "LOOP" << j << ": "<< *iter << std::endl;
    }

  }

//  Teuchos::RCP<MV> curX = problem_->getCurrLHSVec();
//  MVT::MvAddMv( 1.0, *curX, 1.0, *update, *curX );

  delete [] eta;
  delete [] beta;
  delete [] delta;

//  Teuchos::TimeMonitor::summarize();

  // LS is done.
  return 0;
}