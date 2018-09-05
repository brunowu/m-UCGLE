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

typedef std::complex<double>                  ST;
typedef Teuchos::ScalarTraits<ST>             SCT;
typedef SCT::magnitudeType                    MT;

typedef Tpetra::Map<>::local_ordinal_type     LO;
typedef Tpetra::Map<>::global_ordinal_type    GO;
typedef Tpetra::Operator<ST,int>              OP;
typedef Tpetra::CrsMatrix<ST,LO,GO>           MAT;
typedef Tpetra::MultiVector<ST,LO,GO>         MV;
typedef Belos::MultiVecTraits<ST,MV>          MVT;
typedef Belos::OperatorTraits<ST,MV,OP>       OPT;

int main(int argc, char *argv[]){
	
	Teuchos::GlobalMPISession mpisess (&argc, &argv, &std::cout);
	Teuchos::oblackholestream blackhole;

	using Tpetra::global_size_t;
	using Tpetra::Map;
	using Tpetra::Import;
	using Teuchos::RCP;
	using Teuchos::rcp;

	Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
	int myRank = comm->getRank();

	printf("my rank = %d\n", myRank);

	bool allprint = true;

	int numVectors = 1;

	std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
	RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

	std::string filename("utm300_cp.mtx");

	RCP<MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(filename, comm);


    if (myRank == 0){
      std::cout << "GMRES ]> Matrix Loaded ... " << std::endl;
    }

    //A->describe(*fos, Teuchos::VERB_EXTREME);
    const ST ONE  = SCT::one();

    // get the maps
	RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
	RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

	GO nrows = dmnmap->getGlobalNumElements();

	RCP<Map<LO,GO> > root_map
  	= rcp( new Map<LO,GO>(nrows,myRank == 0 ? nrows : 0,0,comm) );
	RCP<MV> Xhat = rcp( new MV(root_map,numVectors) );
	RCP<Import<LO,GO> > importer = rcp( new Import<LO,GO>(dmnmap,root_map) );


	RCP<MV> vec_rhs = rcp(new MV(dmnmap,numVectors));
	vec_rhs->putScalar(1.0);
	//vec_rhs->describe(*fos, Teuchos::VERB_EXTREME);

	RCP<MV> curX = rcp(new MV(dmnmap,numVectors));
	curX->putScalar(0.0);

	double *data_tmp = new double [32];
	
	data_tmp[0] = 10.000000; data_tmp[1] = -0.633583; data_tmp[2] = -1.578017; data_tmp[3] = -0.105339;
	data_tmp[4] = -0.172723; data_tmp[5] = -0.072221; data_tmp[6] = -0.061025; data_tmp[7] = -0.038023;
	data_tmp[8] = -0.009229; data_tmp[9] = -0.007351; data_tmp[10] = -0.023021; data_tmp[11] = -0.009441;
	data_tmp[12] = 0.134621; data_tmp[13] = 1.112780; data_tmp[14] = 0.371290; data_tmp[15] = 0.575088;
	data_tmp[16] = 0.466703; data_tmp[17] = 0.512561; data_tmp[18] = 0.490791; data_tmp[19] = 0.500618;
	data_tmp[20] = 0.496076; data_tmp[21] = 0.498153; data_tmp[22] = 0.000000; data_tmp[23] = -0.843539;
	data_tmp[24] = -0.102049; data_tmp[25] = -0.305847; data_tmp[26] = -0.197462; data_tmp[27] = -0.243319;
	data_tmp[28] = -0.221550; data_tmp[29] = -0.231377; data_tmp[30] = -0.226835; data_tmp[31] = -0.228912;

	RCP<MV> x_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> r0_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> w_1_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> w1_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> w0_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> r1_tmp = MVT::Clone(*curX,numVectors);
	RCP<MV> vec_tmp = MVT::Clone(*curX,numVectors);

	RCP<MV> sol_tmp = MVT::CloneCopy(*curX);

	x_tmp->putScalar(0.0);
	r0_tmp->putScalar(0.0);
	w_1_tmp->putScalar(0.0);
	w1_tmp->putScalar(0.0);
	w0_tmp->putScalar(0.0);
	r1_tmp->putScalar(0.0);
	vec_tmp->putScalar(0.0);

	std::vector<MT> normV( numVectors );
	std::vector<MT> normB( numVectors );

	double *normBB = new double [numVectors];
	double *normVV = new double [numVectors];

	int size = (int) data_tmp[0];

	double alpha;
	double *eta = new double [size];
	double *beta = new double [size];
	double *delta = new double [size];

	std::memcpy(&alpha,&data_tmp[1], 1*sizeof(double));
	std::memcpy(eta,&data_tmp[2], size*sizeof(double));
	std::memcpy(beta,&data_tmp[2 + size], size*sizeof(double));
	std::memcpy(delta,&data_tmp[2 + 2 * size], size*sizeof(double));

	int ls_power = 1, i, j ,k;

	MVT::MvNorm( *vec_rhs, normB);

	normBB = normB.data();

	if(myRank == 0){printf("nromb = %f\n\n", normBB[0]);}
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
	      	MVT::MvAddMv( -delta[i], *w1_tmp, -alpha, *w0_tmp, *w1_tmp );
	      	MVT::MvNorm( *w1_tmp, normV);

			normVV = normV.data();

			if(myRank == 0){
				for(k = 0; k < numVectors; k++){
					printf("r1_tmp_norm[%d] = %f\n", i, normVV[k]/normBB[k]);
				}
			}

	      	/* w1 = w1 - A*w0 */
	      	A->apply(*w0_tmp,*vec_tmp);
	      	/* y = alpha x + y.*/
	      	MVT::MvAddMv(ONE, *w1_tmp, ONE, *vec_tmp, *w1_tmp );
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

	}





	return 0;
}
