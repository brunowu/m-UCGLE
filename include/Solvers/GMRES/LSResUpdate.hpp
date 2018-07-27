#ifndef _LS_RES_UPDATE_H_
#define _LS_RES_UPDATE_H_

#ifndef EIGEN_ALL
#define EIGEN_ALL 10
#endif


#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include <Trilinos_Util_iohb.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_RCP.hpp>

#include "../../Libs/mpi_lsa_com.hpp"
#include <unistd.h>


typedef double Scalar;
typedef Teuchos::ScalarTraits<Scalar>         	SCT;
typedef Tpetra::Map<>::local_ordinal_type 		LO;
typedef Tpetra::Map<>::global_ordinal_type 		GO;
typedef Tpetra::MultiVector<Scalar,LO,GO> 		MV;
typedef Belos::MultiVecTraits<Scalar,MV> MVT;
typedef Tpetra::Operator<Scalar,int>         OP;
typedef Belos::OperatorTraits<Scalar,MV,OP> OPT;
typedef Teuchos::ScalarTraits<Scalar> SCT;
typedef SCT::magnitudeType MT;

using Tpetra::global_size_t;
using Tpetra::Map;
using Tpetra::Import;
using Teuchos::RCP;
using Teuchos::rcp;

static int latency_count = 0;

int LSResUpdate(const Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > &problem, int ls_power, int latency, bool uselsp){
	//simulation of the received data from LS
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  	int grank, gsize;
  	int type = 666;
  	int exit_type = 0;

  	MPI_Comm COMM_FATHER;

  	MPI_Comm_size( MPI_COMM_WORLD, &gsize );
  	MPI_Comm_rank( MPI_COMM_WORLD, &grank );

  	MPI_Comm_get_parent( &COMM_FATHER );

  	bool lspuse = true;

  	double data_tmp[EIGEN_ALL*2*3];
  	double tmp[EIGEN_ALL*2*3];

  	int size_data = EIGEN_ALL;

  	int tmp_size; 

  	int i, j, k;

	if(uselsp){

		printf("Hey LS PPOWER = %d, LS Latency = %d\n", ls_power, latency);
	  	latency_count++;

	  	if(latency_count % latency == 0){
			usleep(100000);
	  	} else{
	  		return 1;
	  	}


		if(!mpi_lsa_com_array_recv(&COMM_FATHER, &size_data, data_tmp)){
			if(grank == 0){
				printf("GMRES ]> GMRES has recived data from LS\n");
				for(i = 0; i < size_data; i++){
	            	printf("GMRES ]>: GMRES rank = %d, data[%d] = %f\n",grank, i, data_tmp[i] );
	            }
			}

			for(i = 0; i < EIGEN_ALL * 3 + 2; i++){
				tmp[i] = data_tmp[i];
			}
			tmp_size = (int)data_tmp[i] * 3 + 2;
			if(((int)data_tmp[0]* 3 + 2) != size_data){
				if(grank == 0){
					printf("GMRES ]> GMRES LS polynomial preconditioned data is not consistent \n");
				}
			}
		} else{
			for(i = 0; i < EIGEN_ALL * 3 + 2; i++){
				data_tmp[i] = tmp[i];
			}
			size_data = tmp_size;
			//return 1;
		}

		Teuchos::RCP<const OP> A = problem->getOperator();
		//std::cout <<*A << std::endl;

		Teuchos::RCP<const MV> vec_rhs = problem->getRHS();
		//std::cout <<*B << std::endl;
		
		const Scalar ONE  = SCT::one();

		int numVectors = vec_rhs->getNumVectors();
		std::cout <<"numVectors = : " << numVectors<< std::endl;

	//	RCP<MV> curX = MVT::Clone(*vec_rhs,numVectors);
	//	curX->putScalar(0.0);

		RCP<MV> curX = problem->getCurrLHSVec();

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

		    MVT::MvNorm( *vec_rhs, normB);

		    std::vector<Scalar>::iterator iter_begin = normV.begin();
		    std::vector<Scalar>::iterator iter_end   = normV.end();
		    std::vector<Scalar>::iterator iter;

		    std::vector<Scalar>::iterator iterB_begin = normB.begin();
		    std::vector<Scalar>::iterator iterB_end   = normB.end();
		    std::vector<Scalar>::iterator iterB;

		    normBB = normB.data();
		    normVV = normV.data();

		    for(k = 0; k < numVectors; k++){
		    	std::cout << "LOOP" << j+1 << ", Rank " << rank << ": "<< normVV[k]/normBB[k] << std::endl;
		    }
		}

		MVT::MvAddMv( ONE, *sol_tmp, -ONE, *curX, *sol_tmp );  
		problem->updateSolution( sol_tmp, true );

		delete [] eta;
		delete [] beta;
		delete [] delta;

	}else{

	  return 1;

	}


	return 0;
}


#endif
