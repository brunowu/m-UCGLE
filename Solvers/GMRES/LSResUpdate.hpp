#ifndef _LS_RES_UPDATE_H_
#define _LS_RES_UPDATE_H_

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

int LSResUpdate(const Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > &problem){
	//simulation of the received data from LS
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  	double data_tmp [32] = {10.000000,-0.720766,-1.395545,-0.660261,-0.527755,-0.457798,-0.368613,-0.319541,-0.259549,-0.189444,-0.168309
                        -0.119417,0.326302,0.576631,0.609613,0.611939,0.612093,0.612104,0.612104,0.612104,0.612104,0.612104,0.000000,
                        0.075973,0.042991,0.040666,0.040511,0.040501,0.040500,0.040500,0.040500,0.040500};

	Teuchos::RCP<const OP> A = problem->getOperator();
	//std::cout <<*A << std::endl;

	Teuchos::RCP<const MV> vec_rhs = problem->getRHS();
	//std::cout <<*B << std::endl;
	
	const Scalar ONE  = SCT::one();

	int numVectors = vec_rhs->getNumVectors();
	std::cout <<"numVectors = : " << numVectors<< std::endl;

	RCP<MV> curX = MVT::Clone(*vec_rhs,numVectors); //RCP<MV> curX = problem->getCurrLHSVec();
	curX->putScalar(0.0);

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
	
	int ls_power = 2;

  	int i, j;

  	std::vector<MT> normV( numVectors );

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

	    std::vector<Scalar>::iterator iter_begin = normV.begin();
	    std::vector<Scalar>::iterator iter_end   = normV.end();
	    std::vector<Scalar>::iterator iter;

	    for(iter = iter_begin; iter != iter_end; ++iter){
	      std::cout << "LOOP" << j+1 << ", Rank " << rank << ": "<< *iter << std::endl;
	    }

	  }

/*
	if (flag && ls){
		perform ls part;
		MVT::MvAddMv( 0.0, *newX, 1.0, *newX, *curX ); //update with the LS polynomial solution
	}else{
		Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
		problem->updateSolution( update, true );
		return 1;
	}
*/
	delete [] eta;
	delete [] beta;
	delete [] delta;

	return 0;
}


#endif