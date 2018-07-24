#ifndef _GMRES_LSA_PRECOND_H_
#define _GMRES_LSA_PRECOND_H_

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

typedef MultiVecTraits<ScalarType,MV> MVT;
typedef OperatorTraits<ScalarType,MV,OP> OPT;
typedef Teuchos::ScalarTraits<ScalarType> SCT;
typedef typename SCT::magnitudeType MagnitudeType;

int GmresLSAPrecond(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem){
	RCP<const OP> Amat = problem->getOperator();

	if (flag && ls){
		perform ls part;
		MVT::MvAddMv( 0.0, *newX, 1.0, *newX, *curX ); //update with the LS polynomial solution
	}else{
		Teuchos::RCP<MV> update = block_gmres_iter->getCurrentUpdate();
		problem->updateSolution( update, true );
		return 1;
	}

}


#endif