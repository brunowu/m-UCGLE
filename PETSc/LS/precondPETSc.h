#include "petsc.h"
#include "petscpc.h"
#include "petscmat.h"
#include "petscksp.h"
#include <stdlib.h>

#define at(a, i, j) (a)[j][i]

#include "../../include/Solvers/Utils/libs.hpp"
#include <complex>

int LSPrecond(PetscReal a_ell, PetscReal d_ell, PetscReal c_ell, PetscReal * eta, PetscReal * alpha, PetscReal * beta,
              PetscReal * delta, PetscComplex * c, PetscComplex * d, PetscInt * mu, PetscInt * nb_eigen, PetscInt * min_eigen,
              PetscInt * nb_eigen_all){

 //////////////////////////////

  PetscInt i, j, k, nu;

  Mat MM, F,fact;

  KSP ksplsqr, kspchol;

  PC  pcchol, pclsqr;

  Vec soln, rhs;

  /* allocate work array gamma and mm_tmp, init it to 0 */
  std::complex<PetscReal> ** gamma  = new std::complex<PetscReal> * [(*nb_eigen_all)+3];
  PetscScalar * fact_tmp, * res ;

  for(i = 0; i < (*nb_eigen_all)+3; i++){
    gamma[i] = new std::complex<PetscReal> [(*nb_eigen_all)+3];
    for(j = 0; j < (*nb_eigen_all)+3; j++){
      gamma[i][j].real(0.0);
      gamma[i][j].imag(0.0);

    }
  }

  PetscReal ** mm_tmp = new PetscReal * [(*nb_eigen_all)+1];

  for(i = 0; i < (*nb_eigen_all)+1; i++){
    mm_tmp[i] = new PetscReal [(*nb_eigen_all)+1];
    for(j = 0; j < (*nb_eigen_all)+1; j++){
      mm_tmp[i][j] = 0.0;
    }
  }

  for(i = 0; i < *nb_eigen_all; i++){
      beta[i]= 0.0;
      delta[i] = 0.0;
  }

  /* begin computations */
  beta[0] = a_ell / 2.0;
  *nb_eigen = *nb_eigen_all;

  for(i = 0; i < *nb_eigen_all; i++){
    delta[i + 1] = (d_ell) / (4.0 * beta[i]);
    beta[i + 1] = (a_ell) - delta[i + 1];
    if(std::abs(beta[i + 1]) < epsilon()){
      *nb_eigen = i;
      if(*nb_eigen < * min_eigen){
        return 0;
      }
    }
  }

  *alpha = c_ell;

  /* computation of gamma */
  for(nu = 0; nu < *mu; nu++){
//    printf("d[%d] = %f+%fi, c[%d] = %f+%fi\n", nu,d[nu].real(), d[nu].imag(),nu,c[nu].real(), c[nu].imag());
    for(i = 0; i < * nb_eigen_all + 1; i++){
      for(j =0; j < * nb_eigen_all + 1; j++){
        gamma[i][j].real(0.0);
        gamma[i][j].imag(0.0);
      }
    }
      gamma[1][1].real(1.0);
      gamma[1][1].imag(0.0);

    for(j = 1; j <= *nb_eigen_all; j++){
      for(i = 1; i <= j + 1; i++){

        gamma[i][j + 1].real(

         ((d[nu]).real()/2.0*((gamma[i+1][j]).real()+(gamma[i-1][j]).real())
        - (d[nu]).imag()/2.0*((gamma[i+1][j]).imag()+(gamma[i-1][j]).imag())
        + ((c[nu]).real()-*alpha)*(gamma[i][j]).real()
        - (c[nu]).imag()*(gamma[i][j]).imag()
        - delta[j-1]*(gamma[i][j-1]).real())/ beta[j-1]

        );

        gamma[i][j + 1].imag(gamma[i][j + 1].imag());

        gamma[i][j + 1].real(gamma[i][j + 1].real());

        gamma[i][j + 1].imag((d[nu].real()/2.0*(gamma[i+1][j].imag()+gamma[i-1][j].imag())
        + d[nu].imag()/2.0*(gamma[i+1][j].real()+gamma[i-1][j].real())
        + ((c[nu]).real() - *alpha)*gamma[i][j].imag()
        + c[nu].imag()*gamma[i][j].real()
        - delta[j-1]*gamma[i][j-1].imag())/ beta[j-1]);

 //       printf("gamma[%d][%d+1] = %f + %fi\n", i, j+1, gamma[i][j+1].real(), gamma[i][j+1].imag());


      }
      gamma[0][j + 1] = gamma[2][j + 1];
    }

    /* computation of MM */
    for(j = 0; j <= *nb_eigen_all; j++){
      for(i = 0; i <= j; i++){
        mm_tmp[i][j] = mm_tmp[i][j] + 4.*(((gamma[1][j + 1].real())*(gamma[1][i + 1].real())
							    +(gamma[1][j + 1].imag())*(gamma[1][i + 1]).imag()));
        if(i > 1){
  					for(k = 1; k < i; k++){
  						mm_tmp[i][j] = mm_tmp[i][j] + 2.*(((gamma[k + 1][j + 1].real())*(gamma[k + 1][i + 1].real())
  						+(gamma[k + 1][j + 1].imag())*(gamma[k + 1][i + 1].imag())));
  					}
  			 }
      }
    }
  }

	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen_all)+1,(*nb_eigen_all)+1,PETSC_NULL,&MM);

	/* Filling of the lower triangular part */

	for(j=0;j<=*nb_eigen_all;j++){
		for(i=0;i<=j;i++){
			MatSetValue(MM,i,j,(PetscScalar)mm_tmp[i][j],INSERT_VALUES);
			MatSetValue(MM,j,i,(PetscScalar)mm_tmp[i][j],INSERT_VALUES);
		}
	}


	MatAssemblyBegin(MM,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(MM,MAT_FINAL_ASSEMBLY);

/*proceed to factorization*/
	KSPCreate(PETSC_COMM_WORLD,&kspchol);
	KSPSetOperators(kspchol,MM,MM);
	KSPSetType(kspchol,KSPPREONLY);
	KSPGetPC(kspchol,&pcchol);
	PCSetType(pcchol,PCCHOLESKY);
	KSPSetUp(kspchol);

/*get factor matrix*/
	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen_all)+1,(*nb_eigen_all)+1,PETSC_NULL,&fact);
	MatAssemblyBegin(fact,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(fact,MAT_FINAL_ASSEMBLY);

	PCFactorGetMatrix(pcchol,&fact); 

	/*Create the matrix operator that will be used in the QR factorization*/
	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen)+1,(*nb_eigen)+1,PETSC_NULL,&F);

	/*get matrix array, fact is of dense format so wouldn't be a problem for addressing*/
	MatDenseGetArray(fact,&fact_tmp);
	MatSetValue(F,0,0,(PetscScalar)(*alpha)*(PetscScalar)fact_tmp[0]+(PetscScalar)beta[0]*(PetscScalar)fact_tmp[1],INSERT_VALUES);
	MatSetValue(F,1,0,(PetscScalar)beta[0]*(PetscScalar)fact_tmp[1+((*nb_eigen_all)+1)],INSERT_VALUES);

	for(j=1;j<*nb_eigen;j++){
	  for(i=0;i<j;i++){
	    MatSetValue(F,i,j,(PetscScalar)delta[j-1]*(PetscScalar)fact_tmp[(j-1)+i*((*nb_eigen_all)+1)]
			      +(PetscScalar)(*alpha)*(PetscScalar)fact_tmp[j+i*((*nb_eigen_all)+1)]
			      +(PetscScalar)(beta[j])*(PetscScalar)fact_tmp[(j+1)+i*((*nb_eigen_all)+1)],INSERT_VALUES);
	  }

	  MatSetValue(F,j,j,(PetscScalar)(*alpha)*(PetscScalar)fact_tmp[j+j*((*nb_eigen_all)+1)]
			    +(PetscScalar)beta[j]*(PetscScalar)fact_tmp[(j+1)+j*((*nb_eigen_all)+1)],INSERT_VALUES);
	  if(j+1<=*nb_eigen){
	    MatSetValue(F,j+1,j,(PetscScalar)beta[j]*(PetscScalar)fact_tmp[(j+1)+(j+1)*((*nb_eigen_all)+1)],INSERT_VALUES);
	  }
	}


/* set the vectors*/
	MatCreateVecs(F,&soln,&rhs);

/*set the solution to zero*/
	VecSet(soln,(PetscScalar)0.0);
 	VecSet(rhs,(PetscScalar)0.0);

/* rhs[0] must be setted to beta*/
	VecSetValue(rhs,0,fact_tmp[0],INSERT_VALUES);

	VecAssemblyBegin(soln);
	VecAssemblyEnd(soln);
	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

/*no longer need to access factored matrix, restore it*/
	MatDenseRestoreArray(fact,&fact_tmp);


	/*assemble F for processing*/
	MatAssemblyBegin(F,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(F,MAT_FINAL_ASSEMBLY);


	/*create the lsqr solver context and set it up*/
	KSPCreate(PETSC_COMM_WORLD,&ksplsqr);
	KSPSetOperators(ksplsqr,F,F);
	KSPGetPC(ksplsqr,&pclsqr);
	PCSetType(pclsqr,PCNONE);
	KSPSetType(ksplsqr,KSPLSQR);
	KSPSetInitialGuessNonzero(ksplsqr,PETSC_TRUE);
	KSPSetUp(ksplsqr);

	KSPSolve(ksplsqr, rhs, soln);

	PetscInt its;
	KSPGetIterationNumber(ksplsqr,&its);

/* extract solution elements and place them into eta*/
	VecGetArray(soln,&res);

	 for(i=0;i<*nb_eigen_all;i++){
	  eta[i]=PetscRealPart(res[i]);
	}

	VecRestoreArray(soln,&res);


	KSPDestroy(&kspchol);
	KSPDestroy(&ksplsqr);
	VecDestroy(&rhs);
	VecDestroy(&soln);
	MatDestroy(&F);
	MatDestroy(&MM);

  //free the allocate memory
  for(i = 0; i < (*nb_eigen_all) + 3; i++){
    delete [] gamma[i];
  }

  delete [] gamma;

  for(i = 0; i < (*nb_eigen_all) + 1; i++){
    delete [] mm_tmp[i];
  }

  delete [] mm_tmp;


  return 0;
}
