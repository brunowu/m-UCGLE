#include "petsc.h"
//#include "../../Libs/mpi_lsa.h"
//#include "../../Libs/mpi_lsa_com.h"
//#include "../../Libs/data_rw.h"
#include "../../include/Solvers/Utils/libs.hpp"
#include "../Utils/convhulPETSc.h"
#include "../Utils/ellipsePETSc.h"
#include "./precondPETSc.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "../../include/Libs/mpi_lsa_com.hpp"

#ifndef EIGEN_MIN
#define EIGEN_MIN 5
#endif

#ifndef EIGEN_ALL
#define EIGEN_ALL 20	
#endif

#ifndef EIGEN_MAX
#define EIGEN_MAX 20	
#endif


int main(){
	/* variables */
	int vector_size;
	PetscInt end,eigen_received,eigen_total,eigen_max;
	int cumul;
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	int i,type=0;
	PetscInt info;
	PetscBool flag,data_load,data_export,continuous_export,data_load_any;
	PetscScalar * data,*eigen_cumul,*eigen_tri,*d,*c;
	PetscReal a_ell,c_ell,d_ell,d_reel;
	int data_size;
	int chsign;
	PetscInt mu1,mu2,mu,result_array_size;
	PetscReal alpha,*beta,*delta;
	PetscReal scalar_tmp;
	PetscReal *eta;
	PetscInt ls_eigen_min, ls_eigen; // use default values
	PetscErrorCode ierr;
 	PetscScalar * result_array,*data_buffer;/*[EIGEN_ALL*3+2]*/
	Vec * v;
	sprintf(load_path,"./lsqr.bin");
	sprintf(export_path,"./lsqr.bin");
	/* check if there is arguments for ls */
	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_eigen_min",&ls_eigen_min,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen_min=EIGEN_MIN;
	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_eigen",&ls_eigen,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen=EIGEN_ALL;
	/* check the number of eigenvalues that one will receive from arnoldi */
	ierr=PetscOptionsGetInt(NULL,PETSC_NULL,"-ksp_ls_k_param",&eigen_max,&flag);CHKERRQ(ierr);
	if(!flag)eigen_max=ls_eigen;

	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_ls_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_ls_load_any",&data_load_any);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_ls_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(NULL,PETSC_NULL,"-ksp_ls_cexport",&continuous_export);CHKERRQ(ierr);

	ierr=PetscMalloc((vector_size)*sizeof(PetscScalar),&eigen_tri);
	ierr=PetscMalloc((vector_size)*sizeof(PetscScalar),&eigen_cumul);
	ierr=PetscMalloc(((vector_size)+1)*sizeof(PetscScalar),&d);
	ierr=PetscMalloc(((vector_size)+1)*sizeof(PetscScalar),&c);
	ierr=PetscMalloc((vector_size)*sizeof(PetscScalar),&data);
	ierr=PetscMalloc(eigen_max*sizeof(PetscScalar),&data_buffer);

	/* data that will be sended to GMRES for it's preconditionning step */
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&eta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&beta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&delta);
	ierr=PetscMalloc((eigen_max*3+2)*sizeof(PetscScalar),&result_array);
	result_array_size=2+3*eigen_max;

	for(i=0;i<(vector_size);i++){
	  eigen_tri[i]=(PetscScalar)0.0;
	  eigen_cumul[i]=(PetscScalar)0.0;
	}
	for(i=0;i<(vector_size)+1;i++){
	  d[i]=(PetscScalar)0.0;
	  c[i]=(PetscScalar)0.0;
	}
	cumul=0;
	eigen_received=0;
	eigen_total=0;
	end=0;
	eigen_total=0;
	ls_eigen=0;
	while(!end){
		/*in any case clear data array*/
		for(i=0;i<eigen_max;i++)
		  data[i]=(PetscScalar)0.0+PETSC_i*(PetscScalar)0.0;
		/* if received something */
			/* we received data or load it depending on the flags (for first step only*/

				/* first we gonna remove some non-needed values */
				epurer(data,&data_size);
				/* add them to the accumulated eigenvalues */
				/* if full renew full eigenvalues */
				if(eigen_total+data_size>vector_size) eigen_total=0;
				/* select eigenvalues */
				for(i=0;i<data_size;i++){
					eigen_cumul[eigen_total+i]=data[i];
				}
				eigen_total+=data_size;
				if(cumul<eigen_total) cumul=eigen_total;
				for(i=0;i<cumul;i++){
					eigen_tri[i]=eigen_cumul[i];
				}
			
			eigen_received+=data_size;
			/* if we didn't received enough eigenvalues */
			if(eigen_received<ls_eigen_min && data) continue;
			else {
				eigen_received=0;
				tri(eigen_tri,cumul,&chsign);
				mu1=0;
				mu2=0;
				/* convex hull computation */
				if(chsign>0){
					//keepPositif(eigen_tri,&cumul);
					convhull(eigen_tri, c, d, chsign, &mu1, 0, 0);
					printf("@} LSQR convhul negatif chsigne %d cumul %d mu1 %d\n",chsign,cumul,mu1);

				}
				if(chsign<cumul){
					convhull(eigen_tri, c, d, cumul-chsign, &mu2, chsign, mu1);
					printf("@} LSQR convhul positif chsigne %d cumul %d mu1 %d mu2 %d\n",chsign,cumul,mu1,mu2);

				}
				mu=mu1+mu2;
				/* Ellipse computation */
				ierr=ellipse(c,  d, mu+1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info);CHKERRQ(ierr);
				if(fabs(d_ell)<epsilon()) d_ell = 1.;
                                        PetscPrintf(PETSC_COMM_WORLD, "!!!!!!!!fabs(a_ell) = %f\n\n", fabs(a_ell));
				if(fabs(a_ell)<epsilon())
				{
					ls_eigen=0;
//					PetscPrintf(PETSC_COMM_WORLD, "!!!!!!!!fabs(a_ell) = %f\n\n", fabs(a_ell));
				}
				else {
					LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta,
					  delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max);
				}
			}

		if(ls_eigen>1){
			/* place the computed results inside the array */
			scalar_tmp=(PetscReal)ls_eigen;
			ierr=PetscMemcpy(&result_array[0],&scalar_tmp,1*sizeof(PetscReal));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[1],&alpha,1*sizeof(PetscReal));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2],eta,ls_eigen*sizeof(PetscReal));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+ls_eigen],beta,ls_eigen*sizeof(PetscReal));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+(2*ls_eigen)],delta,ls_eigen*sizeof(PetscReal));CHKERRQ(ierr);

			/* and send it */
			result_array_size=2+3*ls_eigen;
                        PetscPrintf(PETSC_COMM_WORLD, "LS has sent the parameters\n");
		} //else {                        PetscPrintf(PETSC_COMM_WORLD, "@@@@@@@@@@@> LS CAN NOT PRECPARE TO SEND TO GMRES\n\n");}
		if(ls_eigen>1){
		  ls_eigen=0;
		}
	}

	/* Free the arrays */
	PetscFree(eigen_tri);
	PetscFree(eigen_cumul);
	PetscFree(d);
	PetscFree(c);
	PetscFree(data);
	PetscFree(eta);
	PetscFree(beta);
	PetscFree(delta);
 	ierr=PetscFree(result_array);CHKERRQ(ierr);
	return 0;
}