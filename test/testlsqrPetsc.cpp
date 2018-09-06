#include "testlsqrPetsc.hpp"

static const char help[] = "LS Polynomial";

int main(int argc, char **argv){
	PetscErrorCode ierr;
	Vec x,vec_rhs, vec_sol;
	Mat Amat;

	Vec x_tmp, r0_tmp, w_1_tmp, w1_tmp, w0_tmp, r1_tmp, sol_tmp, vec_tmp;
	PetscScalar alpha;
  	PetscScalar * eta,*delta,*beta;

	ierr=PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Initializing PETSc/SLEPc\n");
	
	/*Load data*/
	ierr=loadInputs(&Amat,&vec_rhs,&x);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Data loaded\n");

	PetscScalar data_tmp[32];

	data_tmp[0] = 10.000000; data_tmp[1] = -0.633583; data_tmp[2] = -1.578017; data_tmp[3] = -0.105339;
	data_tmp[4] = -0.172723; data_tmp[5] = -0.072221; data_tmp[6] = -0.061025; data_tmp[7] = -0.038023;
	data_tmp[8] = -0.009229; data_tmp[9] = -0.007351; data_tmp[10] = -0.023021; data_tmp[11] = -0.009441;
	data_tmp[12] = 0.134621; data_tmp[13] = 1.112780; data_tmp[14] = 0.371290; data_tmp[15] = 0.575088;
	data_tmp[16] = 0.466703; data_tmp[17] = 0.512561; data_tmp[18] = 0.490791; data_tmp[19] = 0.500618;
	data_tmp[20] = 0.496076; data_tmp[21] = 0.498153; data_tmp[22] = 0.000000; data_tmp[23] = -0.843539;
	data_tmp[24] = -0.102049; data_tmp[25] = -0.305847; data_tmp[26] = -0.197462; data_tmp[27] = -0.243319;
	data_tmp[28] = -0.221550; data_tmp[29] = -0.231377; data_tmp[30] = -0.226835; data_tmp[31] = -0.228912;

	PetscInt size = (PetscInt)PetscRealPart(data_tmp[0]);
	PetscMalloc(sizeof(PetscScalar)*size,&eta);
    PetscMalloc(sizeof(PetscScalar)*size,&beta);
    PetscMalloc(sizeof(PetscScalar)*size,&delta);

	ierr=PetscMemcpy(&alpha,&data_tmp[1],1*sizeof(PetscScalar));CHKERRQ(ierr);
	ierr=PetscMemcpy(eta,&data_tmp[2],(size)*sizeof(PetscScalar));CHKERRQ(ierr);
	ierr=PetscMemcpy(beta,&data_tmp[2+(size)],(size)*sizeof(PetscScalar));CHKERRQ(ierr);
	ierr=PetscMemcpy(delta,&data_tmp[2+2*(size)],(size)*sizeof(PetscScalar));CHKERRQ(ierr);

	ierr=VecDuplicate(vec_rhs,&vec_sol);CHKERRQ(ierr);
	ierr=VecSet(vec_sol,(PetscScalar)0.0);CHKERRQ(ierr);

/*
	ierr=VecSet(vec_sol,(PetscScalar)1.0);CHKERRQ(ierr);

	Vec tt;
	ierr=VecDuplicate(vec_sol,&tt);CHKERRQ(ierr);
	ierr = MatMult(Amat,vec_sol,tt);CHKERRQ(ierr);
	ierr=VecCopy(tt,vec_rhs);CHKERRQ(ierr);
*/

	ierr=VecDuplicate(vec_sol,&x_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&r0_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&w_1_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&w1_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&w0_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&r1_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&sol_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(vec_sol,&vec_tmp);CHKERRQ(ierr);

	ierr=VecSet(x_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(r0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(w0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(w1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(r1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecSet(vec_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
	ierr=VecCopy(vec_sol,sol_tmp);CHKERRQ(ierr);

	PetscInt ls_power = 10, i, j;
	PetscReal norm, normb;

	for(i=0;i<size-1;i++){
		//PetscPrintf(PETSC_COMM_WORLD, "delta[%d] = %g\n", i, delta[i]);
	}

	VecNorm(vec_rhs,NORM_2,&normb);
	printf("nromb = %f\n", normb);

	for(j=0;j<ls_power;j++){
		/* r0 = b-Ax*/
		/* put A*x into VEC_TEMP */

		ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
		/* now put residual (-A*x + f) into vec_vv(0) */
		ierr = VecWAXPY(r0_tmp,-1.0,vec_tmp,vec_rhs);CHKERRQ(ierr);
		/* r0 = w0*/
		ierr=VecCopy(r0_tmp,w0_tmp);CHKERRQ(ierr);
		ierr=VecCopy(w0_tmp,x_tmp);CHKERRQ(ierr);
		ierr=VecScale(x_tmp,eta[0]);CHKERRQ(ierr);
		ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);

		/* depending ton the ls polynom size size */
		for(i=0;i<size-1;i++){

		  /* w1=-alpha*w0 - delta[i]*w_1 ((  y = alpha x + delta y. )) (Vec y,PetscScalar alpha,PetscScalar beta,Vec x)*/
		  ierr=VecCopy(w_1_tmp,w1_tmp);CHKERRQ(ierr);
		  ierr=VecAXPBY(w1_tmp,-alpha,-(PetscScalar)delta[i],w0_tmp);CHKERRQ(ierr);
//		  VecNorm(w1_tmp,NORM_2,&norm);
//		  PetscPrintf(PETSC_COMM_WORLD, "r1_tmp_norm[%d] = %f\n", i, norm/normb);
		  /* w1 = w1 - A*w0 */
		  ierr = MatMult(Amat,w0_tmp,vec_tmp);CHKERRQ(ierr);
		  /* y = alpha x + y.  VecAXPY(Vec y,PetscScalar alpha,Vec x)*/
		  ierr = VecAXPY(w1_tmp,1.0,vec_tmp);CHKERRQ(ierr);
		  /* w1 = w1/beta[i] */
		  ierr=VecScale(w1_tmp,1.0/beta[i]);CHKERRQ(ierr);
		  /* w_1 = w0 */
		  ierr=VecCopy(w0_tmp,w_1_tmp);CHKERRQ(ierr);
		  /* w0 = w1 */
		  ierr=VecCopy(w1_tmp,w0_tmp);CHKERRQ(ierr);
		  /* x = x + (w0 * eta[i] ) */
		  ierr=VecAXPY(x_tmp,eta[i+1],w0_tmp);CHKERRQ(ierr);
		}
		/* x1= x1+x*/
		ierr=VecAXPY(sol_tmp,1.0,x_tmp);CHKERRQ(ierr);
		/* put A*x into VEC_TEMP */
		ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
		/* now put residual (-A*x + f) into vec_vv(0) */
		ierr = VecWAXPY(r1_tmp,-1.0,vec_tmp,vec_rhs);CHKERRQ(ierr);
		/* compute norm and see if it's below epsilon */

		VecNorm(r1_tmp,NORM_2,&norm);
		PetscPrintf(PETSC_COMM_WORLD, "r1_tmp_norm[%d] = %f\n", i, norm/normb);


	}

	/*Clean*/
	ierr = VecDestroy(&vec_rhs);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = MatDestroy(&Amat);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"]> Cleaned structures, finalizing\n");

	/*Finalize PETSc*/
	PetscFinalize(); 

	return 0;
}

PetscErrorCode loadInputs(Mat * A, Vec * b, Vec * x){
	PetscErrorCode ierr;
	PetscInt sizex,sizey;
	char bfile[]="-bfile";
	char xfile[]="-xfile";
	
	//load data files
	ierr=loadMatrix(A);CHKERRQ(ierr);
	ierr=loadVector(bfile,b);CHKERRQ(ierr);
	if(*b==NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector b\n");
		ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
		ierr=generateVectorNorm(sizex,b);CHKERRQ(ierr);
	}
	ierr=loadVector(xfile,x);CHKERRQ(ierr);
	if(*x==NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"]> Creating vector x\n");
		ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
		ierr=generateVectorNorm(sizex,x);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode loadMatrix(Mat * A){
	char file[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flag;
	PetscViewer fd;
	PetscInt sizex,sizey;

	/*check args, if no matrix then no work... matrix file is mandatory*/
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-mfile",file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
	if (!flag) {		
		sprintf(err,"Error : mfile is not properly set -> %s\n",file);
		SETERRQ(PETSC_COMM_WORLD,(PetscErrorCode)83,err);
	}

	/* read matrix file */
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",file);

	ierr=MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
	ierr=MatSetType(*A,MATAIJ);CHKERRQ(ierr);

	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr=MatGetSize(*A,&sizex,&sizey);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",sizex,sizey);

	return 0;
}


PetscErrorCode loadVector(char * type_v,Vec * b){
	char file[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flag;
	PetscViewer fd;
	PetscInt size;

	// check if there is a vec file, vector is not mandatory
	ierr=PetscOptionsGetString(NULL,PETSC_NULL,type_v,file,PETSC_MAX_PATH_LEN-1,&flag);CHKERRQ(ierr);
	if (!flag) {		
		PetscPrintf(PETSC_COMM_WORLD,"Error : %s is not properly set\n",type_v);
		*b = NULL;
	}else{
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",file);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr=VecLoad(*b,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		ierr=VecGetSize(*b,&size);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

	return 0;
}


PetscErrorCode generateVectorRandom(int size, Vec * v){
	PetscErrorCode ierr;

	ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
	ierr=generateVector(size,v);CHKERRQ(ierr);
	ierr=VecSetRandom(*v,PETSC_NULL);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Generated Random Vector of size : %d\n",size);	

	return 0;
}


PetscErrorCode generateVectorNorm(int size, Vec * v){
	PetscScalar scal;
	PetscErrorCode ierr;

	ierr=PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");CHKERRQ(ierr);
	ierr=generateVector(size,v);CHKERRQ(ierr);
	scal=1.0;
	ierr=VecSet(*v,scal);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Generated Norm Vector of size : %d\n",size);	

	return 0;
}


PetscErrorCode generateVector(int size, Vec * v){
	PetscErrorCode ierr;

	ierr=VecCreate(PETSC_COMM_WORLD,v);CHKERRQ(ierr);
	ierr=VecSetSizes(*v,PETSC_DECIDE,size);CHKERRQ(ierr);
	ierr=VecSetFromOptions(*v);CHKERRQ(ierr);
	/*initiate the vector to its norm*/

	return 0;
}