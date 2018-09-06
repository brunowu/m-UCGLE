#include "petsc.h"
#include "slepceps.h"
#include <unistd.h>
#include "../../include/Libs/mpi_lsa_com.hpp"
#include <complex>

#ifndef EIGEN_ALL
#define EIGEN_ALL 10
#endif

PetscErrorCode loadInputs(Mat * A, Vec * x);
PetscErrorCode loadMatrix(Mat * A);
PetscErrorCode loadVector(char * type_v,Vec * b);
PetscErrorCode generateVector(PetscInt size, Vec * v);
PetscErrorCode generateVectorRandom(PetscInt size, Vec * v);


