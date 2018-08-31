#ifndef TRIEP_HPP_

#include "petsc.h"
#include <math.h>
#include <stdlib.h>

PetscErrorCode tri(PetscScalar * vp, PetscInt nValues, PetscInt * ch_signe);

PetscErrorCode epurer(PetscScalar * vp, PetscInt * nValues);

void swap(PetscScalar *a, PetscScalar *b);

void sort(PetscScalar arr[], int beg, int end);

#endif
