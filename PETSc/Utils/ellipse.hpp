#ifndef ELLIPSE_HPP_

#include "petsc.h"
#include <math.h>
#include <stdlib.h>

PetscErrorCode ellipse(PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt mu, PetscReal * co, PetscReal * ao2, PetscReal * do2, PetscReal * dr, PetscInt * info);

#endif