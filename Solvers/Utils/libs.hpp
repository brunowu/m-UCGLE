#ifndef _LIBS_H_
#define _LIBS_H_

/* calcul de l'epsilon machine */
double epsilon(void){
	double one=(PetscReal)1.0, two=(PetscReal)2.0,temp=(PetscReal)1.0;

	while(one+temp>one)
		temp/=two;

	return two*temp;
}

#endif
