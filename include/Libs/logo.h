#ifndef __LOGO_H__
#define __LOGO_H__

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void border_print(void)
{
	printf(
	"==================================================================="
	"=============\n");
}

void border_print2(void)
{
	printf(
	"-------------------------------------------------------------------"
	"-------------\n");
}

void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}


void logo()
{
	border_print();

	printf("\n");
	printf("\n");
	printf("\n");
	center_print("MUCGLE: Multiple Unite and Conquer GMRES/LS-ERAM Method", 79);

	printf("\n");
	printf("\n");
	printf("\n");
	border_print();
    center_print("Developed by Xinzhe WU at Maison de la Simulation, France", 79);
	border_print();

	printf("\n\n\n");

}

#endif
