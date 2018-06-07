#ifndef _ARGS_PARSER_H_
#define _ARGS_PARSER_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

int argsParser(int argc, char ** argv, int *gmres_nb, int *arnoldi_nb, int *gmres_proc, int *arnoldi_proc);
char **argsParserGmresExec(int argc, char ** argv, int gmres_nb);
char **argsParserArnoldiExec(int argc, char ** argv, int arnoldi_nb);
char *argsParserLsqrExec(int argc, char ** argv);

#endif
