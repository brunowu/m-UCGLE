#ifndef _ARGS_PARSER_H_
#define _ARGS_PARSER_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

// -nb_gmres ${GNB} -proc_gmres ${GPNB}
/*function to parser the args to setup the method*/
int argsParser(int argc, char ** argv, int *gmres_nb, int *arnoldi_nb, int *gmres_proc, int *arnoldi_proc){

  int i, j;

  int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;

  for(i = 0; i < argc; i++){
    //nb gmres executable
    if(strcmp(argv[i],"-nb_gmres") == 0){
      cnt1 = 1;
        if(i + 1 < argc){
          i++;
          *gmres_nb = atoi(argv[i]);
        }
      }
      //nb each gmres proc
      else if(strcmp(argv[i],"-proc_gmres") == 0){
        cnt2 = 1;
        if(i + 1 < argc){
          i++;
          *gmres_proc = atoi(argv[i]);
        }
      }
      //nb arnoldi executable
      else if(strcmp(argv[i],"-nb_arnoldi") == 0){
        cnt3 = 1;
        if(i + 1 < argc){
          i++;
          *arnoldi_nb = atoi(argv[i]);
        }
      }
      //nb of proc for each arnoldi
      else if (strcmp(argv[i],"-proc_arnoldi") == 0){
        cnt4 = 1;
        if(i + 1 < argc){
          i++;
          *arnoldi_proc = atoi(argv[i]);
        }
      }
    }

    //if well parsered, show the initailization information
    if(cnt1 && cnt2 && cnt3 && cnt4){
      printf("Initialization ]> GMRES executable number is %d \n", *gmres_nb);
      printf("Initialization ]> Each GMRES Proc number is %d \n", *gmres_proc);
      printf("Initialization ]> ERAM executable number is %d \n", *arnoldi_nb);
      printf("Initialization ]> Each ERAM Proc number is %d \n", *arnoldi_proc);
      printf("Initialization ]> Thus Total proc number = 1 + 1 + %d * %d + %d * %d = %d \n", *gmres_nb, *gmres_proc, *arnoldi_nb, *arnoldi_proc,2 + *gmres_nb * *gmres_proc + *arnoldi_nb * *arnoldi_proc);
    }

  return 0;
}

char **argsParserGmresExec(int argc, char ** argv, int gmres_nb){

  int cnt5 = 0;
  int i, j;
  char **gmres_cmds = new char* [gmres_nb];

  //string parser

    for(i =0; i < argc; i++){
      //parse gmres executable name char
      if (strcmp(argv[i],"-gmres_exec") == 0){
        cnt5 = 1;
        if(i + gmres_nb < argc){
          for(j = 1; j <= gmres_nb; j++){
            gmres_cmds[j - 1] = new char [20];
            strcpy(gmres_cmds[j - 1], argv[i + j]);
          }
        }
      }
    }
    //if well parsered, show the initailization information
    if(cnt5){

      printf("Initialization ]> GMRES executable name list is: [ ");
      for(j = 0; j < gmres_nb; j++){
        printf("%s  ", gmres_cmds[j]);
        if(j == gmres_nb - 1){
          printf("]\n");
        }
      }
    }
  return gmres_cmds;
}

char **argsParserArnoldiExec(int argc, char ** argv, int arnoldi_nb){

  int cnt6 = 0;
  int i, j;
  char **arnoldi_cmds = new char* [arnoldi_nb]; //arnoldi executable name

  //string parser

    for(i =0; i < argc; i++){
      //parse arnoldi executable name char
      if (strcmp(argv[i],"-arnoldi_exec") == 0){
        cnt6 = 1;
        if(i + arnoldi_nb < argc){
          for(j = 1; j <= arnoldi_nb; j++){
            arnoldi_cmds[j - 1] = new char [20];
            strcpy(arnoldi_cmds[j - 1], argv[i + j]);
          }
        }
      }
    }
    //if well parsered, show the initailization information
    if(cnt6){

      printf("Initialization ]> ERAM executable name list is: [ ");
      for(j = 0; j < arnoldi_nb; j++){
        printf("%s  ", arnoldi_cmds[j]);
        if(j == arnoldi_nb - 1){
          printf("]\n");
        }
      }
    }
  return arnoldi_cmds;
}

char *argsParserLsqrExec(int argc, char ** argv){
  int cnt7 = 0;
  int i;
  char *lsqr_cmd = new char [20];//lsqr executable name

  //string parser

    for(i =0; i < argc; i++){
      //parse lsqr executable name char
      if (strcmp(argv[i],"-lsqr_exec") == 0){
        cnt7 = 1;
        if(i + 1 < argc){
          strcpy(lsqr_cmd, argv[i + 1]);
        }
      }
    }
    //if well parsered, show the initailization information
    if(cnt7){

      printf("Initialization ]> LS executable name list is: [ ");
      printf("%s  ", lsqr_cmd);
      printf("]\n");
    }

  return lsqr_cmd;

}

char **argsParserArnoldiRuntime(int argc, char ** argv){
  
  int i, j = 0;

  char **arnoldi_runtime_args = new char* [30];
  char *skr = new char [100];
  char *eps = new char [100];

  for(i = 0; i < argc; i++){
    //parse matrix market file name
    memcpy(skr, argv[i] + 0, 10);
    if (strcmp(skr,"--filename") == 0){
      arnoldi_runtime_args[j] = new char [100];
      strcpy(arnoldi_runtime_args[j], argv[i]);
      printf("Parser Matrix filename Flag = %s  \n", arnoldi_runtime_args[j]);
      j++;
    }
  }

  for(i = 0; i < argc; i++){
    //parse arnoldi runtime parameters
    memcpy(eps, argv[i] + 0, 5);
    if (strcmp(eps,"--eps") == 0){
      arnoldi_runtime_args[j] = new char [100];
      strcpy(arnoldi_runtime_args[j], argv[i]);
      printf("Parser Arnoldi runtime Flag = %s  \n", arnoldi_runtime_args[j]);
      j++;
    }
  }
  

  delete [] skr;
  delete [] eps;

  return arnoldi_runtime_args;

}

char **argsParserGMRESRuntime(int argc, char ** argv){
  
  int i, j = 0;

  char **gmres_runtime_args = new char* [30];
  char *skr = new char [100];
  char *ksp = new char [100];

  for(i = 0; i < argc; i++){
    //parse arnoldi runtime parameters
    memcpy(ksp, argv[i] + 0, 6);
    if (strcmp(ksp,"--ksp-") == 0){
      gmres_runtime_args[j] = new char [100];
      strcpy(gmres_runtime_args[j], argv[i]);
      printf("Parser GMRES runtime Flag = %s  \n", gmres_runtime_args[j]);
      j++;
    }
  }
  

  delete [] skr;
  delete [] ksp;

  return gmres_runtime_args;
  
}

char **argsParserLSPRuntime(int argc, char ** argv){
  
  int i, j = 0;

  char **lsp_runtime_args = new char* [30];
  char *lsp = new char [100];

  for(i = 0; i < argc; i++){
    //parse arnoldi runtime parameters
    memcpy(lsp, argv[i] + 0, 5);
    if (strcmp(lsp,"--lsp") == 0){
      lsp_runtime_args[j] = new char [100];
      strcpy(lsp_runtime_args[j], argv[i]);
      printf("Parser Least Square Polynomial runtime Flag = %s  \n", lsp_runtime_args[j]);
      j++;
    }
  }
  

  delete [] lsp;

  return lsp_runtime_args;
  
}

#endif
