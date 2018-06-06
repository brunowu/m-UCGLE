#include "args_parser.h"

// -nb_gmres ${GNB} -proc_gmres ${GPNB}
/*function to parser the args to setup the method*/
int argsParser(int argc, char ** argv, int *gmres_nb, int *arnoldi_nb, int *gmres_proc, int *arnoldi_proc, char **cmd){

  int i, j;
  //nb of gmres and arnoldi spawned by the father process

  int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0, cnt5 = 0, cnt6 = 0;

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

//string parser
  char gmres_cmds[*gmres_nb][20]; //gmres executable name
  char arnoldi_cmds[*arnoldi_nb][20]; //arnoldi executable name

    for(i =0; i < argc; i++){
      //parse gmres executable name char
      if (strcmp(argv[i],"-gmres_exec") == 0){
        cnt5 = 1;
        if(i + *gmres_nb < argc){
          for(j = 1; j <= *gmres_nb; j++){
              strcpy(gmres_cmds[j - 1], argv[i + j]);
          }
        }
      }
      else if (strcmp(argv[i],"-arnoldi_exec") == 0){
        cnt6 = 1;
        if(i + *arnoldi_nb < argc){
          for(j = 1; j <= *arnoldi_nb; j++){
            strcpy(arnoldi_cmds[j - 1], argv[j + i]);
          }
        }
      }
    }

    //if well parsered, show the initailization information
    if(cnt1 && cnt2 && cnt3 && cnt4 && cnt5 && cnt6){
      printf("Initialization ]> GMRES executable number is %d \n", *gmres_nb);
      printf("Initialization ]> Each GMRES Proc number is %d \n", *gmres_proc);
      printf("Initialization ]> ERAM executable number is %d \n", *arnoldi_nb);
      printf("Initialization ]> Each ERAM Proc number is %d \n", *arnoldi_proc);
      printf("Initialization ]> Thus Total proc number = 1 + 1 + %d * %d + %d * %d = %d \n", *gmres_nb, *gmres_proc, *arnoldi_nb, *arnoldi_proc,2 + *gmres_nb * *gmres_proc + *arnoldi_nb * *arnoldi_proc);

      printf("Initialization ]> GMRES executable name list is: [ ");
      for(j = 0; j < *gmres_nb; j++){
        printf("%s  ", gmres_cmds[j]);
        if(j == *gmres_nb - 1){
          printf("]\n");
        }
      }

      printf("Initialization ]> ERAM executable name list is: [ ");
      for(j = 0; j < *arnoldi_nb; j++){
        printf("%s  ", arnoldi_cmds[j]);
        if(j == *arnoldi_nb - 1){
          printf("]\n");
        }
      }

    }

    cmd = (char **)malloc(1*sizeof(char[20]));

  return 0;
}

int main(int argc, char ** argv){

  int gmres_nb, arnoldi_nb, gmres_proc, arnoldi_proc;
  char *cmd[20];

  MPI_Init (&argc, &argv);

  argsParser(argc, argv, &gmres_nb, &arnoldi_nb, &gmres_proc, &arnoldi_proc, cmd);
  printf("main: gmres_nb = %d, arnoldi_nb = %d, gmres_proc = %d, arnoldi_proc = %d\n",gmres_nb, arnoldi_nb, gmres_proc, arnoldi_proc);


  MPI_Finalize();

  return 0;
}
