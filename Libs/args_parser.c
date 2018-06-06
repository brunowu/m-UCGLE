#include "args_parser.h"

// -nb_gmres ${GNB} -proc_gmres ${GPNB}
/*function to parser the args to setup the method*/
int argsParser(int argc, char ** argv){

  int i, j, rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //nb of gmres and arnoldi spawned by the father process
  int gmres_nb, arnoldi_nb;
  int gmres_proc, arnoldi_proc;

  for(i =0; i < argc; i++){
    //nb gmres executable
    if(strcmp(argv[i],"-nb_gmres") == 0){
        if(i + 1 < argc - 1){
          i++;
          gmres_nb = atoi(argv[i]);
        }
      }
      //nb each gmres proc
      else if(strcmp(argv[i],"-proc_gmres") == 0){
        if(i + 1 < argc - 1){
          i++;
          gmres_proc = atoi(argv[i]);
        }
      }
      //nb arnoldi executable
      else if(strcmp(argv[i],"-nb_arnoldi") == 0){
        if(i + 1 < argc - 1){
          i++;
          arnoldi_nb = atoi(argv[i]);
        }
      }
      //nb of proc for each arnoldi
      else if (strcmp(argv[i],"-proc_arnoldi") == 0){
        if(i + 1 < argc - 1){
          i++;
          arnoldi_proc = atoi(argv[i]);
        }
      }
    }

//string parser
  char gmres_cmds[][20]; //gmres executable name
  char arnoldi_cmds[][20]; //arnoldi executable name

    for(i =0; i < argc; i++){
      //parse gmres executable name char
      if (strcmp(argv[i],"-gmres_exec") == 0){
        if(i + gmres_nb < argc - 1){
          for(j = 1; j <= gmres_nb; j++){
            gmres_cmds[j - 1] = argv[j + gmres_nb];
          }
        }
      }
      else if (strcmp(argv[i],"-arnoldi_exec") == 0){
        if(i + arnoldi_nb < argc - 1){
          for(j = 1; j <= gmres_nb; j++){
            arnoldi_cmds[j - 1] = argv[j + arnoldi_nb];
          }
        }
      }
    }
  return 0;
}
