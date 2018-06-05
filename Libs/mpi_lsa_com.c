#include "mpi_lsa_com.h"

/*send type information functionality*/
int mpi_lsa_com_type_recv(MPI_Comm * com, int * type){
  //check if any type to receive
  int flag = 0, count;
  MPI_Status status;
  MPI_Request request;

  while(!flag){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, *com, &flag, &status);
  }

  if(flag){
    MPI_Get_count(&status,MPI_INT,&count);

    if(count == 1){
      MPI_Recv(type, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, * com, &status);
    }
  }

  return 0;
}
/*receive type information functionality*/
int mpi_lsa_com_type_send(MPI_Comm * com, int * type){
  MPI_Request aReq[2];
  MPI_Status status;
  for(int i = 0; i < 2; i++){
    MPI_Isend(type, 1, MPI_INT, i, i,*com, &aReq[i]);
    MPI_Wait(&aReq[i], &status);
  }

  return 0;
}
