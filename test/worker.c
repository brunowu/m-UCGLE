#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{  
   MPI_Init(&argc, &argv);

   MPI_Comm intercomm; 
   MPI_Comm_get_parent(&intercomm);

   int sendbuf[2]; // redundant for worker.
   int recvbuf;

   MPI_Scatter(sendbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, 0, intercomm);
   printf("recvbuf = %d\n", recvbuf);

   MPI_Finalize();
   return 0;
}

