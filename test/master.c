#include "mpi.h"

int main(int argc, char *argv[])
{ 
   int n_spawns = 2;
   MPI_Comm intercomm;

   MPI_Init(&argc, &argv);

   MPI_Comm_spawn("worker_program", MPI_ARGV_NULL, n_spawns, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE); 

   int sendbuf[2] = {3, 5};
   int recvbuf; // redundant for master.

   MPI_Scatter(sendbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, MPI_ROOT, intercomm);

   MPI_Finalize();
   return 0;
}
