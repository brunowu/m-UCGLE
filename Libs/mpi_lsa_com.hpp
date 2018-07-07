#ifndef _MPI_LSA_COM_H_
#define _MPI_LSA_COM_H_

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <complex>

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
int mpi_lsa_com_type_send(MPI_Comm * com, int * type, int *count){
  int remote_size, flag;
  MPI_Comm_remote_size(* com, &remote_size);

  MPI_Request aReq[remote_size];
  MPI_Status status;

  int i;

  //check if previous requests where completed, if no, cancel it then
  for(i = 0; i < *count; i++){
    MPI_Test(&aReq[i],&flag,&status);
    // if not cancel it
    if(!flag){
      MPI_Cancel(&aReq[i]);
    }
  }
  for(i = 0; i < remote_size; i++){
    MPI_Isend(type, 1, MPI_INT, i, i,*com, &aReq[i]);
    MPI_Wait(&aReq[i], &status);
    *count = *count+1;
  }

  return 0;
}

/*send double array functionality*/
int mpi_lsa_com_array_send(MPI_Comm * com, int * size, double * data){
  int remote_size;
  MPI_Comm_remote_size(* com, &remote_size);

  int i;

  MPI_Request array_Req[remote_size];
  MPI_Status status;

  double *array_out_sended_buffer = new double [*size];

  for(i = 0; i < * size; i++){
    array_out_sended_buffer[i] = data[i];
  }

  for(i = 0; i < remote_size; i++){
    MPI_Isend(array_out_sended_buffer, *size, MPI_DOUBLE, i, i, *com, &array_Req[i]);
    MPI_Wait(&array_Req[i], &status);
  }

  return 0;
}

/*receive double array functionality*/
int mpi_lsa_com_array_recv(MPI_Comm * com, int * size, double * data){
  //check if any array to receive
  int flag = 0;
  MPI_Status status;
  MPI_Request request;
  int i;


  while(!flag){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, *com, &flag, &status);
  }

  if(flag){

    MPI_Get_count(&status,MPI_INT,size);
    if(*size == 1){
      return 1;
    }
  
    //how large the array to receive is
    MPI_Get_count(&status,MPI_DOUBLE,size);

    if(*size >= 1){
      MPI_Recv(data, *size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, * com, &status);
    }
  }

  return 0;
}


/*send complex double array functionality*/
int mpi_lsa_com_cplx_array_send(MPI_Comm * com, int * size, std::complex<double> * data){

  int remote_size;

  MPI_Comm_remote_size(* com, &remote_size);

  int i;

  MPI_Request array_Req[remote_size];
  MPI_Status status;

  double *array_out_sended_buffer = new double [*size*2];

  for(i = 0; i < * size*2; i++){
    array_out_sended_buffer[i] = data[i].real();
    array_out_sended_buffer[i + 1] = data[i].imag();
    i = i + 2;
  }

  for(i = 0; i < remote_size; i++){
    MPI_Isend(array_out_sended_buffer, *size*2, MPI_DOUBLE, i, i, *com, &array_Req[i]);
    MPI_Wait(&array_Req[i], &status);
  }

  return 0;
}

/*receive complex double array functionality*/
int mpi_lsa_com_cplx_array_recv(MPI_Comm * com, int * size, std::complex<double> * data){
  //check if any array to receive
  int flag = 0;

  MPI_Status status;
  MPI_Request request;
  int i;

  while(!flag){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, *com, &flag, &status);
  }

  if(flag){

    MPI_Get_count(&status,MPI_DOUBLE,size);
    if(*size == 1){
      return 1;
    }
  
    //how large the array to receive is
    MPI_Get_count(&status,MPI_COMPLEX,size);

    if(*size >= 1){
      MPI_Recv(data, *size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, * com, &status);
    }
  }

  return 0;
}

#endif
