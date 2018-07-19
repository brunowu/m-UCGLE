#ifndef _MPI_LSA_COM_H_
#define _MPI_LSA_COM_H_

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <complex>

/*recv type information functionality*/
int mpi_lsa_com_type_recv(MPI_Comm * com, int * type){
  //check if any type to receive
  int flag, count = 0;
  MPI_Status status;
  MPI_Request request;

  MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, *com, &flag, &status);

  if(!flag)
    return 1;

  MPI_Get_count(&status,MPI_INT,&count);

  printf("count = %d\n",count);

  if(count != 1) return 1;

  MPI_Recv(type, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, * com, &status);

  printf("type == %d\n", *type);

  return 0;
}
/*send type information functionality*/
int mpi_lsa_com_type_send(MPI_Comm * com, int * type){
  int remote_size, flag;

  MPI_Comm_remote_size(* com, &remote_size);

  MPI_Request aReq[remote_size];
  MPI_Status status;

  int i;

  //check if previous requests where completed, if no, cancel it then
  

/*
  for(i = 0; i < remote_size; i++){
    MPI_Test(&aReq[i],&flag,&status);
    // if not cancel it
    if(!flag){
      MPI_Cancel(&aReq[i]);
    }
  }
 */ 
  for(i = 0; i < remote_size; i++){
    MPI_Isend(type, 1, MPI_INT, i, i, *com, &aReq[i]);
    MPI_Wait(&aReq[i], &status);
  }

  printf("hello2\n");

  return 0;
}

/*send double array functionality*/
int mpi_lsa_com_array_send(MPI_Comm * com, int * size, double * data, MPI_Request * array_Req){
  int remote_size;
  MPI_Comm_remote_size(* com, &remote_size);

  int i;

//  MPI_Request array_Req[remote_size];
  MPI_Status status;

  double *array_out_sended_buffer = new double [*size];

  for(i = 0; i < * size; i++){
    array_out_sended_buffer[i] = data[i];
  }

  for(i = 0; i < remote_size; i++){
    MPI_Isend(array_out_sended_buffer, *size, MPI_DOUBLE, i, i, *com, &array_Req[i]);
    //MPI_Wait(&array_Req[i], &status);
  }

  return 0;
}

/*receive double array functionality*/
int mpi_lsa_com_array_recv(MPI_Comm * com, int * size, double * data){
  /*
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

  */

  MPI_Status status;
  int flag;
  MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,*com,&flag,&status);
  MPI_Get_count(&status,MPI_INT,size);

  if(!flag || *size==1){
    return 1;
  }

  MPI_Get_count(&status,MPI_DOUBLE,size);

  MPI_Recv(data, *size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, * com, &status);

  return 0;
}


/*send complex double array functionality*/
int mpi_lsa_com_cplx_array_send(MPI_Comm * com, int * size, std::complex<double> * data){

  int remote_size;

  MPI_Comm_remote_size(* com, &remote_size);

  int i , j;

  MPI_Request array_Req[remote_size];
  MPI_Status status;

  double *array_out_sended_buffer = new double [*size*2];

  for(i = 0; i < * size*2; i = i+2){
    j = (int) i / 2;
    array_out_sended_buffer[i] = data[j].real();
    array_out_sended_buffer[i + 1] = data[j].imag();
  }

  for(i = 0; i < remote_size; i++){
    MPI_Isend(array_out_sended_buffer, *size*2, MPI_DOUBLE, i, i, *com, &array_Req[i]);
    //MPI_Wait(&array_Req[i], &status);
  }

  return 0;
}

/*receive complex double array functionality*/
int mpi_lsa_com_cplx_array_recv(MPI_Comm * com, int * size, std::complex<double> * data){
  int flag = 0;

  int double_size;

  MPI_Status status;
  MPI_Request request;

  int i, j;

  while(!flag){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, *com, &flag, &status);
  }

  if(flag){

    MPI_Get_count(&status,MPI_DOUBLE,&double_size);

    *size = (int) double_size/2;

    if(*size == 1){
      return 1;
    }

    //how large the array to receive is
    MPI_Get_count(&status,MPI_DOUBLE,&double_size);

    *size = (int) double_size/2;

    double *array_out_recv_buffer = new double [double_size];

    if(*size >= 1){
      MPI_Recv(array_out_recv_buffer,double_size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, * com, &status);
      
      for(i = 0; i < *size; i++ ){
        j = (int) 2 * i;
        data[i].real(array_out_recv_buffer[j]);
        data[i].imag(array_out_recv_buffer[j+1]);
      }
    }
  }

  return 0;
}

#endif
