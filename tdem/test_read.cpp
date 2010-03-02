/*
 * =====================================================================================
 *
 *       Filename:  multithread.cpp
 *
 *    Description:  multi thread
 *
 *        Version:  1.0
 *        Created:  24/02/10 15:08:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
*/
#include "hdf5.h"
#include "hdf5_hl.h"
#include <iostream>

#define RANK 2

using namespace std;

int main( void )
{

 hid_t       file_id;
 double      *data;
 hsize_t     dims[1];
 herr_t      status;
 size_t     i, j, nrow, n_values;

 /* open file from ex_lite1.c */
 file_id = H5Fopen ("muu.msys", H5F_ACC_RDONLY, H5P_DEFAULT);

 /* read dataset */
 //status = H5LTread_dataset_double(file_id,"/Particle/cod",data);

 /* get the dimensions of the dataset */
 status = H5LTget_dataset_info(file_id,"/Particle/cod",dims,NULL,NULL);

 cout << dims[0] << endl;

 data = new double[5];
 status = H5LTread_dataset_double(file_id,"/Particle/cod",data);
 cout << data[2] << endl;


 /* print it by rows */
 //n_values = (size_t)(dims[0] * dims[1]);
 //nrow = (size_t)dims[0];
 //cout << nrow << endl;

 //for (i=0; i<n_values/nrow; i++ )
 //{
  //for (j=0; j<nrow; j++)
   //printf ("  %d", data[i*nrow + j]);
  //printf ("\n");
 //}

 /* close file */
 status = H5Fclose (file_id);

 return 0;

}
