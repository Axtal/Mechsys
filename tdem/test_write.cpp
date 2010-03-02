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
#include <vector>

using namespace std;

class numero
{
public:
    int n;
    vector<int> V;
};

int main( void )
{

hid_t       file_id, group_id;
 hsize_t     dims[1]={1},dims2[1]={5};
 int         data[1]={5};
 double      cod[5]={1.5,123.1,10.5,1.00003,2.0};
 herr_t      status;

 /* create a HDF5 file */
 file_id = H5Fcreate ("muu.msys", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

 group_id = H5Gcreate(file_id, "/Particle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

 /* create and write an integer type dataset named "dset" */
 status = H5LTmake_dataset(file_id,"/Particle/nv",1,dims,H5T_NATIVE_INT,data);
 status = H5LTmake_dataset(file_id,"/Particle/cod",1,dims2,H5T_NATIVE_DOUBLE,cod);

 /* close file */
 H5Gclose(group_id);
 status = H5Fclose (file_id);

 return 0;
}
