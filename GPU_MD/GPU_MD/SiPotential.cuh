
#ifndef SIPOTENTIAL_CU_H_
#define SIPOTENTIAL_CU_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include "Calculations.cuh"
#include "SiPotential.cuh"
#include <math.h>
#include <fstream>
#include <cuda.h>//cuda
#include <cuda_runtime.h>//cuda
#include <device_launch_parameters.h>//cuda

//------------- potential between two Si ------------------// 
__device__ double d_f2(double);

__device__ double d_v2(double);
//--------------------------------------------------------//

__device__ double d_f2_derivative_of_rij_tag(double);

__device__ double d_v2_derivative_of_rix(real3, real3, double);
__device__ double d_v2_derivative_of_riy(real3, real3, double);
__device__ double d_v2_derivative_of_riz(real3, real3, double);
//--------------------------------------------------------//

//--------- force between three Si particles --------------//
__device__ double d_hi_derivative_of_rij_tag(double, double, double);
__device__ double d_hi_derivative_of_rik_tag(double, double, double);
__device__ double d_hj_derivative_of_rij_tag(double, double, double);
__device__ double d_hj_derivative_of_rik_tag(double, double, double);
__device__ double d_hk_derivative_of_rij_tag(double, double, double);
__device__ double d_hk_derivative_of_rik_tag(double, double, double);

__device__ double d_f3_derivative_of_rij_tag(double, double, double);
__device__ double d_f3_derivative_of_rik_tag(double, double, double);

__device__ double d_v3_derivative_of_rix(real3,real3,real3, double, double, double);
__device__ double d_v3_derivative_of_riy(real3,real3,real3, double, double, double);
__device__ double d_v3_derivative_of_riz(real3,real3,real3, double, double, double);
//--------------------------------------------------------//
//------------- potential between two Si ------------------// 
__device__ double d_hi(double, double, double);
__device__ double d_hj(double, double, double);
__device__ double d_hk(double, double, double);

__device__ double d_f3(double, double, double);

__device__ double d_v3(double, double, double);
//--------------------------------------------------------//

#endif //SIPOTENTIAL_CU_H_