/*
#ifndef SIPOTENTIAL_H_
#define SIPOTENTIAL_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include "Calculations.h"
#include "SiPotential.h"
#include <fstream>
#include <cuda.h>//cuda
#include <cuda_runtime.h>//cuda
#include <device_launch_parameters.h>//cuda

//------------- potential between two Si ------------------// 
__host__ __device__ double f2(double);

__host__ __device__ double v2(double);
//--------------------------------------------------------//

__host__ __device__ double f2_derivative_of_rij_tag(double);

__host__ __device__ double v2_derivative_of_rix(real3, real3, double);
__host__ __device__ double v2_derivative_of_riy(real3, real3, double);
__host__ __device__ double v2_derivative_of_riz(real3, real3, double);
//--------------------------------------------------------//

//--------- force between three Si particles --------------//
__host__ __device__ double hi_derivative_of_rij_tag(double, double, double);
__host__ __device__ double hi_derivative_of_rik_tag(double, double, double);
__host__ __device__ double hj_derivative_of_rij_tag(double, double, double);
__host__ __device__ double hj_derivative_of_rik_tag(double, double, double);
__host__ __device__ double hk_derivative_of_rij_tag(double, double, double);
__host__ __device__ double hk_derivative_of_rik_tag(double, double, double);

__host__ __device__ double f3_derivative_of_rij_tag(double, double, double);
__host__ __device__ double f3_derivative_of_rik_tag(double, double, double);

__host__ __device__ double v3_derivative_of_rix(real3,real3,real3, double, double, double);
__host__ __device__ double v3_derivative_of_riy(real3,real3,real3, double, double, double);
__host__ __device__ double v3_derivative_of_riz(real3,real3,real3, double, double, double);
//--------------------------------------------------------//
//------------- potential between two Si ------------------// 
__host__ __device__ double hi(double, double, double);
__host__ __device__ double hj(double, double, double);
__host__ __device__ double hk(double, double, double);

__host__ __device__ double f3(double, double, double);

__host__ __device__ double v3(double, double, double);
//--------------------------------------------------------//

#endif //SIPOTENTIAL_H_
*/