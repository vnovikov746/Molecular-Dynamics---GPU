/*
 *
 */

#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include "LennardJones.h"
#include "SiPotential.h"
#include <fstream>
#include <cuda.h>//cuda
#include <cuda_runtime.h>//cuda
#include <device_launch_parameters.h>//cuda

//---------- force between two Si particles ---------------//
__host__ __device__ double distance2(real3, real3);//r(i,j)

//-------------- calculate forces between Si --------------//
void calculateForce_Si(int, int, particleStruct*, particleStruct*, int, int, bool, bool);
//--------------------------------------------------------//

void calculateForce_Xe(int, int, particleStruct*, particleStruct*, int, int, bool);

//---------calculate total potential of Si ---------------//
double V_total_Si(int, particleStruct*);
double V_total_Si_Xe(int, particleStruct*, int, particleStruct*);
double V_total_Xe(int, particleStruct*);
void V_total(int, particleStruct*, int, particleStruct*);
//--------------------------------------------------------//

void initiateAcceleration(particleStruct*, int, double);

void predict(particleStruct*, int, float);

void correct(particleStruct*, double, int, double);

//-------------- calculate forces between Si --------------//
__global__ void d_calculateForce_Si(int, int, particleStruct*, particleStruct*, particleStruct*, int, int, bool, bool);
//--------------------------------------------------------//

__global__ void d_calculateForce_Xe(int, int, particleStruct*, particleStruct*, particleStruct*, int, int, bool);

//---------calculate total potential of Si ---------------//
//__global__ void d_V_total_Si(int, particleStruct*, double*);
//__global__ void d_V_total_Si_Xe(int, particleStruct*, int, particleStruct*, double*);
//__global__ void d_V_total_Xe(int, particleStruct*, double*);
//__global__ void d_V_total(int, particleStruct*, int, particleStruct*);
//--------------------------------------------------------//

__global__ void d_initiateAcceleration(particleStruct*, int, double);

__global__ void d_predict(particleStruct*, int, float);

__global__ void d_correct(particleStruct*, double, int, double);


//-----------------Lennard Jones Force -------------------//
__host__ __device__ double lennardJonesForce(double, double, double);
//--------------------------------------------------------//
__host__ __device__ double lennardJonesPotential(double, double, double);



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


#endif	//CALCULATIONS_H_