#ifndef LENNARDJONES_CU_H_
#define LENNARDJONES_CU_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include <math.h>
#include <fstream>
#include <cuda.h>//cuda
#include <cuda_runtime.h>//cuda
#include <device_launch_parameters.h>//cuda

//-----------------Lennard Jones Force -------------------//
__device__ double d_lennardJonesForce(double, double, double);
//--------------------------------------------------------//
__device__ double d_lennardJonesPotential(double, double, double);

#endif //LENNARDJONES_CU_H_