/*
 *
 */

#include "Calculations.h"

__host__ __device__ double distance2(real3 i, real3 j)
{
	return sqrt((i.x-j.x)*(i.x-j.x) + (i.y-j.y)*(i.y-j.y) + (i.z-j.z)*(i.z-j.z));
}

//-------------------------- calculate Force Si ------------------------------//
__global__ void
__launch_bounds__(1024, 4)
d_calculateForce_Si(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, particleStruct* siParticles, particleStruct* siParticles2, particleStruct* xeParticles, int numOfSi, int numOfXe, bool USE_NEIGHBOR_LISTS, bool useLennardJonesPotentialForSi)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < numOfXe)
	{
		real3 iPosition;
		real3 jPosition;
		real3 kPosition;
		double r_ij = 0.0;
		double r_ik = 0.0;
		double r_jk = 0.0;
		iPosition = siParticles[idx].position;
		siParticles[idx].force.x = 0.0;
		siParticles[idx].force.y = 0.0;
		siParticles[idx].force.z = 0.0;
		for(int j = 0; j < numOfSi; j++)
		{
			if(j != idx)
			{
				jPosition = siParticles2[j].position;
				r_ij = distance2(iPosition, jPosition);
				if(r_ij/sigma_Si < a_Si)
				{
					siParticles[idx].force.x -= v2_derivative_of_rix(iPosition, jPosition, r_ij);
					siParticles[idx].force.y -= v2_derivative_of_riy(iPosition, jPosition, r_ij);
					siParticles[idx].force.z -= v2_derivative_of_riz(iPosition, jPosition, r_ij);
				}
				for(int k = 0; k < numOfSi; k++)
				{
					if(k != idx && k != j)
					{
						kPosition = siParticles2[k].position;							

						r_ik = distance2(iPosition, kPosition);
						r_jk = distance2(jPosition, kPosition);
						if((r_ij/sigma_Si < a_Si && r_ik/sigma_Si < a_Si) || (r_ij/sigma_Si < a_Si && r_jk/sigma_Si < a_Si) || (r_ik/sigma_Si < a_Si && r_jk/sigma_Si < a_Si))
						{
							siParticles[idx].force.x -= v3_derivative_of_rix(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[idx].force.x -= v3_derivative_of_rix(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
							siParticles[idx].force.y -= v3_derivative_of_riy(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[idx].force.y -= v3_derivative_of_riy(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
							siParticles[idx].force.z -= v3_derivative_of_riz(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[idx].force.z -= v3_derivative_of_riz(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
						}
					}
				}
			}
		}
	}
	/*	for(int j = 0; j < countXe; j++)
		{
			jPosition = xeParticles[j].position;

			jPosition.x += 0.25*config->SI_LENGTH*space_Si;
			jPosition.y += 0.25*config->SI_LENGTH*space_Si;
			jPosition.z += config->SI_HEIGHT+config->LA_SPACE;
			r_ij = distance2(iPosition, jPosition);
			if(r_ij/sigma_Si_Xe < a_Si_Xe)
			{
				siParticles[i].force.x += ((iPosition.x-jPosition.x)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
				siParticles[i].force.y += ((iPosition.y-jPosition.y)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
				siParticles[i].force.z += ((iPosition.z-jPosition.z)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
			}
		}*/
//	}
}
//----------------------------------------------------------------------------//

__global__ void 
__launch_bounds__(1024, 4)
d_calculateForce_Xe(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, particleStruct* xeParticles, particleStruct* xeParticles2, particleStruct* siParticles, int numOfSi, int numOfXe, bool USE_NEIGHBOR_LISTS)
{
//	extern __shared__ particleStruct* sharedXe[];

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < numOfXe)
	{
		real3 iPosition;
		real3 jPosition;
		double r_ij = 0.0;
		iPosition = xeParticles[idx].position;
		xeParticles[idx].force.x = 0.0;
		xeParticles[idx].force.y = 0.0;
		xeParticles[idx].force.z = 0.0;
		for(int j = 0; j < numOfXe; j++)
		{
			if(j != idx)
			{
				jPosition = xeParticles2[j].position;
				r_ij = distance2(iPosition, jPosition);
				if(r_ij/sigma_Xe_Xe < xe_Cluster)
				{
					xeParticles[idx].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
					xeParticles[idx].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
					xeParticles[idx].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
				}
			}
		}
	}
	/*	for(int j = 0; j < countSi; j++)
		{
			jPosition = siParticles[j].position;

			r_ij = distance2(iPosition, jPosition);
			if(r_ij/sigma_Si_Xe < a_Si_Xe)
			{
				xeParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
				xeParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
				xeParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
			}
		}*/
//	}
}

//////////////////////////////////////////////////////////////////////////////////////////////!!!!
//----------------------calculate total potential of Si -----------------------//
/*__global__ double d_V_total_Si(int SI_PARTICLES, particleStruct* siParticles)
{

}
//----------------------------------------------------------------------------//
__global__ double d_V_total_Si_Xe(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{

}

__global__ double d_V_total_Xe(int XE_PARTICLES, particleStruct* xeParticles)
{

}

__global__ void d_V_total(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{

}*/
//////////////////////////////////////////////////////////////////////////////////////////////!!!!
__global__ void d_initiateAcceleration(particleStruct *particles, int listSize, double mass)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < listSize)
	{
		particles[idx].aAcc.x = particles[idx].force.x/mass;
		particles[idx].aAcc.y = particles[idx].force.y/mass;
		particles[idx].aAcc.z = particles[idx].force.z/mass;
	}
}

__global__ void d_predict(particleStruct *particles, int listSize, float dt)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;

	__shared__ double c1;
    __shared__ double c2;
    __shared__ double c3;
    __shared__ double c4;

	c1 = dt;
	c2 = c1*dt/2.0;
	c3 = c2*dt/3.0;
	c4 = c3*dt/4.0;

	if(idx < listSize)
	{
		particles[idx].position.x += c1*particles[idx].velocity.x + c2*particles[idx].aAcc.x + c3*particles[idx].bAcc.x + c4*particles[idx].cAcc.x;
		particles[idx].position.y += c1*particles[idx].velocity.y + c2*particles[idx].aAcc.y + c3*particles[idx].bAcc.y + c4*particles[idx].cAcc.y;
		particles[idx].position.z += c1*particles[idx].velocity.z + c2*particles[idx].aAcc.z + c3*particles[idx].bAcc.z + c4*particles[idx].cAcc.z;
		particles[idx].velocity.x += c1*particles[idx].aAcc.x + c2*particles[idx].bAcc.x + c3*particles[idx].cAcc.x;
		particles[idx].velocity.y += c1*particles[idx].aAcc.y + c2*particles[idx].bAcc.y + c3*particles[idx].cAcc.y;
		particles[idx].velocity.z += c1*particles[idx].aAcc.z + c2*particles[idx].bAcc.z + c3*particles[idx].cAcc.z;
		particles[idx].aAcc.x += c1*particles[idx].bAcc.x + c2*particles[idx].cAcc.x;
		particles[idx].aAcc.y += c1*particles[idx].bAcc.y + c2*particles[idx].cAcc.y;
		particles[idx].aAcc.z += c1*particles[idx].bAcc.z + c2*particles[idx].cAcc.z;
		particles[idx].bAcc.x += c1*particles[idx].cAcc.x;
		particles[idx].bAcc.y += c1*particles[idx].cAcc.y;
		particles[idx].bAcc.z += c1*particles[idx].cAcc.z;
	}
}

__global__ void d_correct(particleStruct *particles, double dt, int listSize, double mass)
{
	__shared__ double c1;
	__shared__ double c2;
	__shared__ double c3;
	__shared__ double c4;

	__shared__ double cr;
	__shared__ double cv;
	__shared__ double cb;
	__shared__ double cc;

	c1 = dt ;
	c2 = c1*dt/2.0;
	c3 = c2*dt/3.0; 
	c4 = c3*dt/4.0;

	cr = GEAR1*c2;
	cv = GEAR2*c2/c1;
	cb = GEAR3*c2/c3;
	cc = GEAR4*c2/c4;

	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx < listSize)
	{
		double axi = particles[idx].force.x/mass;
      	double ayi = particles[idx].force.y/mass;
      	double azi = particles[idx].force.z/mass;

		double corrx = axi - particles[idx].aAcc.x;
    	double corry = ayi - particles[idx].aAcc.y;
    	double corrz = azi - particles[idx].aAcc.z;

		particles[idx].position.x += cr*corrx;
		particles[idx].position.y += cr*corry;
		particles[idx].position.z += cr*corrz;
		particles[idx].velocity.x += cv*corrx;
		particles[idx].velocity.y += cv*corry;
		particles[idx].velocity.z += cv*corrz;
		particles[idx].aAcc.x = axi;
		particles[idx].aAcc.y = ayi;
		particles[idx].aAcc.z = azi;
		particles[idx].bAcc.x += cb*corrx;
		particles[idx].bAcc.y += cb*corry;
		particles[idx].bAcc.z += cb*corrz;
		particles[idx].cAcc.x += cc*corrx;
		particles[idx].cAcc.y += cc*corry;
    	particles[idx].cAcc.z += cc*corrz;
	}
}


//------------------------Lennard Jones Potential -----------------------------//
__host__ __device__ double lennardJonesForce(double dist, double sig, double eps)
{
	double sigsq = sig*sig;
	double con = 24.0*eps/sigsq;
	double dist2 = dist * dist;
	dist2 /= sigsq;
	
	double dist4 = dist2*dist2;
	double dist8 = dist4*dist4;
	double dist14 = dist2*dist4*dist8;
	double invdist8= 1.0/dist8;
	double invdist14= 1.0/dist14;
	double s = 2.0f*invdist14-invdist8;
	return s * con;
}
//----------------------------------------------------------------------------//

__host__ __device__ double lennardJonesPotential(double dist, double sig, double eps)
{
	double expr = sig/dist;
	double expr2 = expr*expr;
	double expr4 = expr2*expr2;
	double expr6 = expr4*expr2;
	double expr12 = expr6*expr6;

	return 4.0*eps*(expr12-expr6);
}




//-------------------- force between two Si particles ---------------------//
__host__ __device__ double f2_derivative_of_rij_tag(double r_ij_tag)
{
	double r_ij_tag_minus_a = r_ij_tag - a_Si;//r'ij-a
	double r_ij_tag_minus_a2 = r_ij_tag_minus_a*r_ij_tag_minus_a;

	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a);//(r'ij-a)^(-1)
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a2));//(r'ij-a)^(-2)

	double r_ij_tag2 = r_ij_tag*r_ij_tag;
	double r_ij_tag4 = r_ij_tag2*r_ij_tag2;
	double r_ij_tag5 = r_ij_tag4*r_ij_tag;

	double r_ij_tag_in_mFive = (1.0/(r_ij_tag5));//r'ij^(-5)
	double r_ij_tag_in_mFour = (1.0/(r_ij_tag4));//r'ij^(-4)
	double expression = B_Si * r_ij_tag_in_mFour;
	expression = expression - 1.0;//(B*r'ij^(-4) - 1)

	double exponent = exp(r_ij_tag_minus_a_in_mOne);

	double f2_derivative_part_1 = -4.0 * B_Si * r_ij_tag_in_mFive;
	double f2_derivative_part_2 = expression * r_ij_tag_minus_a_in_mTwo;

	return A_Si*exponent*(f2_derivative_part_1 - f2_derivative_part_2);
}

__host__ __device__ double v2_derivative_of_rix(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative * (epsilon_Si/sigma_Si);//v2 derivative of distance
	double dist_x = (i.x - j.x);
	dist_x = dist_x / (r_ij);

	double v2_derivative = f2_derivative * dist_x;

	return v2_derivative;
}

__host__ __device__ double v2_derivative_of_riy(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative * (epsilon_Si/sigma_Si);
	double dist_y = i.y - j.y;
	dist_y = dist_y / (r_ij);

	double v2_derivative = f2_derivative * dist_y;

	return v2_derivative;
}

__host__ __device__ double v2_derivative_of_riz(real3 i, real3 j, double r_ij)
{
	if(r_ij/sigma_Si == pow(2.0,1.0/6.0))
	{
		return 0;
	}

	double f2_derivative = f2_derivative_of_rij_tag(r_ij/sigma_Si);

	f2_derivative = f2_derivative  * (epsilon_Si/sigma_Si);
	double dist_z = i.z - j.z;
	dist_z = dist_z / (r_ij);

	double v2_derivative = f2_derivative * dist_z;

	return v2_derivative;
}
//----------------------------------------------------------------------------//

//-------------------- potential between two Si particles ---------------------//
__host__ __device__ double f2(double r_ij_tag)
{
	if(r_ij_tag >= a_Si)
	{
		return 0;
	}
	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);

	double exponent = exp(r_ij_tag_minus_a_in_mOne);
	double r_ij_tag2 = r_ij_tag*r_ij_tag;
	double r_ij_tag4 = r_ij_tag2*r_ij_tag2;

	double expression = (1.0/(r_ij_tag4));
	expression *= B_Si;
	expression -= 1.0;

	return A_Si*expression*exponent;
}

__host__ __device__ double v2(double r_ij_tag)
{
	if(r_ij_tag == pow(2.0,1.0/6.0))
	{
		return -epsilon_Si;
	}
	return f2(r_ij_tag)*epsilon_Si;
}
//----------------------------------------------------------------------------//

//------------------------ force between three Si particles -------------------//
__host__ __device__ double hi_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a = r_ij_tag - a_Si;
	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a) * gama_Si;
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a*r_ij_tag_minus_a)) * gama_Si;

	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag - a_Si) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_ik_tag_minus_a_in_mOne);

	double expression = (r_ij_tag*r_ij_tag - r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag) / (r_ij_tag*r_ij_tag * r_ik_tag);

	expression -= (r_ij_tag_minus_a_in_mTwo*cosJik_plus_oneThird);

	return lamda_Si*exponent*cosJik_plus_oneThird*expression;
}

__host__ __device__ double hi_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = (1.0/(r_ij_tag - a_Si)) * gama_Si;

	double r_ik_tag_minus_a = r_ik_tag - a_Si;
	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a) * gama_Si;
	double r_ik_tag_minus_a_in_mTwo = (1.0/(r_ik_tag_minus_a*r_ik_tag_minus_a)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_ik_tag_minus_a_in_mOne);

	double expression = (r_ik_tag*r_ik_tag - r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag) / (r_ik_tag*r_ik_tag * r_ij_tag);

	expression -= (r_ik_tag_minus_a_in_mTwo*cosJik_plus_oneThird);

	return lamda_Si*exponent*cosJik_plus_oneThird*expression;
}

__host__ __device__ double hj_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a = r_ij_tag - a_Si;
	double r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a) * gama_Si;
	double r_ij_tag_minus_a_in_mTwo = (1.0/(r_ij_tag_minus_a*r_ij_tag_minus_a)) * gama_Si;

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (r_ij_tag*r_ij_tag - r_jk_tag*r_jk_tag + r_ik_tag*r_ik_tag) / (r_ij_tag*r_ij_tag * r_jk_tag);

	expression -= (r_ij_tag_minus_a_in_mTwo*cosIjk_plus_oneThird);

	return lamda_Si*exponent*cosIjk_plus_oneThird*expression;
}

__host__ __device__ double hj_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = (1.0/(r_ij_tag - a_Si)) * gama_Si;
	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (-r_ik_tag) / (r_ij_tag * r_jk_tag);

	return lamda_Si*exponent*2.0*cosIjk_plus_oneThird*expression;
}

__host__ __device__ double hk_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_ik_tag_minus_a_in_mOne = (1.0/(r_ik_tag - a_Si)) * gama_Si;

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (-r_ij_tag) / (r_ik_tag * r_jk_tag);

	return lamda_Si*exponent*2.0*cosIkj_plus_oneThird*expression;
}

__host__ __device__ double hk_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2.0 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_jk_tag_minus_a_in_mOne = (1.0/(r_jk_tag - a_Si)) * gama_Si;

	double r_ik_tag_minus_a = r_ik_tag - a_Si;
	double r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a) * gama_Si;
	double r_ik_tag_minus_a_in_mTwo = (1.0/(r_ik_tag_minus_a*r_ik_tag_minus_a)) * gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne+r_jk_tag_minus_a_in_mOne);

	double expression = (r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag + r_ij_tag*r_ij_tag) / (r_ik_tag*r_ik_tag * r_jk_tag);

	expression -= (r_ik_tag_minus_a_in_mTwo*cosIkj_plus_oneThird);

	return lamda_Si*exponent*cosIkj_plus_oneThird*expression;
}

__host__ __device__ double f3_derivative_of_rij_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double hi_derivative_of_rij = 0.0;
	double hj_derivative_of_rij = 0.0;
	double hk_derivative_of_rij = 0.0;
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		hi_derivative_of_rij = hi_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ij_tag < a_Si)
	{
		hj_derivative_of_rij = hj_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		hk_derivative_of_rij = hk_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}

	return hi_derivative_of_rij + hj_derivative_of_rij + hk_derivative_of_rij;
}

__host__ __device__ double f3_derivative_of_rik_tag(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double hi_derivative_of_rik = 0.0;
	double hj_derivative_of_rik = 0.0;
	double hk_derivative_of_rik = 0.0;
	if(r_ik_tag < a_Si && r_ij_tag < a_Si)
	{
		hi_derivative_of_rik = hi_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ij_tag < a_Si)
	{
		hj_derivative_of_rik = hj_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		hk_derivative_of_rik = hk_derivative_of_rij_tag(r_ij_tag,r_ik_tag,r_jk_tag);
	}

	return hi_derivative_of_rik + hj_derivative_of_rik + hk_derivative_of_rik;
}

__host__ __device__ double v3_derivative_of_rix(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijx = (i.x-j.x);
	double dist_ikx = (i.x-k.x);
	double expression1 = (dist_ijx/(r_ij));
	double expression2 = (dist_ikx/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}

__host__ __device__ double v3_derivative_of_riy(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijy = (i.y-j.y);
	double dist_iky = (i.y-k.y);
	double expression1 = (dist_ijy/(r_ij));
	double expression2 = (dist_iky/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}

__host__ __device__ double v3_derivative_of_riz(real3 i, real3 j, real3 k, double r_ij, double r_ik, double r_jk)
{
	double v3_derived_by_rij = (f3_derivative_of_rij_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);
	double v3_derived_by_rik = (f3_derivative_of_rik_tag(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si))*(epsilon_Si/sigma_Si);

	double dist_ijz = (i.z-j.z);
	double dist_ikz = (i.z-k.z);
	double expression1 = (dist_ijz/(r_ij));
	double expression2 = (dist_ikz/(r_ik));

	return v3_derived_by_rij*expression1 + v3_derived_by_rik*expression2;
}
//----------------------------------------------------------------------------//


//-------------------- potential between three Si particles -------------------//
__host__ __device__ double hi(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosJik_plus_oneThird = ((r_ij_tag*r_ij_tag + r_ik_tag*r_ik_tag - r_jk_tag*r_jk_tag)/(2.0 * r_ij_tag * r_ik_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);
	r_ij_tag_minus_a_in_mOne *= gama_Si;

	double r_ik_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a_in_mOne);
	r_ik_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne + r_ik_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosJik_plus_oneThird*cosJik_plus_oneThird;
}

__host__ __device__ double hj(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIjk_plus_oneThird = ((r_ij_tag*r_ij_tag + r_jk_tag*r_jk_tag - r_ik_tag*r_ik_tag)/(2.0 * r_ij_tag * r_jk_tag)) + (1.0/3.0);

	double r_ij_tag_minus_a_in_mOne = r_ij_tag - a_Si;
	r_ij_tag_minus_a_in_mOne = (1.0/r_ij_tag_minus_a_in_mOne);
	r_ij_tag_minus_a_in_mOne *= gama_Si;

	double r_jk_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_jk_tag_minus_a_in_mOne = (1/r_jk_tag_minus_a_in_mOne);
	r_jk_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ij_tag_minus_a_in_mOne + r_jk_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosIjk_plus_oneThird*cosIjk_plus_oneThird;
}

__host__ __device__ double hk(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double cosIkj_plus_oneThird = ((r_ik_tag*r_ik_tag + r_jk_tag*r_jk_tag - r_ij_tag*r_ij_tag)/(2.0 * r_ik_tag * r_jk_tag)) + (1.0/3.0);

	double r_ik_tag_minus_a_in_mOne = r_ik_tag - a_Si;
	r_ik_tag_minus_a_in_mOne = (1.0/r_ik_tag_minus_a_in_mOne);
	r_ik_tag_minus_a_in_mOne *= gama_Si;

	double r_jk_tag_minus_a_in_mOne = r_jk_tag - a_Si;
	r_jk_tag_minus_a_in_mOne = (1.0/r_jk_tag_minus_a_in_mOne);
	r_jk_tag_minus_a_in_mOne *= gama_Si;

	double exponent = exp(r_ik_tag_minus_a_in_mOne + r_jk_tag_minus_a_in_mOne);

	return lamda_Si*exponent*cosIkj_plus_oneThird*cosIkj_plus_oneThird;
}

__host__ __device__ double f3(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	double h_i = 0.0;
	double h_j = 0.0;
	double h_k = 0.0;
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		h_i = hi(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_ij_tag < a_Si && r_ik_tag < a_Si)
	{
		h_j = hj(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	if(r_jk_tag < a_Si && r_ik_tag < a_Si)
	{
		h_k = hk(r_ij_tag,r_ik_tag,r_jk_tag);
	}
	return h_i + h_j + h_k;
}

__host__ __device__ double v3(double r_ij_tag, double r_ik_tag, double r_jk_tag)
{
	return f3(r_ij_tag,r_ik_tag,r_jk_tag)*epsilon_Si;
}
//----------------------------------------------------------------------------//
