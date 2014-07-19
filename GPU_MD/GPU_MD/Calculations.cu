/*
 *
 */

#include "Calculations.h"

__host__ __device__ double distance2(real3 i, real3 j)
{
	return sqrt((i.x-j.x)*(i.x-j.x) + (i.y-j.y)*(i.y-j.y) + (i.z-j.z)*(i.z-j.z));
}

//-------------------------- calculate Force Si ------------------------------//
__global__ void d_calculateForce_Si(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, int* a)
{
	int idx = threadIdx.x + blockIdx.x*blockDim.x;
	a[idx] *= 3;
/*	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		particleStruct particle = siParticles[i];

		real3 iPosition = particle.position;
		int countSi = 0;
		int countXe = 0;
		double r_ij = 0.0;
		double r_ik = 0.0;
		double r_jk = 0.0;
		siParticles[i].force.x = 0.0;
		siParticles[i].force.y = 0.0;
		siParticles[i].force.z = 0.0;
		real3 jPosition;
		real3 kPosition;

		if(config.USE_NEIGHBOR_LISTS)
		{
			for(countSi = 0; countSi < MAX_SI_NEIGHBORS && particle.siNeighbors[countSi] != -1; countSi++);
			for(countXe = 0; countXe < MAX_XE_NEIGHBORS && particle.xeNeighbors[countXe] != -1; countXe++);
		}

		if(!config.USE_NEIGHBOR_LISTS)
		{
			countSi = config.SI_PARTICLES;
			countXe = config.XE_PARTICLES;
		}

		if(config.useLennardJonesPotentialForSi)
		{
			for(int j = 0; j < countSi; j++)
			{
				if(j != i)
				{
					if(!config.USE_NEIGHBOR_LISTS)
						jPosition = siParticles[j].position;
					else
						jPosition = siParticles[particle.siNeighbors[j]].position;

					r_ij = distance2(iPosition, jPosition);	
		//			if(r_ij/sigma_Si < 1.8)
					{
						siParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
						siParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
						siParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
					}
				}
			}
		}

		else
		{	
			for(int j = 0; j < countSi; j++)
			{
				if(j != i)
				{
					if(!config.USE_NEIGHBOR_LISTS)
						jPosition = siParticles[j].position;
					else
						jPosition = siParticles[particle.siNeighbors[j]].position;
			
					r_ij = distance2(iPosition, jPosition);

					if(r_ij/sigma_Si < a_Si)
					{
						siParticles[i].force.x -= v2_derivative_of_rix(iPosition, jPosition, r_ij);
						siParticles[i].force.y -= v2_derivative_of_riy(iPosition, jPosition, r_ij);
						siParticles[i].force.z -= v2_derivative_of_riz(iPosition, jPosition, r_ij);
					}
					for(int k = 0; k < countSi; k++)
					{
						if(k != i && k != j)
						{
							if(!config.USE_NEIGHBOR_LISTS)
							{
								kPosition = siParticles[k].position;							
							}
							else
							{
								kPosition = siParticles[particle.siNeighbors[k]].position;
							}

							r_ik = distance2(iPosition, kPosition);
							r_jk = distance2(jPosition, kPosition);
							if((r_ij/sigma_Si < a_Si && r_ik/sigma_Si < a_Si) || (r_ij/sigma_Si < a_Si && r_jk/sigma_Si < a_Si) || (r_ik/sigma_Si < a_Si && r_jk/sigma_Si < a_Si))
							{
								siParticles[i].force.x -= v3_derivative_of_rix(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
								siParticles[i].force.x -= v3_derivative_of_rix(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
								siParticles[i].force.y -= v3_derivative_of_riy(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
								siParticles[i].force.y -= v3_derivative_of_riy(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
								siParticles[i].force.z -= v3_derivative_of_riz(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
								siParticles[i].force.z -= v3_derivative_of_riz(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
							}
						}
					}
				}
			}
		}*/

	/*	for(int j = 0; j < countXe; j++)
		{
				if(!config.USE_NEIGHBOR_LISTS)
					jPosition = xeParticles[j].position;
				else
					jPosition = xeParticles[particle.xeNeighbors[j]];

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

__global__ void d_calculateForce_Xe(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, particleStruct* xeParticles, particleStruct* siParticles, configurations config)
{
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		particleStruct particle = xeParticles[i];

		real3 iPosition = particle.position;

		int countSi = 0;
		int countXe = 0;
		double r_ij = 0.0;
		xeParticles[i].force.x = 0.0;
		xeParticles[i].force.y = 0.0;
		xeParticles[i].force.z = 0.0;
		real3 jPosition;

		if(config.USE_NEIGHBOR_LISTS)
		{
			for(countSi = 0; countSi < MAX_SI_NEIGHBORS && particle.siNeighbors[countSi] != -1; countSi++);
			for(countXe = 0; countXe < MAX_XE_NEIGHBORS && particle.xeNeighbors[countXe] != -1; countXe++);
		}

		if(!config.USE_NEIGHBOR_LISTS)
		{
			countSi = config.SI_PARTICLES;
			countXe = config.XE_PARTICLES;
		}

		for(int j = 0; j < countXe; j++)
		{
			if(j != i)
			{
				if(!config.USE_NEIGHBOR_LISTS)
					jPosition = xeParticles[j].position;
				else
					jPosition = xeParticles[particle.xeNeighbors[j]].position;

				r_ij = distance2(iPosition, jPosition);
				if(r_ij/sigma_Xe_Xe < xe_Cluster)
				{
					xeParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
					xeParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
					xeParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
				}
			}
		}

	/*	for(int j = 0; j < countSi; j++)
		{
			if(!config.USE_NEIGHBOR_LISTS)
				jPosition = siParticles[j].position;
			else
				jPosition = siParticles[particle.siNeighbors[j]];

			r_ij = distance2(iPosition, jPosition);
			if(r_ij/sigma_Si_Xe < a_Si_Xe)
			{
				xeParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
				xeParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
				xeParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
			}
		}*/
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////!!!!
//----------------------calculate total potential of Si -----------------------//
/*__global__ double d_V_total_Si(int SI_PARTICLES, particleStruct* siParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = i+1; j < SI_PARTICLES; j++)
		{
			vTotal += v2(distance2(siParticles[i].position, siParticles[j].position)/sigma_Si);
		}
	}
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = i+1; j < SI_PARTICLES; j++)
		{
			for(int k = j+1; k < SI_PARTICLES; k++)
			{
				vTotal += v3(distance2(siParticles[i].position,siParticles[j].position)/sigma_Si, distance2(siParticles[i].position,siParticles[k].position)/sigma_Si, distance2(siParticles[j].position,siParticles[k].position)/sigma_Si);
			}
		}
	}
	for(int i = 0; i < SI_PARTICLES; i++)
	{
//		cout<<(pow((sqrt((siParticles[i].velocity.x*siParticles[i].velocity.x)+(siParticles[i].velocity.y*siParticles[i].velocity.y)+(siParticles[i].velocity.z*siParticles[i].velocity.z))),2.0)*SiMass)/2<<endl;
//		vTotal += (pow((sqrt((siParticles[i].velocity.x*siParticles[i].velocity.x)+(siParticles[i].velocity.y*siParticles[i].velocity.y)+(siParticles[i].velocity.z*siParticles[i].velocity.z))),2.0)*SiMass)/2;
	}
	return vTotal;
}
//----------------------------------------------------------------------------//
__global__ double d_V_total_Si_Xe(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = 0; j < XE_PARTICLES; j++)
		{
//			if(distance2(siParticles[i].position, xeParticles[j].position)/sigma_Si_Xe < a_Si_Xe)
				vTotal += lennardJonesPotential(distance2(siParticles[i].position, xeParticles[j].position), sigma_Si_Xe, epsilon_Si_Xe);
		}
	}
	return vTotal;
}

__global__ double d_V_total_Xe(int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < XE_PARTICLES-1; i++)
	{
		for(int j = i+1; j < XE_PARTICLES; j++)
		{
			if(distance2(xeParticles[i].position, xeParticles[j].position)/sigma_Xe_Xe < a_Si_Xe)
				vTotal += lennardJonesPotential(distance2(xeParticles[i].position, xeParticles[j].position), sigma_Xe_Xe, epsilon_Xe_Xe);
		}
	}
	for(int i = 0; i < XE_PARTICLES; i++)
	{
		vTotal += (pow((sqrt((xeParticles[i].velocity.x*xeParticles[i].velocity.x)+(xeParticles[i].velocity.y*xeParticles[i].velocity.y)+(xeParticles[i].velocity.z*xeParticles[i].velocity.z))),2.0)*XeMass)/2;
	}
	return vTotal;
}

__global__ void d_V_total(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	vTotal += V_total_Si(SI_PARTICLES, siParticles);
	vTotal += V_total_Si_Xe(SI_PARTICLES, siParticles, XE_PARTICLES, xeParticles);
	vTotal += V_total_Xe(XE_PARTICLES, xeParticles);
	cout<<vTotal<<endl;
}*/
//////////////////////////////////////////////////////////////////////////////////////////////!!!!
__global__ void d_initiateAcceleration(particleStruct *particles, int listSize, double mass)
{
	for(int i = 0; i < listSize; i++)
	{
		particles[i].aAcc.x = particles[i].force.x/mass;
		particles[i].aAcc.y = particles[i].force.y/mass;
		particles[i].aAcc.z = particles[i].force.z/mass;
	}
}

__global__ void d_predict(particleStruct *particles, int listSize, float dt)
{
	double c1 = dt;
    double c2 = c1*dt/2.0;
    double c3 = c2*dt/3.0;
    double c4 = c3*dt/4.0;

	for(int i = 0; i < listSize; i++)
	{
		particles[i].position.x += c1*particles[i].velocity.x + c2*particles[i].aAcc.x + c3*particles[i].bAcc.x + c4*particles[i].cAcc.x;
		particles[i].position.y += c1*particles[i].velocity.y + c2*particles[i].aAcc.y + c3*particles[i].bAcc.y + c4*particles[i].cAcc.y;
		particles[i].position.z += c1*particles[i].velocity.z + c2*particles[i].aAcc.z + c3*particles[i].bAcc.z + c4*particles[i].cAcc.z;
		particles[i].velocity.x += c1*particles[i].aAcc.x + c2*particles[i].bAcc.x + c3*particles[i].cAcc.x;
		particles[i].velocity.y += c1*particles[i].aAcc.y + c2*particles[i].bAcc.y + c3*particles[i].cAcc.y;
		particles[i].velocity.z += c1*particles[i].aAcc.z + c2*particles[i].bAcc.z + c3*particles[i].cAcc.z;
		particles[i].aAcc.x += c1*particles[i].bAcc.x + c2*particles[i].cAcc.x;
		particles[i].aAcc.y += c1*particles[i].bAcc.y + c2*particles[i].cAcc.y;
		particles[i].aAcc.z += c1*particles[i].bAcc.z + c2*particles[i].cAcc.z;
		particles[i].bAcc.x += c1*particles[i].cAcc.x;
		particles[i].bAcc.y += c1*particles[i].cAcc.y;
		particles[i].bAcc.z += c1*particles[i].cAcc.z;
	}
}

__global__ void d_correct(particleStruct *particles, double dt, int listSize, double mass)
{
	double c1 = dt ;
	double c2 = c1*dt/2.0;
	double c3 = c2*dt/3.0; 
	double c4 = c3*dt/4.0;

	double cr = GEAR1*c2;
	double cv = GEAR2*c2/c1;
	double cb = GEAR3*c2/c3;
	double cc = GEAR4*c2/c4;

	for(int i = 0; i < listSize; i++)
	{
		double axi = particles[i].force.x/mass;
      	double ayi = particles[i].force.y/mass;
      	double azi = particles[i].force.z/mass;

		double corrx = axi - particles[i].aAcc.x;
    	double corry = ayi - particles[i].aAcc.y;
    	double corrz = azi - particles[i].aAcc.z;

		particles[i].position.x += cr*corrx;
		particles[i].position.y += cr*corry;
		particles[i].position.z += cr*corrz;
		particles[i].velocity.x += cv*corrx;
		particles[i].velocity.y += cv*corry;
		particles[i].velocity.z += cv*corrz;
		particles[i].aAcc.x = axi;
		particles[i].aAcc.y = ayi;
		particles[i].aAcc.z = azi;
		particles[i].bAcc.x += cb*corrx;
		particles[i].bAcc.y += cb*corry;
		particles[i].bAcc.z += cb*corrz;
		particles[i].cAcc.x += cc*corrx;
		particles[i].cAcc.y += cc*corry;
    	particles[i].cAcc.z += cc*corrz;
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
	double s = 2.0*invdist14-invdist8;
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
	double first = -4*B_Si*(1.0/(r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag));
	double second = ((B_Si*(1.0/(r_ij_tag*r_ij_tag*r_ij_tag*r_ij_tag)))-1)*(1.0/((r_ij_tag-a_Si)*(r_ij_tag-a_Si)));
	double print = first-second;

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
