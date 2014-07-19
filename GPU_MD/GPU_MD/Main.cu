/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#include "Configurations.h"
#include "InitiateSystem.h"
#include "Calculations.h"
#include "Calculations.cuh"
#include "NeighborLists.h"
#include <time.h>
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Equi.h"
#include <math.h>
#include <glut.h> 
#include "DebugPrints.h"
#include <cuda.h>//cuda
#include <cuda_runtime.h>//cuda
#include <device_launch_parameters.h>//cuda
#include <ctime>

GLfloat *xp;
GLfloat *yp;
GLfloat *zp;
configurations config;
allLists lists;
ofstream file;

int place;
int resolution = 150;
int iterationNum;

time_t startT, endT;

double prevVPsi = 0.0;
double prevVKsi = 0.0;
double prevVPxe = 0.0;
double prevVKxe = 0.0;

bool debugSi;
bool debugXe;
bool chooseVeloSi;
bool chooseVeloXe;
bool useNeighborList;
bool animation;
bool sysPause;
bool useLennardJonesPotentialForSi;

void Animate(int,char**);
void Initialize();
void Draw();
void DrawPoints();
void calculatePositions();

ofstream si_place_and_velocity; 
ofstream xe_place_and_velocity; 
ofstream si_potential; 
ofstream xe_potential;

unsigned t0;
unsigned elapsed;
double averageTime = 0.0;

int main(int argc, char *argv[])
{	
	cout.precision(19);
	si_place_and_velocity.open("si_place_and_velocity");
	xe_place_and_velocity.open("xe_place_and_velocity");
	si_potential.open("si_potential"); 
	xe_potential.open("xe_potential"); 
	
	if(!si_place_and_velocity.good() || !xe_place_and_velocity.good() || !si_potential.good() || !xe_potential.good())
	{
		cout<<"Can not open the files for debug"<<endl;
		return EXIT_FAILURE;
	}

	si_place_and_velocity.precision(19);
	xe_place_and_velocity.precision(19);
	si_potential.precision(19);
	xe_potential.precision(19);
	iterationNum = 0;

	cout<<"read configurations\n";
	readConfig(argv[1], &config);
	cout<<"read configurations\n";

	///// debug outputs //////
	if(config.PRINT_GRAPHS)
	{
		printDebug(config.OUT_FOR_GRAPHS);
	}
	//////////////////////////	

	////////////////////// DEBUG CONFIGURATIONS ////////////////////
	debugSi = config.DEBUG_SI;
	debugXe = config.DEBUG_XE;
	chooseVeloSi = config.CHOOSE_VELO_SI;
	chooseVeloXe = config.CHOOSE_VELO_XE;
	useNeighborList = config.USE_NEIGHBOR_LISTS;
	animation = config.ANIMATE;
	sysPause = config.SYS_PAUSE;
	useLennardJonesPotentialForSi = config.useLennardJonesPotentialForSi;
	/////////////////////////////////////////////////////////////////

	cout<<"write input file\n";
	writeParticlesInput(&config);

	cout<<"lists allocation\n";
	if(listsAlloc(&config, &lists) == EXIT_FAILURE)
	{
		system("pause");
		exit(EXIT_FAILURE);
	}

	cout<<"read input\n";
	readInput(*(config.INPUT_FILE), &lists, &config);

	if(chooseVeloSi)
	{
		cout<<"choose Si velocities\n";
		if(config.SI_PARTICLES > 0 && chooseVeloSi)
			chooseVelocities(lists.siParticles, config.SI_PARTICLES, config.TEMPERATURE, &config, SiMass);
	}

	if(chooseVeloXe)
	{
		cout<<"choose Si velocities\n";
		if(config.XE_PARTICLES > 0 && chooseVeloXe)
			chooseVelocities(lists.xeParticles, config.XE_PARTICLES, config.TEMPERATURE, &config, XeMass);
	}
	
	if(useNeighborList)
	{
		cout<<"build neighbors\n";
		divideToCells(lists.siParticles, lists.xeParticles, config.SI_PARTICLES, config.XE_PARTICLES, si_xe_Cluster);

		buildNeighbors(lists.siParticles, lists.xeParticles, &config, si_xe_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.SI_HEIGHT+config.XE_HEIGHT+config.LA_SPACE)*2, XeSi, config.SI_PARTICLES, config.XE_PARTICLES);
		buildNeighbors(lists.xeParticles, lists.siParticles, &config, si_xe_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.SI_HEIGHT+config.XE_HEIGHT+config.LA_SPACE)*2, XeSi, config.XE_PARTICLES, config.SI_PARTICLES);

		divideToCells(lists.siParticles, config.SI_PARTICLES, si_Cluster);
		buildNeighbors(lists.siParticles, NULL, &config, si_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.XE_HEIGHT+config.SI_HEIGHT+config.LA_SPACE)*2, Si, config.SI_PARTICLES, 0);

		divideToCells(lists.xeParticles, config.XE_PARTICLES, xe_Cluster);
		buildNeighbors(lists.xeParticles, NULL, &config, xe_Cluster,  max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.XE_HEIGHT+config.SI_HEIGHT+config.LA_SPACE)*2, Xe, config.XE_PARTICLES, 0);
	}
	
	cout<<"calculate forces for first time\n";
//	calculateForce_Si(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, lists.siParticles, lists.xeParticles, config);//!
//	calculateForce_Xe(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, lists.xeParticles, lists.siParticles, config);//!
	
//////////////////////////////////////////////////CUDA
	particleStruct* d_siParticles;
	particleStruct* d_xeParticles;
	int* a = new int[5];
	a[0] = 0;
	a[1] = 1;
	a[2] = 2;
	a[3] = 3;
	a[4] = 4;
	int* d_a;
	cudaMalloc(&d_a, sizeof(int)*5);
	cudaMemcpy(d_a, a, sizeof(int)*5, cudaMemcpyHostToDevice);
	cudaMalloc(&d_siParticles, sizeof(particleStruct)*config.SI_PARTICLES);
	cudaMalloc(&d_xeParticles, sizeof(particleStruct)*config.XE_PARTICLES);
	cudaMemcpy(d_siParticles, lists.siParticles, sizeof(particleStruct)*config.SI_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy(d_xeParticles, lists.xeParticles, sizeof(particleStruct)*config.XE_PARTICLES, cudaMemcpyHostToDevice);
	d_calculateForce_Si<<<1,1024>>>(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, d_a);//!
	// wait until tasks are completed
	cudaDeviceSynchronize();
	// check for errors
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) 
	{
	  fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
	}
	cudaMemcpy(a, d_a, sizeof(int)*5, cudaMemcpyDeviceToHost);
	cout<<a[0]<<endl;
	cout<<a[1]<<endl;
	cout<<a[2]<<endl;
	cout<<a[3]<<endl;
	cout<<a[4]<<endl;

	d_calculateForce_Xe<<<100,1024>>>(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, d_xeParticles, d_siParticles, config);//!
	cudaDeviceSynchronize();
	// check for errors
	error = cudaGetLastError();
	if (error != cudaSuccess) 
	{
	  fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
	}
	
	cudaMemcpy(lists.siParticles, d_siParticles, sizeof(particleStruct)*config.SI_PARTICLES, cudaMemcpyDeviceToHost);
	cudaMemcpy( lists.xeParticles, d_xeParticles,sizeof(particleStruct)*config.XE_PARTICLES, cudaMemcpyDeviceToHost);
//////////////////////////////////////////////////CUDA

	if(config.SI_PARTICLES > 0)
		initiateAcceleration(lists.siParticles, config.SI_PARTICLES, SiMass);
	if(config.XE_PARTICLES > 0)
		initiateAcceleration(lists.xeParticles, config.XE_PARTICLES, XeMass);

	if(animation)
	{
		time(&startT);
		Animate(argc, argv);
	}

	else
	{
		for(int i = 0; i < config.STEPS; i++)
		{
			calculatePositions();
		}
		si_place_and_velocity.close();
		xe_place_and_velocity.close();
		si_potential.close();
		xe_potential.close();
	}
//	system("pause");
	cudaDeviceReset();
	cout<<averageTime/(iterationNum)<<endl;
	return EXIT_SUCCESS;
}

void calculatePositions()
{
	cout<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
	if(debugSi)
	{
//		cout<<"SI: "<<endl<<endl;
		double vP = 0.0;
		double vK = 0.0;
		double V3 = 0.0;
		real3 iPosition;
		real3 jPosition;
		real3 kPosition;

		si_place_and_velocity<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		for(int i = 0; i < config.SI_PARTICLES; i++)
		{
			double potentialForOneAtom = 0.0;
			double potential = 0.0;
			iPosition = lists.siParticles[i].position;
			//print place and velocity
			si_place_and_velocity<<"\n"<<i<<" particle\n"<<"place: "<<iPosition.x<<"  "<<iPosition.y<<"  "<<iPosition.z<<"\n"<<
				"velocity: "<<lists.siParticles[i].velocity.x<<"  "<<lists.siParticles[i].velocity.y<<"  "<<lists.siParticles[i].velocity.z<<"\n"<<
				"k energy: "<<(((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z))*SiMass)/2.0<<"\n";

			for(int j = 0; j < config.SI_PARTICLES; j++)
			{
				if(i != j)
				{
					jPosition = lists.siParticles[j].position;
					double r_ij = distance2(iPosition,jPosition);
					if(useLennardJonesPotentialForSi)
					{
						potential = lennardJonesPotential(r_ij,sigma_Si,epsilon_Si);
						potentialForOneAtom += potential;
						if(i < j)
						{
							vP += potential;
						}
					}
					else
					{
						potential = v2(r_ij/sigma_Si);
						potentialForOneAtom += potential;
						vP += potential;
						for(int k = 0; k < config.SI_PARTICLES; k++)
						{
							if(k != i && k != j)
							{
								kPosition = lists.siParticles[k].position;
								double r_ik = distance2(iPosition,kPosition);
								double r_jk = distance2(jPosition,kPosition);
								if((r_ij/sigma_Si < a_Si && r_ik/sigma_Si < a_Si) || (r_ij/sigma_Si < a_Si && r_jk/sigma_Si < a_Si) || (r_ik/sigma_Si < a_Si && r_jk/sigma_Si < a_Si))
								{
									potential = v3(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si);
									potentialForOneAtom += potential;
									if(i < j && j < k)
									{
										vP += potential;
									}
								}
							}
						}
					}
				}
			}
			si_place_and_velocity<<"potential: "<<potentialForOneAtom<<"\n";
		}
		for(int i = 0; i < config.SI_PARTICLES; i++)
		{
			vK += (((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z))*SiMass)/2.0;
		}

		vP /= 2;
		si_potential<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		si_potential<<endl<<"Ep	   	    ";
		si_potential<<"Ek		            ";
		si_potential<<"E_total"<<"\n\n";
		si_potential<<endl<<vP<<" "; // Ep
		si_potential<<vK<<" "; // Ek
		si_potential<<vP+vK<<"\n\n"; // E_total = Ep + Ek
		si_potential<<endl<<"Change in Ep    ";
		si_potential<<"Change in Ek   ";
		si_potential<<"Change in E_total"<<"\n\n";
		si_potential<<endl<<(vP-prevVPsi)<<" "; // change in Ep
		si_potential<<(vK-prevVKsi)<<" "; // change in Ek
		si_potential<<((vP-prevVPsi)+(vK-prevVKsi))<<"\n\n"; // change in E_total
		// determine which energy changed more (non if equals)
		if(abs(vK-prevVKsi) < abs(vP-prevVPsi))
			si_potential<<"vP"<<endl<<endl;
		else if(abs(vK-prevVKsi) > abs(vP-prevVPsi))
			si_potential<<"vK"<<endl<<endl;
		else
			si_potential<<"non"<<endl<<endl;
		prevVPsi = vP;
		prevVKsi = vK;
		si_place_and_velocity<<"\n\n\n\n";
	}
	/////////////////////////////////////////////////////////////////////////////////

	if(debugXe)
	{
//		cout<<"XE: "<<endl<<endl;
		double vP = 0.0;
		double vK = 0.0;
		real3 iPosition;
		real3 jPosition;

		xe_place_and_velocity<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		for(int i = 0; i < config.XE_PARTICLES; i++)
		{
			double potentialForOneAtom = 0.0;
			double potential = 0.0;
			iPosition = lists.xeParticles[i].position;
			//print place and velocity
			xe_place_and_velocity<<"\n"<<i<<" particle\n"<<"place: "<<iPosition.x<<"  "<<iPosition.y<<"  "<<iPosition.z<<"\n"<<
				"velocity: "<<lists.xeParticles[i].velocity.x<<"  "<<lists.xeParticles[i].velocity.y<<"  "<<lists.xeParticles[i].velocity.z<<"\n"<<
				"k energy: "<<(((lists.xeParticles[i].velocity.x*lists.xeParticles[i].velocity.x)+(lists.xeParticles[i].velocity.y*lists.xeParticles[i].velocity.y)+(lists.xeParticles[i].velocity.z*lists.xeParticles[i].velocity.z))*XeMass)/2.0<<"\n";

			for(int j = 0; j < config.XE_PARTICLES; j++)
			{
				if(i != j)
				{
					jPosition = lists.xeParticles[j].position;
					potential = lennardJonesPotential(distance2(iPosition, jPosition),sigma_Xe_Xe,epsilon_Xe_Xe);
					potentialForOneAtom += potential;
					if(i < j)
					{
						vP += potential;
					}
				}
			}
			xe_place_and_velocity<<"potential: "<<potentialForOneAtom<<"\n";
		}
		for(int i = 0; i < config.XE_PARTICLES; i++)
		{
			vK += (((lists.xeParticles[i].velocity.x*lists.xeParticles[i].velocity.x)+(lists.xeParticles[i].velocity.y*lists.xeParticles[i].velocity.y)+(lists.xeParticles[i].velocity.z*lists.xeParticles[i].velocity.z))*XeMass)/2.0;
		}

		vP /= 2;
		xe_potential<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		xe_potential<<endl<<"Ep	   	    ";
		xe_potential<<"Ek		            ";
		xe_potential<<"E_total"<<"\n\n";
		xe_potential<<endl<<vP<<" "; // Ep
		xe_potential<<vK<<" "; // Ek
		xe_potential<<vP+vK<<"\n\n"; // E_total = Ep + Ek
		xe_potential<<endl<<(vP-prevVPxe)<<" "; // change in Ep
		xe_potential<<(vK-prevVKxe)<<" "; // change in Ek
		xe_potential<<((vP-prevVPxe)+(vK-prevVKxe))<<"\n\n"; // change in E_total
		// determine which energy changed more (non if equals)
		if(abs(vK-prevVKxe) < abs(vP-prevVPxe))
			xe_potential<<"vP"<<endl<<endl;
		else if(abs(vK-prevVKxe) > abs(vP-prevVPxe))
			xe_potential<<"vK"<<endl<<endl;
		else
			xe_potential<<"non"<<endl<<endl;
		prevVPxe = vP;
		prevVKxe = vK;
		xe_place_and_velocity<<"\n\n\n\n";
	}

	if((debugSi || debugXe) && sysPause)
	{
		system("pause");
	}

	if(config.SI_PARTICLES > 0)
		predict(lists.siParticles, config.SI_PARTICLES, config.TIMESTEPS);
	if(config.XE_PARTICLES > 0)
		predict(lists.xeParticles, config.XE_PARTICLES, config.TIMESTEPS);
	
	t0 = clock();
//	calculateForce_Si(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, lists.siParticles, lists.xeParticles, config);//!
//	calculateForce_Xe(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, lists.xeParticles, lists.siParticles, config);//!

//////////////////////////////////////////////////CUDA
	particleStruct* d_siParticles;
	particleStruct* d_xeParticles;
	int* a = new int[5];
	a[0] = 0;
	a[1] = 1;
	a[2] = 2;
	a[3] = 3;
	a[4] = 4;
	int* d_a;
	cudaMalloc(&d_a, sizeof(int)*5);
	cudaMemcpy(d_a, a, sizeof(int)*5, cudaMemcpyHostToDevice);
	cudaMalloc(&d_siParticles, sizeof(particleStruct)*config.SI_PARTICLES);
	cudaMalloc(&d_xeParticles, sizeof(particleStruct)*config.XE_PARTICLES);
	cudaMemcpy(d_siParticles, lists.siParticles, sizeof(particleStruct)*config.SI_PARTICLES, cudaMemcpyHostToDevice);
	cudaMemcpy(d_xeParticles, lists.xeParticles, sizeof(particleStruct)*config.XE_PARTICLES, cudaMemcpyHostToDevice);
	d_calculateForce_Si<<<1,1024>>>(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, d_a);//!
	// wait until tasks are completed
	cudaDeviceSynchronize();
	// check for errors
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) 
	{
	  fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
	}
	cudaMemcpy(a, d_a, sizeof(int)*5, cudaMemcpyDeviceToHost);
	cout<<a[0]<<endl;
	cout<<a[1]<<endl;
	cout<<a[2]<<endl;
	cout<<a[3]<<endl;
	cout<<a[4]<<endl;

	d_calculateForce_Xe<<<1,1024>>>(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS, d_xeParticles, d_siParticles, config);//!
	cudaDeviceSynchronize();
	// check for errors
	error = cudaGetLastError();
	if (error != cudaSuccess) 
	{
	  fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
	}

	cudaMemcpy(lists.siParticles, d_siParticles, sizeof(particleStruct)*config.SI_PARTICLES, cudaMemcpyDeviceToHost);
	cudaMemcpy( lists.xeParticles, d_xeParticles,sizeof(particleStruct)*config.XE_PARTICLES, cudaMemcpyDeviceToHost);
//////////////////////////////////////////////////CUDA

	elapsed = clock()-t0;
	cout<<elapsed<<endl;
	averageTime += elapsed;

	if(config.SI_PARTICLES > 0)
		correct(lists.siParticles, config.TIMESTEPS, config.SI_PARTICLES, SiMass);
	if(config.XE_PARTICLES > 0)
		correct(lists.xeParticles, config.TIMESTEPS, config.XE_PARTICLES, XeMass);

	iterationNum++;
}

void Animate(int argc, char *argv[])
{
	xp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	yp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	zp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(950, 650);
	glutInitWindowPosition(50, 50);
	glutCreateWindow("Molecular Dynamics");
	Initialize();
	glutDisplayFunc(Draw);
	glutMainLoop();
	free(xp);
	free(yp);
	free(zp);
}

void Initialize()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}

void Draw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	DrawPoints();
	glFlush();

	calculatePositions();

	glutPostRedisplay();
}

void DrawPoints()
{
	for(place = 0; place < config.SI_PARTICLES; place++)
	{
		xp[place] = lists.siParticles[place].position.x;
		yp[place] = lists.siParticles[place].position.y;
		zp[place] = lists.siParticles[place].position.z;
	}
	for(int i = 0; i < config.XE_PARTICLES; i++, place++)
	{
		xp[place] = (lists.xeParticles[i].position.x)+0.25*config.SI_LENGTH*space_Si;
		yp[place] = (lists.xeParticles[i].position.y)+0.25*config.SI_LENGTH*space_Si;
		zp[place] = (lists.xeParticles[i].position.z)+config.SI_HEIGHT+config.LA_SPACE;
	}
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		glColor3f(1.0, 1.0, 1.0);
		glPointSize(2.5);
		glBegin(GL_POINTS);
			glVertex3f(((xp[i]/resolution)+0.15), ((zp[i]/resolution)+0.25), ((yp[i]/resolution))+0.25);
		glEnd();
	}
	for(int i = config.SI_PARTICLES; i < config.SI_PARTICLES+config.XE_PARTICLES; i++)
	{
		glColor3f(0.246, 0.554, 0.51908);
		glPointSize(2.5);
		glBegin(GL_POINTS);
			glVertex3f(((xp[i]/resolution)+0.15), ((zp[i]/resolution)+0.25), ((yp[i]/resolution))+0.25);
		glEnd();
	}
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		glColor3f(1.0, 1.0, 1.0);
		glPointSize(2.5);
		glBegin(GL_POINTS);
			glVertex3f(((xp[i]/resolution)+0.55), ((yp[i]/resolution)+0.25), ((zp[i]/resolution)+0.25));
		glEnd();
	}
	for(int i = config.SI_PARTICLES; i < config.SI_PARTICLES+config.XE_PARTICLES; i++)
	{
		glColor3f(0.246, 0.554, 0.51908);
		glPointSize(2.5);
		glBegin(GL_POINTS);
			glVertex3f(((xp[i]/resolution)+0.55), ((yp[i]/resolution)+0.25), ((zp[i]/resolution)+0.25));
		glEnd();
	}
}
