
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>
using namespace std;
/////////////////////////
#ifndef _MYCODE_H_
#define _MYCODE_H_
#ifdef CUDADLL64_EXPORTS
#define CUDADLL64_API _declspec( dllexport )
#else
#define CUDADLL64_API _declspec(dllimport)
#endif
extern "C" CUDADLL64_API void HH_Entrance(float *population, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA, int POPULATION_SIZE);
#endif
///////////////////////////////////////////////Ca
class CCa
{
protected:
	double g_Ca;
	double m_c;
public:
	__device__ CCa(double gca, double n);
	__device__ double Compute_n(double tempV, double TimeStep);
	__device__ double Compute_KI(double tempV, double TimeStep);
	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
	__device__ CCa();
	__device__ ~CCa();
};
__device__ CCa::CCa(double gca, double n)
{
	g_Ca = gca;
	m_c = n;
}
__device__ double CCa::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double CCa::Compute_n(double tempV, double TimeStep)
{
	double an = 0.3*(tempV + 13) / (1 - exp(-(tempV + 13)*1.0 / 10));
	double bn = 10.0 * exp(-(tempV + 38) *1.0 / 18);
	m_c = Iterative(m_c, an, bn, TimeStep);
	return m_c;
}
__device__ double CCa::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double Ca_I = g_Ca * pow(tmp_n, 3)*(tempV - 120);
	return Ca_I;
}
__device__ CCa::CCa()
{}
__device__ CCa::~CCa()
{
}

/////////////////////////////////////////////////////K
class CK
{
protected:
	double g_K;
	double m_Kn;
public:
	__device__ CK(double gk, double n);
	__device__ double Compute_n(double tempV, double TimeStep);
	__device__ double Compute_KI(double tempV, double TimeStep);
	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
	__device__ CK();
	__device__ ~CK();
};
__device__ CK::CK(double gk, double n)
{
	g_K = gk;
	m_Kn = n;
}
__device__ double CK::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double CK::Compute_n(double tempV, double TimeStep)
{
	double an = 0.01*(tempV + 55) / (1 - exp(-0.1*(tempV + 55)));
	double bn = 0.125*exp(-(tempV + 65) *1.0 / 80);
	m_Kn = Iterative(m_Kn, an, bn, TimeStep);
	return m_Kn;
}
__device__ double CK::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double K_I = g_K * pow(tmp_n, 4)*(tempV + 77);
	return K_I;
}
__device__ CK::CK()
{}
__device__ CK::~CK()
{
}
////////////////////////////////////////////////////KM
class CKM
{
protected:
	double g_KM;
	double m_z;
public:
	__device__ CKM(double gk, double n);
	__device__ double Compute_n(double tempV, double TimeStep);
	__device__ double Compute_KI(double tempV, double TimeStep);
	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
	__device__ CKM();
	__device__ ~CKM();
};

__device__ CKM::CKM(double gk, double n)
{
	g_KM = gk;
	m_z = n;
}
__device__ double CKM::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double CKM::Compute_n(double tempV, double TimeStep)
{
	double an = 1.0 / (1 + exp(-0.2*(tempV + 39))) / 75;
	double bn = 1.0 / 75 - an;
	m_z = Iterative(m_z, an, bn, TimeStep);
	return m_z;
}
__device__ double CKM::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double K_I = g_KM * pow(tmp_n, 1)*(tempV + 77);
	return K_I;
}
__device__ CKM::CKM()
{}
__device__ CKM::~CKM()
{
}
////////////////////////////////////////////////////Kv
class CKv
{
protected:
	double g_Kv;
	double m_p;
public:
	__device__ CKv(double gk, double n);
	__device__ double Compute_n(double tempV, double TimeStep);
	__device__ double Compute_KI(double tempV, double TimeStep);
	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
	__device__ CKv();
	__device__ ~CKv();
};

__device__ CKv::CKv(double gk, double n)
{
	g_Kv = gk;
	m_p = n;
}
__device__ double CKv::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double CKv::Compute_n(double tempV, double TimeStep)
{
	double an = (tempV - 95)*1.0 / (1 - exp(-(tempV - 95) / 11.8));
	double bn = 0.025*exp(-tempV / 22.222);
	m_p = Iterative(m_p, an, bn, TimeStep);
	return m_p;
}
__device__ double CKv::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double K_I = g_Kv * pow(tmp_n, 2)*(tempV + 77);
	return K_I;
}
__device__ CKv::CKv()
{}
__device__ CKv::~CKv()
{
}

//////////////////////////////////////////////////////Na

class CNa
{
protected:
	double g_Na;
	double m_Nam;
	double m_Nah;
public:
	__device__ CNa(double gna, double m, double h);
	__device__ double Compute_m(double tempV, double TimeStep);
	__device__ double Compute_h(double tempV, double TimeStep);
	__device__ double Compute_NaI(double tempV, double TimeStep);
	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
	__device__ CNa();
	__device__ ~CNa();
};

__device__ CNa::CNa(double gna, double m, double h)
{
	g_Na = gna;
	m_Nam = m;
	m_Nah = h;
}
__device__ double CNa::Compute_m(double tempV, double TimeStep)
{
	double am = 0.1*(tempV + 40) / (1 - exp(-0.1*(tempV + 40)));
	double bm = 4 * exp(-(tempV + 65) *1.0 / 18);//////////18
	m_Nam = Iterative(m_Nam, am, bm, TimeStep);
	return m_Nam;
}
__device__ double CNa::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double CNa::Compute_h(double tempV, double TimeStep)
{
	double ah = 0.07*exp(-(tempV + 65)*1.0 / 20);//////////20
	double bh = 1.0 / (exp(-0.1*(tempV + 35)) + 1);
	m_Nah = Iterative(m_Nah, ah, bh, TimeStep);
	return m_Nah;
}
__device__ double CNa::Compute_NaI(double tempV, double TimeStep)
{
	double tmp_m = Compute_m(tempV, TimeStep);
	double tmp_h = Compute_h(tempV, TimeStep);
	double Na_I = g_Na * pow(tmp_m, 3)*tmp_h*(tempV - 55);
	return Na_I;
}
__device__ CNa::CNa()
{}
__device__ CNa::~CNa()
{}

////////////////////////////////////////////////////Cell
class CCell
{
protected:
	double m_start;
	double m_Current;
	double m_TimeStep;
public:
	CNa CNa_a;
	CK CK_a;
	CKM CKM_a;
	CKv CKv_a;
	CCa CCa_a;
	float *m_pArray_V;

public:
	__device__ CCell(double start, double Current, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a);
	__device__ double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep);

	__device__ double  CCell::GetStart();
	__device__ double  CCell::GetTimeStep();
	__device__ double  CCell::GetCurrent();
	//CStringArray CCell::GetPArray_V();

	__device__ ~CCell();
};

__device__ CCell::CCell(double start, double Current, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a)
{
	m_start = start;
	CNa_a = Na_a;
	CK_a = K_a;
	CKM_a = KM_a;
	CKv_a = Kv_a;
	CCa_a = Ca_a;
	//m_pArray_V = pArray_V;
	m_Current = Current;
	m_TimeStep = TimeStep;
}
__device__ double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
__device__ double  CCell::GetStart()
{
	double start = m_start;
	return start;
}

__device__ double  CCell::GetTimeStep()
{
	double 	 TimeStep = m_TimeStep;
	return TimeStep;
}
__device__ double  CCell::GetCurrent()
{
	double 	 Current = m_Current;
	return Current;
}


__device__ CCell::~CCell()
{
}

///////////////////////////////////////////////////
void cudasafe(cudaError_t error, char* message = "Error occured") {
	if (error != cudaSuccess) {
		fprintf(stderr, "ERROR: %s : %i\n", message, error);
		exit(-1);
	}
}
__device__ void Calculation_I(CCell Cell_B, double Mtime, int FlagParameter[], float gL, float C, float *pArrayA)
{

	double temp, y_preB, y_nextB;
	double cell_B_K_I;
	double cell_B_Na_I;
	double cell_B_Ca_I;
	double cell_B_KM_I;
	double cell_B_Kv_I;
	//double cell_pre_Nmda_a = Cell_pre.CNmda_a.CNMDA_I(tempV_a, tempVB, Time, Ta, Tb, Times_b, Cell_B.GetTimeStep());
	///////
	//////////////////////////////////
	//Cell_B.m_pArray_V = new float[int(Mtime) * 100];
	double tempVB = Cell_B.GetStart();
	for (int k = 0; k < Mtime / Cell_B.GetTimeStep(); k++)
	{
		if (FlagParameter[1] == 0)
			cell_B_K_I = 0;
		else
			cell_B_K_I = Cell_B.CK_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
		if (FlagParameter[0] == 0)
			cell_B_Na_I = 0;
		else
			cell_B_Na_I = Cell_B.CNa_a.Compute_NaI(tempVB, Cell_B.GetTimeStep());
		if (FlagParameter[4] == 0)
			cell_B_Ca_I = 0;
		else
			cell_B_Ca_I = Cell_B.CCa_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
		if (FlagParameter[2] == 0)
			cell_B_KM_I = 0;
		else
			cell_B_KM_I = Cell_B.CKM_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
		if (FlagParameter[3] == 0)
			cell_B_Kv_I = 0;
		else
			cell_B_Kv_I = Cell_B.CKv_a.Compute_KI(tempVB, Cell_B.GetTimeStep());

		y_preB = tempVB + Cell_B.GetTimeStep()*(gL*(Cell_B.GetStart() - tempVB) + Cell_B.GetCurrent()
			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I) / C;

		cell_B_K_I = Cell_B.CK_a.Compute_KI(y_preB, Cell_B.GetTimeStep());

		cell_B_Na_I = Cell_B.CNa_a.Compute_NaI(y_preB, Cell_B.GetTimeStep());
		cell_B_Ca_I = Cell_B.CCa_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
		cell_B_KM_I = Cell_B.CKM_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
		cell_B_Kv_I = Cell_B.CKv_a.Compute_KI(y_preB, Cell_B.GetTimeStep());

		y_nextB = tempVB + Cell_B.GetTimeStep()*(gL*(Cell_B.GetStart() - y_preB) + Cell_B.GetCurrent()
			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I) / C;

		temp = (y_preB + y_nextB)*1.0 / 2;
		pArrayA[k] = temp;
		tempVB = temp;
	}
}


__global__ void HH_EntranceA(float *population, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA, int VAR_NUMBER)
{
	//freopen("V_output.txt", "w", stdout);
	//double Mtime = 100;
	//double tempVB = -69;
	//double TimeStep = 0.01;
	//printf("v_%f ", tempVB);
	int spike_length = int(Mtime / TimeStep);
	const int gid = blockDim.x * blockIdx.x + threadIdx.x;
	//const int tid = threadIdx.x;
	printf("gid_%d ", gid);
	float gNa = population[gid*VAR_NUMBER + 0];
	float gK = population[gid*VAR_NUMBER + 1];
	float gKM = population[gid*VAR_NUMBER + 2];
	float gKv = population[gid*VAR_NUMBER + 3];
	float gCa = population[gid*VAR_NUMBER + 4];
	//printf("na_%f ", gNa);
	CK K_a(gK, 0.318);//36
	CNa Na_a(gNa, 0.053, 0.596);///(120, 0.053, 0.596)
	CKM KM_a(gKM, 0);///(120, 0.053, 0.596)
	CKv Kv_a(gKv, 0);///(120, 0.053, 0.596)
	CCa Ca_a(gCa, 0);///(120, 0.053, 0.596)
	//double m_I = 10;
	//printf("s_%d ", 2);
	CCell Cell_D(tempVB, m_I, TimeStep, Na_a, K_a, KM_a, Kv_a, Ca_a);
	Calculation_I(Cell_D, Mtime, FlagParameter, gL, C, &pArrayA[gid*spike_length]);
	__syncthreads();
	//for (int i = 0; i < int(Mtime) * 100; i++)
	//{
	//	printf("%f\n", pArrayA[i]);
	//}
	////////////////
	//delete pArrayA;

}
extern "C" void HH_Entrance(float *population, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA, int POPULATION_SIZE)
{
	int VAR_NUMBER = 5;
	float *devicePopulation = 0;
	float *devicepArrayA = 0;
	int *deviceFlagParameter = 0;
	const int THREADS_PER_BLOCK = 256;
	const int BLOCKS_NUMBER = (POPULATION_SIZE + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
	int spike_length = int(Mtime / TimeStep);
	cudasafe(cudaMalloc((void **)&devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float)), "Could not allocate memory for devicePopulation");
	cudasafe(cudaMalloc((void **)&devicepArrayA, POPULATION_SIZE * spike_length * sizeof(float)), "Could not allocate memory for devicepArrayA");
	cudasafe(cudaMalloc((void **)&deviceFlagParameter, VAR_NUMBER* sizeof(int)), "Could not allocate memory for deviceFlagParameter");

	//cudasafe(cudaMemcpy(devicepArrayA, pArrayA, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy pArrayA to device");
	cudasafe(cudaMemcpy(devicePopulation, population, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy population to device");
	cudasafe(cudaMemcpy(deviceFlagParameter, FlagParameter, VAR_NUMBER * sizeof(int), cudaMemcpyHostToDevice), "Could not copy FlagParameter to device");

	HH_EntranceA << <BLOCKS_NUMBER, THREADS_PER_BLOCK >> >(devicePopulation, Mtime, tempVB, TimeStep, m_I, deviceFlagParameter, gL, C, devicepArrayA, VAR_NUMBER);
	cudasafe(cudaMemcpy(pArrayA, devicepArrayA, POPULATION_SIZE * spike_length * sizeof(float), cudaMemcpyDeviceToHost), "Could not copy pArrayA from device");

	cudasafe(cudaFree(devicepArrayA), "Could not free devicepArrayA");
	cudasafe(cudaFree(devicePopulation), "Could not free devicePopulation");
}


