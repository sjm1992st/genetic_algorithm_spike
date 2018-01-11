#include <iostream>
#include <ctime>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <fstream>  
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math_constants.h>
#include <sstream>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <time.h> 
#include <string>
#include "struct_cu.h"
#include "constants.h"
//#include "buffer.h"
using namespace std;
// assume block size equal population size

void cudasafe(cudaError_t error, char* message = "Error occured") {
	if(error != cudaSuccess) {
		fprintf(stderr,"ERROR: %s : %i\n", message, error);
		exit(-1);
	}
}

__global__ void randomInit(curandState* state, unsigned long seed) {
	int tid = threadIdx.x;
	curand_init(seed, tid, 0, state + tid);
}

__device__ float fitness(M_args deviceParameter, M_args_Tset deviceParameter_Tset,float tau) {
    float result = 0;
	//printf("%d_a",deviceParameter.spike_data_num);
	//printf("%d_b", deviceParameter_Tset.length);
	for (size_t i = 0; i<deviceParameter_Tset.length; ++i)
		for (size_t j = 0; j<deviceParameter_Tset.length; ++j)
	{
		result += expf(-fabsf(deviceParameter_Tset.spike_TestData[i] - deviceParameter_Tset.spike_TestData[j])*1.0 / tau);
		//printf("%f_3 ", result);
		// ++curPos;
	}

	for (size_t i = 0; i<deviceParameter.spike_data_num; ++i)
		for (size_t j = 0; j<deviceParameter.spike_data_num; ++j)
		{
			result += expf(-fabsf(deviceParameter.spike_data[i] - deviceParameter.spike_data[j])*1.0 / tau);
			//printf("%f_2 ", result);
			// ++curPos;
		}
	for (size_t i = 0; i<deviceParameter.spike_data_num; ++i)
		for (size_t j = 0; j<deviceParameter_Tset.length; ++j)
		{
			//printf("%f_c ", deviceParameter.spike_data[i]);
			//printf("%f_d ", deviceParameter_Tset.spike_TestData[j]);
			result -= 2*expf(-fabsf(deviceParameter.spike_data[i] - deviceParameter_Tset.spike_TestData[j])*1.0 / tau);
			//printf("%f_1 ", result);
			// ++curPos;
		}
    //printf("%f_4 ", result);
    return result;
}


__global__ void GAKernel_GenEach(float* population, ScoreWithId* score, curandState* randomStates, M_args deviceParameter, M_args_Tset *deviceParameter_Tset, float tau, int genindex, int MaxGeneration, M_args_Bound *deviceParameter_Bound, const int POPULATION_SIZE, float crossver,float mutations) {
	// we first have to calculate the score for the first half of threads
	//const float *curPos = sharedPopulation[tid];
	__shared__ float sharedPopulation[THREADS_PER_BLOCK*2][VAR_NUMBER];
	//__shared__ float sharedScore[THREADS_PER_BLOCK * 2];
	//__shared__ float **sharedPopulation = new float*[THREADS_PER_BLOCK * 2];
	__shared__ float sharedScore[THREADS_PER_BLOCK*2];
	//const float SIGN[2] = { -1.0f, 1.0f };
	//const float MULT[2] = { 1.0f, 0.0f };

	const int gid = blockDim.x * blockIdx.x + threadIdx.x + THREADS_PER_BLOCK*blockIdx.x;
	const int tid = threadIdx.x;
	printf("i_%f ", crossver);
	// loading initial random population into shared memory
	if ((gid + THREADS_PER_BLOCK)< POPULATION_SIZE) {
		for (int i = 0; i < VAR_NUMBER; ++i){
			__syncthreads();

			sharedPopulation[tid][i] = population[gid * VAR_NUMBER + i];
			sharedPopulation[tid + THREADS_PER_BLOCK][i] = population[(gid + THREADS_PER_BLOCK) * VAR_NUMBER + i];
		}
		M_args_Tset curPos = deviceParameter_Tset[gid];
		sharedScore[tid] = fitness(deviceParameter, curPos, tau);
		M_args_Tset curPos_b = deviceParameter_Tset[(gid + THREADS_PER_BLOCK)];
		sharedScore[tid + THREADS_PER_BLOCK] = fitness(deviceParameter, curPos_b, tau);
		__syncthreads();
		}

	//sharedScore[tid + THREADS_PER_BLOCK] = 123123.0;
	curandState &localState = randomStates[(tid*genindex)%THREADS_PER_BLOCK];
	__syncthreads();


	//__syncthreads();
	if (genindex <= MaxGeneration && (gid) < POPULATION_SIZE)
	{// selection
		// first half of threads writes best individual into its position
			if (sharedScore[tid] > sharedScore[tid + THREADS_PER_BLOCK]) {
				for (int i = 0; i < VAR_NUMBER; ++i)
					sharedPopulation[tid][i] = sharedPopulation[tid + THREADS_PER_BLOCK][i];
				sharedScore[tid] = sharedScore[tid + THREADS_PER_BLOCK];
			}
		
		__syncthreads();
		////int temp_size = THREADS_PER_BLOCK*(1 - crossver);
		////if (tid<temp_size && sharedScore[tid] > sharedScore[tid + temp_size]) {
		////	for (int i = 0; i < VAR_NUMBER; ++i)
		////		sharedPopulation[tid][i] = sharedPopulation[tid + temp_size][i];
		////	sharedScore[tid] = sharedScore[tid + temp_size];
		////}
		__syncthreads();
		// crossover
		////printf("i_%f ", crossver);
		////printf("gi_%f ", curand_uniform(&localState));
		if (curand_uniform(&localState) < crossver) {
			curandState &localState1 = randomStates[tid * 3];
			curandState &localState2 = randomStates[tid * 4];
			const int first = curand_uniform(&localState1) * THREADS_PER_BLOCK;
			const int second = curand_uniform(&localState2) * THREADS_PER_BLOCK;
			const float weight = curand_uniform(&localState1);
			for (int i = 0; i < VAR_NUMBER; ++i) {
				////const float Temp_weight = sharedPopulation[first][i] * weight + sharedPopulation[second][i] * (1.0f - weight);
				//////printf("i_%f ", Temp_weight);
				////if (Temp_weight>deviceParameter_Bound[i].g[1] || Temp_weight < deviceParameter_Bound[i].g[0])
				////	sharedPopulation[tid + THREADS_PER_BLOCK][i] = deviceParameter_Bound[i].g[0] + weight*(deviceParameter_Bound[i].g[1] - deviceParameter_Bound[i].g[0]);
				////else
				////	sharedPopulation[tid + THREADS_PER_BLOCK][i] = Temp_weight;
				//printf("i_%f %f", deviceParameter_Bound[i].g[0], deviceParameter_Bound[i].g[1]);
				//printf("i_%f ", Temp_weight);
				sharedPopulation[tid + THREADS_PER_BLOCK][i] = sharedPopulation[first][i] * weight + sharedPopulation[second][i] * (1.0f - weight);
				
			}
		}
		__syncthreads();

		// mutations on second half of population
		curandState &localState3= randomStates[tid*5];
		if (curand_uniform(&localState3) < mutations) {
			//const float order = (curand_uniform(&localState) - mutations*1.0/2);
			float guass = curand_normal(&localState3);
			for (int i = 0; i < VAR_NUMBER; ++i) {
				////const float mult = MULT[order < 0.0f];
				////if (sharedPopulation[tid + THREADS_PER_BLOCK][i]>deviceParameter_Bound[i].g[1] || sharedPopulation[tid + THREADS_PER_BLOCK][i] < deviceParameter_Bound[i].g[0])
				sharedPopulation[tid + THREADS_PER_BLOCK][i] = (deviceParameter_Bound[i].g[0] + deviceParameter_Bound[i].g[1]) / 2 + guass;
			}
		}
		__syncthreads();
		const int third = curand_uniform(&localState)*THREADS_PER_BLOCK;
		//////sharing a part of population with others
		////if (curand_uniform(&localState)<0.1) {
		////	for (int i = 0; i < VAR_NUMBER; ++i)
		////		population[gid * VAR_NUMBER + i] = sharedPopulation[tid + THREADS_PER_BLOCK][i];
		////}

		// take some best individuals from neighbour
		if ((blockIdx.x + third) % 3 == 0) {
			if (curand_uniform(&localState) < 0.11) {
				const int anotherBlock = curand_uniform(&localState) * (POPULATION_SIZE / THREADS_PER_BLOCK);
				const int ngid = blockDim.x * anotherBlock + threadIdx.x;
				for (int i = 0; i < VAR_NUMBER; ++i)
					sharedPopulation[tid][i] = population[ngid * VAR_NUMBER + i];
				//sharedScore[tid] = fitness(sharedPopulation[tid], deviceParameter);
				//sharedScore[tid]=fitness(deviceParameter, curPos_b, tau);
			}
		}
	}
	__syncthreads();
	////// output current population back
	if ((gid ) < POPULATION_SIZE) {
		for (int i = 0; i < VAR_NUMBER; ++i)
			{
				population[gid * VAR_NUMBER + i] = sharedPopulation[tid][i];
				population[(gid + THREADS_PER_BLOCK) * VAR_NUMBER + i] = sharedPopulation[tid + THREADS_PER_BLOCK][i];
				printf("%f ", population[gid * VAR_NUMBER + i]);
			}
	
	if (genindex <= MaxGeneration)
		{
			score[gid].score = sharedScore[tid];
		}
	}
}


void printFinalPopulation(float* population, const ScoreWithId* deviceScore, const int POPULATION_SIZE) {
	////float **population = new float *[POPULATION_SIZE];
	////for (int i = 0; i < POPULATION_SIZE; i++)
	////	population[i] = new float[VAR_NUMBER];
	////cudasafe(cudaMemcpy(population, devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyDeviceToHost), "Could not copy population from device");

	ScoreWithId *score = new ScoreWithId[POPULATION_SIZE];
	cudasafe(cudaMemcpy(score, deviceScore, POPULATION_SIZE * sizeof (ScoreWithId), cudaMemcpyDeviceToHost), "Could not copy score to host");

	///std::cout.cetf(std::ios::fixed);
	std::cout.precision(12);

	for (int i = 0; i<1; ++i) {
		std::cout <<"gen_index"<< std::setw(15) << i << ' ';
	}
	std::cout << std::endl;

	for (int i=0; i<VAR_NUMBER; i++) {
			for (int u=0; u<POPULATION_SIZE; ++u) {
				std::cout << std::setw(15) << population[u*VAR_NUMBER + i] << ' ';
			}
			std::cout << std::endl;
		}
	std::cout << "Score: " << std::endl;
	for (int i = 0; i<POPULATION_SIZE; ++i) {
		std::cout << std::setw(15) << score[i].score << ' ';
	}
	std::cout << std::endl;
	//delete population;
}
extern "C" float solveGPU(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, vector<float> m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult) {
    cudasafe(cudaSetDevice(0), "Could not set device 0");
	float ans = 0;
	clock_t start, finish, startMy, finishMy;
	float tau = 12;
	//M_args *IndexParameter_ = new M_args[MaxGeneration];
	//IndexParameter_ = 0;
	//IndexParameter_[0] = Parameter_;
	/////////////////////////////////////
	////cout << Parameter_.spike_data_num << endl;
	////for (int i = 0; i < Parameter_.spike_data_num; i++)
	////	cout << Parameter_.spike_data[i] << " ";
	////cout << endl;
	////cout << mutations << " " << crossver<<endl;
	////std::cout << std::setw(15) << Parameter_Bound[0].g[0] << " " << Parameter_Bound[0].g[1] << " " << tempVB << endl;
	std::cout << FlagParameter[0] << ' ' << FlagParameter[1] << ' ' << FlagParameter[2] << ' ' << FlagParameter[3] << ' ' << FlagParameter[4] << endl;
	std::cout << Parameter_Bound[0].g[0] << ' ' << Parameter_Bound[1].g[0] << ' ' << Parameter_Bound[2].g[0] <<' ' << Parameter_Bound[3].g[0]<<endl;
	//////////////////////////////////////
	float *population = new float[POPULATION_SIZE*VAR_NUMBER];
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<POPULATION_SIZE; i++) {
		for (int j=0; j<VAR_NUMBER; j++) {
			if (FlagParameter[j] == 0)
			{
				population[i*VAR_NUMBER + j] = 0;
			}
			else
			{
				population[i*VAR_NUMBER + j] = (float_random(Parameter_Bound[j]));
				//std::cout << std::setw(15) << population[i*VAR_NUMBER + j] << ' ';
			}
		}
		//std::cout << endl;
	}

	M_args_Tset *Parameter_Tset=new M_args_Tset[POPULATION_SIZE];

	// copying population to device
	float *devicePopulation = 0;
	float *nextGeneration = 0;
	M_args_Tset *deviceParameter_Tset = 0;
	
	curandState* randomStates;
	M_args deviceParameter_;
	//deviceParameter_.current_data_num = Parameter_.current_data_num;
	deviceParameter_.spike_data_num = Parameter_.spike_data_num;
	M_args_Bound *deviceParameter_Bound = 0;
	//int DataLength = getArrayLen(Parameter_.spike_data);


	//int DataLengthC = getArrayLen(Parameter_.current_data);
	//Parameter_.length = DataLength;
	cudasafe(cudaMalloc((void **)&deviceParameter_Bound, VAR_NUMBER * sizeof(M_args_Bound)), "Could not allocate memory for deviceParameter_Bound");

	cudasafe(cudaMalloc(&randomStates, THREADS_PER_BLOCK * sizeof(curandState)), "Could not allocate memory for randomStates");
	cudasafe(cudaMalloc((void **)&devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float)), "Could not allocate memory for devicePopulation");
	cudasafe(cudaMalloc((void **)&nextGeneration, POPULATION_SIZE * VAR_NUMBER * sizeof(float)), "Could not allocate memory for nextGeneration");
	cudasafe(cudaMalloc((void **)&deviceParameter_Tset, 2*POPULATION_SIZE * sizeof (M_args_Tset)), "Could not allocate memory for deviceParameter_Tset");
	//cudasafe(cudaMalloc((void **)&deviceParameter_.current_data, Parameter_.current_data_num*sizeof(float)), "Could not allocate memory for deviceParameter_");
	cudasafe(cudaMalloc((void **)&deviceParameter_.spike_data, Parameter_.spike_data_num*sizeof(float)), "Could not allocate memory for deviceParameter_");
	//cudasafe(cudaMalloc((void **)&deviceParameter_, sizeof(M_args)), "Could not allocate memory for deviceParameter_");
	//cudasafe(cudaMalloc((void **)&deviceParameter_.spike_TestData, DataLength*sizeof(float)), "Could not allocate memory for deviceParameter_");

	///cudasafe(cudaMemcpy(devicePopulation, population, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy population to device");
	//cudasafe(cudaMemcpy(deviceParameter_.current_data, Parameter_.current_data, Parameter_.current_data_num*sizeof(float), cudaMemcpyHostToDevice), "Could not copy Parameter_current_data to device");
	cudasafe(cudaMemcpy(deviceParameter_.spike_data, Parameter_.spike_data, Parameter_.spike_data_num*sizeof(float), cudaMemcpyHostToDevice), "Could not copy Parameter_spike_data to device");

	//cudasafe(cudaMemcpy(deviceParameter_.spike_TestData, Parameter_.spike_TestData, DataLength*sizeof(float), cudaMemcpyHostToDevice), "Could not copy Parameter_ to device");
	for (int kb = 0; kb<VAR_NUMBER;kb++)
		cudasafe(cudaMemcpy(&deviceParameter_Bound[kb], &Parameter_Bound[kb], 2 * sizeof(float),cudaMemcpyHostToDevice), "Could not allocate memory for deviceParameter_Bound");


	// invoking random init
	randomInit<<<1, THREADS_PER_BLOCK>>>(randomStates, 900);
	cudasafe(cudaGetLastError(), "Could not invoke kernel randomInit");
	cudasafe(cudaDeviceSynchronize(), "Failed to syncrhonize device after calling randomInit");

	const int BLOCKS_NUMBER = (POPULATION_SIZE + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

	
	//__shared__ float sharedPopulation[THREADS_PER_BLOCK * 2][VAR_NUMBER];
	//__shared__ float sharedScore[THREADS_PER_BLOCK * 2];
    //for (int i=0; i<1115; i++) {
	//void GAKernel_GenEach(float* population, ScoreWithId* score, curandState* randomStates, M_args deviceParameter, M_args_Tset *deviceParameter_Tset, float tau)
	for (int k = 1; k <=MaxGeneration; k++) 
	{
		m_gen = k;
		hDll = LoadLibrary("CHdll.dll");
		startMy = clock();
		start = clock();
		int spike_length = int(Mtime / TimeStep);
		float *temp_data = new float[spike_length*POPULATION_SIZE];
		HH_return(population, VAR_NUMBER, Mtime, tempVB, TimeStep, m_I, FlagParameter, gL, C, temp_data,POPULATION_SIZE);
		FreeLibrary(hDll);
		//Parameter_Tset[j].spike_TestData = HH_SpikeTime(temp_data, Mtime, TimeStep, Parameter_Tset[j].length, Parameter_.spike_data_num);
		/////////////
		////for (int j = 0; j < POPULATION_SIZE; j++)
		////	for (int i = 0; i < spike_length; i++)
		////		cout << temp_data[j*spike_length+i] << endl;
		//////////////////
		//#pragma omp parallel for
		for (int j = 0; j < POPULATION_SIZE; j++)
		{
			float *temp_spike_TestData;
			vector<float>mdata_tmp = HH_SpikeTime(&temp_data[j*spike_length], Mtime, TimeStep, Parameter_Tset[j].length, Parameter_.spike_data_num);
			Parameter_Tset[j].spike_TestData = new float[mdata_tmp.size()];
			convert_data_two(mdata_tmp, Parameter_Tset[j].spike_TestData);
			cudaMalloc(&temp_spike_TestData, Parameter_Tset[j].length*sizeof(float));
			//std::cout << Parameter_Tset[j].length << std::endl;
			cudasafe(cudaMemcpy(&deviceParameter_Tset[j], &Parameter_Tset[j], sizeof (M_args_Tset), cudaMemcpyHostToDevice), "Could not copy deviceParameter_Tset1 to device");
			cudasafe(cudaMemcpy(temp_spike_TestData, Parameter_Tset[j].spike_TestData, (Parameter_Tset[j].length*sizeof(float)), cudaMemcpyHostToDevice), "Could not copy deviceParameter_Tset_spike_TestData2 to device");
			cudasafe(cudaMemcpy(&deviceParameter_Tset[j].spike_TestData, &temp_spike_TestData, sizeof(float*), cudaMemcpyHostToDevice), "Could not copy deviceParameter_Tset_spike_TestData to device");
			cudasafe(cudaFree(temp_spike_TestData), "Could not free temp_spike_TestData");
			delete Parameter_Tset[j].spike_TestData;
			//cout << j << endl;
		}
		delete temp_data;
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time1 " << duration << endl;
		ScoreWithId *deviceScore = 0;
		cudasafe(cudaMalloc((void **)&deviceScore, POPULATION_SIZE * sizeof (ScoreWithId)), "Could not allocate memory for deviceScore");
		cudasafe(cudaMemcpy(devicePopulation, population, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy population to device");

		start = clock();
		int A = 0;

		GAKernel_GenEach << <BLOCKS_NUMBER, THREADS_PER_BLOCK >> >(devicePopulation, deviceScore, randomStates, deviceParameter_, deviceParameter_Tset, tau, k, MaxGeneration, deviceParameter_Bound, POPULATION_SIZE, crossver, mutations);
		//delete[]population;
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time2 " << duration << endl;
		//float *population = new float[POPULATION_SIZE * VAR_NUMBER];
		cudasafe(cudaMemcpy(population, devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyDeviceToHost), "Could not copy population from device");
		cudasafe(cudaFree(deviceScore), "Failed to free deviceScore");
		//printf("%d_1111\n", k);
		//printFinalPopulation(population, deviceScore, POPULATION_SIZE);
		//printf("%d_2222\n", k);
		//while (MFlage)
		/////////////////////exit-code
		if (myexit)
		{
			cudasafe(cudaFree(devicePopulation), "Failed to free devicePopulation");

			cudasafe(cudaFree(randomStates), "Could not free randomStates");
			cudasafe(cudaFree(nextGeneration), "Could not free nextGeneration");
			cudasafe(cudaFree(deviceParameter_Tset), "Could not free deviceParameter_Tset");
			delete Parameter_Tset;
			delete population;
			return 0;
		}
		///////////////////////
		cudasafe(cudaGetLastError(), "Could not invoke GAKernel");
		cudasafe(cudaDeviceSynchronize(), "Failed to syncrhonize device after calling GAKernel");

		//printPopulation(devicePopulation, deviceScore);
		//}


		mtx.lock();
		start = clock();
		ans = FinallResult(k, tau, population, Mtime, tempVB, TimeStep, m_I, Parameter_, FlagParameter, gL, C, Parameter_Bound, POPULATION_SIZE, strResult);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time3 " << duration << endl;
		//MFlage = 3;
		mtx.unlock();
		finishMy = clock();
		
		durationMy = (double)(finishMy - startMy) / CLOCKS_PER_SEC;
		cout << "timeAll " << durationMy << endl;
		/////////////////////exit-code
		if (myexit)
		{
			cudasafe(cudaFree(devicePopulation), "Failed to free devicePopulation");

			cudasafe(cudaFree(randomStates), "Could not free randomStates");
			cudasafe(cudaFree(nextGeneration), "Could not free nextGeneration");
			cudasafe(cudaFree(deviceParameter_Tset), "Could not free deviceParameter_Tset");
			delete Parameter_Tset;
			delete population;
			return 0;
		}
		
		//////////////////////
		////GAKernel_gen << <BLOCKS_NUMBER, THREADS_PER_BLOCK >> >(devicePopulation, sharedPopulation, sharedScore, deviceScore, randomStates, deviceParameter_, deviceParameter_Tset, tau);
	}
	
	//ans = FinallResult(tau, population, Mtime, tempVB, TimeStep, m_I, Parameter_, FlagParameter, gL, C, Parameter_Bound, POPULATION_SIZE, strResult);
	// freeing memory
	cudasafe(cudaFree(devicePopulation), "Failed to free devicePopulation");

	cudasafe(cudaFree(randomStates), "Could not free randomStates");
	cudasafe(cudaFree(nextGeneration), "Could not free nextGeneration");
	cudasafe(cudaFree(deviceParameter_Tset), "Could not free deviceParameter_Tset");
	
	delete Parameter_Tset;
	//////////////////////
	delete population;
	myexit = true;
	return ans;
}











////////////////////////////////////////////////////////////\\gpu_HH
/////////////////////////////////////////////////////////////\\
///////////////////////////////////////////////////////////////
/////////////////////////
///////////////////////////////////////////////Ca
////class CCa
////{
////protected:
////	double g_Ca;
////	double m_c;
////public:
////	__device__ CCa(double gca, double n);
////	__device__ double Compute_n(double tempV, double TimeStep);
////	__device__ double Compute_KI(double tempV, double TimeStep);
////	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
////	__device__ CCa();
////	__device__ ~CCa();
////};
////__device__ CCa::CCa(double gca, double n)
////{
////	g_Ca = gca;
////	m_c = n;
////}
////__device__ double CCa::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double CCa::Compute_n(double tempV, double TimeStep)
////{
////	double an = 0.3*(tempV + 13) / (1 - exp(-(tempV + 13)*1.0 / 10));
////	double bn = 10.0 * exp(-(tempV + 38) *1.0 / 18);
////	m_c = Iterative(m_c, an, bn, TimeStep);
////	return m_c;
////}
////__device__ double CCa::Compute_KI(double tempV, double TimeStep)
////{
////	double tmp_n = Compute_n(tempV, TimeStep);
////	double Ca_I = g_Ca * pow(tmp_n, 3)*(tempV - 120);
////	return Ca_I;
////}
////__device__ CCa::CCa()
////{}
////__device__ CCa::~CCa()
////{
////}
////
/////////////////////////////////////////////////////////K
////class CK
////{
////protected:
////	double g_K;
////	double m_Kn;
////public:
////	__device__ CK(double gk, double n);
////	__device__ double Compute_n(double tempV, double TimeStep);
////	__device__ double Compute_KI(double tempV, double TimeStep);
////	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
////	__device__ CK();
////	__device__ ~CK();
////};
////__device__ CK::CK(double gk, double n)
////{
////	g_K = gk;
////	m_Kn = n;
////}
////__device__ double CK::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double CK::Compute_n(double tempV, double TimeStep)
////{
////	double an = 0.01*(tempV + 55) / (1 - exp(-0.1*(tempV + 55)));
////	double bn = 0.125*exp(-(tempV + 65) *1.0 / 80);
////	m_Kn = Iterative(m_Kn, an, bn, TimeStep);
////	return m_Kn;
////}
////__device__ double CK::Compute_KI(double tempV, double TimeStep)
////{
////	double tmp_n = Compute_n(tempV, TimeStep);
////	double K_I = g_K * pow(tmp_n, 4)*(tempV + 77);
////	return K_I;
////}
////__device__ CK::CK()
////{}
////__device__ CK::~CK()
////{
////}
////////////////////////////////////////////////////////KM
////class CKM
////{
////protected:
////	double g_KM;
////	double m_z;
////public:
////	__device__ CKM(double gk, double n);
////	__device__ double Compute_n(double tempV, double TimeStep);
////	__device__ double Compute_KI(double tempV, double TimeStep);
////	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
////	__device__ CKM();
////	__device__ ~CKM();
////};
////
////__device__ CKM::CKM(double gk, double n)
////{
////	g_KM = gk;
////	m_z = n;
////}
////__device__ double CKM::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double CKM::Compute_n(double tempV, double TimeStep)
////{
////	double an = 1.0 / (1 + exp(-0.2*(tempV + 39))) / 75;
////	double bn = 1.0 / 75 - an;
////	m_z = Iterative(m_z, an, bn, TimeStep);
////	return m_z;
////}
////__device__ double CKM::Compute_KI(double tempV, double TimeStep)
////{
////	double tmp_n = Compute_n(tempV, TimeStep);
////	double K_I = g_KM * pow(tmp_n, 1)*(tempV + 77);
////	return K_I;
////}
////__device__ CKM::CKM()
////{}
////__device__ CKM::~CKM()
////{
////}
////////////////////////////////////////////////////////Kv
////class CKv
////{
////protected:
////	double g_Kv;
////	double m_p;
////public:
////	__device__ CKv(double gk, double n);
////	__device__ double Compute_n(double tempV, double TimeStep);
////	__device__ double Compute_KI(double tempV, double TimeStep);
////	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
////	__device__ CKv();
////	__device__ ~CKv();
////};
////
////__device__ CKv::CKv(double gk, double n)
////{
////	g_Kv = gk;
////	m_p = n;
////}
////__device__ double CKv::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double CKv::Compute_n(double tempV, double TimeStep)
////{
////	double an = (tempV - 95)*1.0 / (1 - exp(-(tempV - 95) / 11.8));
////	double bn = 0.025*exp(-tempV / 22.222);
////	m_p = Iterative(m_p, an, bn, TimeStep);
////	return m_p;
////}
////__device__ double CKv::Compute_KI(double tempV, double TimeStep)
////{
////	double tmp_n = Compute_n(tempV, TimeStep);
////	double K_I = g_Kv * pow(tmp_n, 2)*(tempV + 77);
////	return K_I;
////}
////__device__ CKv::CKv()
////{}
////__device__ CKv::~CKv()
////{
////}
////
//////////////////////////////////////////////////////////Na
////
////class CNa
////{
////protected:
////	double g_Na;
////	double m_Nam;
////	double m_Nah;
////public:
////	__device__ CNa(double gna, double m, double h);
////	__device__ double Compute_m(double tempV, double TimeStep);
////	__device__ double Compute_h(double tempV, double TimeStep);
////	__device__ double Compute_NaI(double tempV, double TimeStep);
////	__device__ double Iterative(double tmp_n, double an, double bn, double TimeStep);
////	__device__ CNa();
////	__device__ ~CNa();
////};
////
////__device__ CNa::CNa(double gna, double m, double h)
////{
////	g_Na = gna;
////	m_Nam = m;
////	m_Nah = h;
////}
////__device__ double CNa::Compute_m(double tempV, double TimeStep)
////{
////	double am = 0.1*(tempV + 40) / (1 - exp(-0.1*(tempV + 40)));
////	double bm = 4 * exp(-(tempV + 65) *1.0 / 18);//////////18
////	m_Nam = Iterative(m_Nam, am, bm, TimeStep);
////	return m_Nam;
////}
////__device__ double CNa::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double CNa::Compute_h(double tempV, double TimeStep)
////{
////	double ah = 0.07*exp(-(tempV + 65)*1.0 / 20);//////////20
////	double bh = 1.0 / (exp(-0.1*(tempV + 35)) + 1);
////	m_Nah = Iterative(m_Nah, ah, bh, TimeStep);
////	return m_Nah;
////}
////__device__ double CNa::Compute_NaI(double tempV, double TimeStep)
////{
////	double tmp_m = Compute_m(tempV, TimeStep);
////	double tmp_h = Compute_h(tempV, TimeStep);
////	double Na_I = g_Na * pow(tmp_m, 3)*tmp_h*(tempV - 55);
////	return Na_I;
////}
////__device__ CNa::CNa()
////{}
////__device__ CNa::~CNa()
////{}
////
////////////////////////////////////////////////////////Cell
////class CCell
////{
////protected:
////	double m_start;
////	double m_TimeStep;
////public:
////	CNa CNa_a;
////	CK CK_a;
////	CKM CKM_a;
////	CKv CKv_a;
////	CCa CCa_a;
////	float *m_pArray_V;
////
////public:
////	__device__ CCell(double start, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a);
////	__device__ double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep);
////
////	__device__ double  CCell::GetStart();
////	__device__ double  CCell::GetTimeStep();
////	//CStringArray CCell::GetPArray_V();
////
////	__device__ ~CCell();
////};
////
////__device__ CCell::CCell(double start, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a)
////{
////	m_start = start;
////	CNa_a = Na_a;
////	CK_a = K_a;
////	CKM_a = KM_a;
////	CKv_a = Kv_a;
////	CCa_a = Ca_a;
////	//m_pArray_V = pArray_V;
////	m_TimeStep = TimeStep;
////}
////__device__ double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep)
////{
////	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
////	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
////	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
////	return tmp_avgn;
////}
////__device__ double  CCell::GetStart()
////{
////	double start = m_start;
////	return start;
////}
////
////__device__ double  CCell::GetTimeStep()
////{
////	double 	 TimeStep = m_TimeStep;
////	return TimeStep;
////}
////
////
////
////__device__ CCell::~CCell()
////{
////}
////
///////////////////////////////////////////////////////
////__device__ void Calculation_I(CCell Cell_B, float* pCurrent, double Mtime, int FlagParameter[], float gL, float C, float *pArrayA)
////{
////
////	double temp, y_preB, y_nextB;
////	double cell_B_K_I;
////	double cell_B_Na_I;
////	double cell_B_Ca_I;
////	double cell_B_KM_I;
////	double cell_B_Kv_I;
////	//double cell_pre_Nmda_a = Cell_pre.CNmda_a.CNMDA_I(tempV_a, tempVB, Time, Ta, Tb, Times_b, Cell_B.GetTimeStep());
////	///////
////	//////////////////////////////////
////	//Cell_B.m_pArray_V = new float[int(Mtime) * 100];
////	double tempVB = Cell_B.GetStart();
////	for (int k = 0; k < Mtime / Cell_B.GetTimeStep(); k++)
////	{
////		if (FlagParameter[1] == 0)
////			cell_B_K_I = 0;
////		else
////			cell_B_K_I = Cell_B.CK_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
////		if (FlagParameter[0] == 0)
////			cell_B_Na_I = 0;
////		else
////			cell_B_Na_I = Cell_B.CNa_a.Compute_NaI(tempVB, Cell_B.GetTimeStep());
////		if (FlagParameter[4] == 0)
////			cell_B_Ca_I = 0;
////		else
////			cell_B_Ca_I = Cell_B.CCa_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
////		if (FlagParameter[2] == 0)
////			cell_B_KM_I = 0;
////		else
////			cell_B_KM_I = Cell_B.CKM_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
////		if (FlagParameter[3] == 0)
////			cell_B_Kv_I = 0;
////		else
////			cell_B_Kv_I = Cell_B.CKv_a.Compute_KI(tempVB, Cell_B.GetTimeStep());
////
////		y_preB = tempVB + Cell_B.GetTimeStep()*(gL*(Cell_B.GetStart() - tempVB) + pCurrent[k]
////			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I) / C;
////
////		cell_B_K_I = Cell_B.CK_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
////
////		cell_B_Na_I = Cell_B.CNa_a.Compute_NaI(y_preB, Cell_B.GetTimeStep());
////		cell_B_Ca_I = Cell_B.CCa_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
////		cell_B_KM_I = Cell_B.CKM_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
////		cell_B_Kv_I = Cell_B.CKv_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
////
////		y_nextB = tempVB + Cell_B.GetTimeStep()*(gL*(Cell_B.GetStart() - y_preB) + pCurrent[k]
////			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I) / C;
////
////		temp = (y_preB + y_nextB)*1.0 / 2;
////		pArrayA[k] = temp;
////		tempVB = temp;
////	}
////}
////
////
////__global__ void HH_EntranceA(float *population, int POPULATION_SIZE, double Mtime, double tempVB, double TimeStep, float* pCurrent, int FlagParameter[], float gL, float C, float *pArrayA, int VAR_NUMBER)
////{
////	//freopen("V_output.txt", "w", stdout);
////	//double Mtime = 100;
////	//double tempVB = -69;
////	//double TimeStep = 0.01;
////	//printf("v_%f ", tempVB);
////	int spike_length = int(Mtime / TimeStep);
////	const int gid = blockDim.x * blockIdx.x + threadIdx.x;
////	//const int tid = threadIdx.x;
////	//if (gid < POPULATION_SIZE)
////	//{
////	//printf("gid_%d \n", gid);
////	float gNa = population[gid*VAR_NUMBER + 0];
////	float gK = population[gid*VAR_NUMBER + 1];
////	float gKM = population[gid*VAR_NUMBER + 2];
////	float gKv = population[gid*VAR_NUMBER + 3];
////	float gCa = population[gid*VAR_NUMBER + 4];
////	//printf("na_%f ", gNa);
////	CK K_a(gK, 0.318);//36
////	CNa Na_a(gNa, 0.053, 0.596);///(120, 0.053, 0.596)
////	CKM KM_a(gKM, 0);///(120, 0.053, 0.596)
////	CKv Kv_a(gKv, 0);///(120, 0.053, 0.596)
////	CCa Ca_a(gCa, 0);///(120, 0.053, 0.596)
////	//double m_I = 10;
////	//printf("s_%d ", 2);
////	CCell Cell_D(tempVB, TimeStep, Na_a, K_a, KM_a, Kv_a, Ca_a);
////	Calculation_I(Cell_D, pCurrent, Mtime, FlagParameter, gL, C, &pArrayA[gid*spike_length]);
////	//}
////	__syncthreads();
////	//for (int i = 0; i < int(Mtime) * 100; i++)
////	//{
////	//	printf("%f\n", pArrayA[i]);
////	//}
////	////////////////
////	//delete pArrayA;
////
////}
////void HH_Entrance(float *population, double Mtime, double tempVB, double TimeStep, float* m_I, int FlagParameter[], float gL, float C, float *pArrayA, int POPULATION_SIZE)
////{
////	int VAR_NUMBER = 5;
////	float *devicePopulation = 0;
////	float *deviceM_I = 0;
////	float *devicepArrayA = 0;
////	int *deviceFlagParameter = 0;
////	const int THREADS_PER_BLOCK = 1024;
////	const int BLOCKS_NUMBER = (POPULATION_SIZE + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
////	int spike_length = int(Mtime / TimeStep);
////	cudasafe(cudaMalloc((void **)&devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float)), "Could not allocate memory for devicePopulation");
////	cudasafe(cudaMalloc((void **)&deviceM_I, spike_length * sizeof(float)), "Could not allocate memory for deviceM_I");
////
////	cudasafe(cudaMalloc((void **)&devicepArrayA, POPULATION_SIZE * spike_length * sizeof(float)), "Could not allocate memory for devicepArrayA");
////	cudasafe(cudaMalloc((void **)&deviceFlagParameter, VAR_NUMBER* sizeof(int)), "Could not allocate memory for deviceFlagParameter");
////
////	//cudasafe(cudaMemcpy(devicepArrayA, pArrayA, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy pArrayA to device");
////	cudasafe(cudaMemcpy(devicePopulation, population, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyHostToDevice), "Could not copy population to device");
////	cudasafe(cudaMemcpy(deviceM_I, m_I, spike_length* sizeof(float), cudaMemcpyHostToDevice), "Could not copy m_I to device");
////
////	cudasafe(cudaMemcpy(deviceFlagParameter, FlagParameter, VAR_NUMBER * sizeof(int), cudaMemcpyHostToDevice), "Could not copy FlagParameter to device");
////
////	HH_EntranceA << <BLOCKS_NUMBER, THREADS_PER_BLOCK >> >(devicePopulation, POPULATION_SIZE, Mtime, tempVB, TimeStep, deviceM_I, deviceFlagParameter, gL, C, devicepArrayA, VAR_NUMBER);
////	cudasafe(cudaMemcpy(pArrayA, devicepArrayA, POPULATION_SIZE * spike_length * sizeof(float), cudaMemcpyDeviceToHost), "Could not copy pArrayA from device");
////
////	cudasafe(cudaFree(devicepArrayA), "Could not free devicepArrayA");
////	cudasafe(cudaFree(devicePopulation), "Could not free devicePopulation");
////	cudasafe(cudaFree(deviceM_I), "Could not free deviceM_I");
////	cudasafe(cudaFree(deviceFlagParameter), "Could not free deviceFlagParameter");
////}

////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////