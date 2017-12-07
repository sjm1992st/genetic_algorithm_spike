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

	const int gid = blockDim.x * blockIdx.x + threadIdx.x + THREADS_PER_BLOCK;
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
				////const float sign = SIGN[order < 0.0f];
				////sharedPopulation[tid + THREADS_PER_BLOCK][i] += powf(10.0f, order + order_deviation) * sign * mult;
				///sharedPopulation[tid + THREADS_PER_BLOCK][i] += sharedPopulation[tid + THREADS_PER_BLOCK][i] * order_deviation;
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
extern "C" float solveGPU(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult) {
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
		startMy = clock();
		start = clock();
		hDll = LoadLibrary("CHdll.dll");
		int spike_length = int(Mtime / TimeStep);
		float *temp_data = new float[spike_length*POPULATION_SIZE];
		HH_return(population, VAR_NUMBER, Mtime, tempVB, TimeStep, m_I, FlagParameter, gL, C, temp_data,POPULATION_SIZE);
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
		GAKernel_GenEach << <BLOCKS_NUMBER, THREADS_PER_BLOCK >> >(devicePopulation, deviceScore, randomStates, deviceParameter_, deviceParameter_Tset, tau, k, MaxGeneration, deviceParameter_Bound, POPULATION_SIZE, crossver, mutations);
		//delete[]population;
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time2 " << duration << endl;
		//float *population = new float[POPULATION_SIZE * VAR_NUMBER];
		cudasafe(cudaMemcpy(population, devicePopulation, POPULATION_SIZE * VAR_NUMBER * sizeof(float), cudaMemcpyDeviceToHost), "Could not copy population from device");

		//printf("%d_1111\n", k);
		//printFinalPopulation(population, deviceScore, POPULATION_SIZE);
		//printf("%d_2222\n", k);
		//while (MFlage)
		/////////////////////exit-code
		if (myexit)
		{
			cudasafe(cudaFree(deviceScore), "Failed to free deviceScore");
			cudasafe(cudaFree(devicePopulation), "Failed to free devicePopulation");

			cudasafe(cudaFree(randomStates), "Could not free randomStates");
			cudasafe(cudaFree(nextGeneration), "Could not free nextGeneration");
			cudasafe(cudaFree(deviceParameter_Tset), "Could not free deviceParameter_Tset");
			delete Parameter_Tset;
			delete population;
			return 0;
		}
		//////////////////////
		while (MFlage == 2){}
		MFlage = 1;
		start = clock();
		ans = FinallResult(k, tau, population, Mtime, tempVB, TimeStep, m_I, Parameter_, FlagParameter, gL, C, Parameter_Bound, POPULATION_SIZE, strResult);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time3 " << duration << endl;
		MFlage = 3;
		cudasafe(cudaFree(deviceScore), "Failed to free deviceScore");
		FreeLibrary(hDll);
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
    cudasafe(cudaGetLastError(), "Could not invoke GAKernel");
    cudasafe(cudaDeviceSynchronize(), "Failed to syncrhonize device after calling GAKernel");

    //printPopulation(devicePopulation, deviceScore);
    //}

	// freeing memory
	cudasafe(cudaFree(devicePopulation), "Failed to free devicePopulation");
	
	cudasafe(cudaFree(randomStates), "Could not free randomStates");
	cudasafe(cudaFree(nextGeneration), "Could not free nextGeneration");
	cudasafe(cudaFree(deviceParameter_Tset), "Could not free deviceParameter_Tset");
	delete Parameter_Tset;
	delete population;
	myexit = true;
	return ans;
}


