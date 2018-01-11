#include <windows.h>
#include<vector>
#include <tchar.h>
#include <conio.h>  
#include <stdio.h>
#include <algorithm>
#include <map>
#include "omp.h"
#include <mutex> 

using namespace std;

///void HH_Entrance(float *population, double Mtime, double tempVB, double TimeStep, float* m_I, int FlagParameter[], float gL, float C, float *pArrayA, int POPULATION_SIZE);
typedef void(*AddFunc)(float *population, double Mtime, double tempVB, double TimeStep, float* m_I, int FlagParameter[], float gL, float C, float *pArrayA, int POPULATION_SIZE);
//HMODULE hDll = LoadLibrary("CHdll.dll");
HMODULE hDll;
///extern "C" _declspec(dllimport) void HH_Entrance(float gNa, float gK, float gKM, float gKv, float gCa, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *Array_Data);
#define CLOCKS_PER_SEC ((clock_t)1000) 
//#pragma comment(lib,"CHH_dll.lib")
const double KNOWN_ANSWER = 0;
const int THREADS_PER_BLOCK = 512;
 // should be a multiple of block size and POPULATION_SIZE>THREADS_PER_BLOCK
const unsigned U_RAND_MAX = static_cast<unsigned>(RAND_MAX) + 1;
const int VAR_NUMBER = 5;
float best_population[VAR_NUMBER];
float best_score = 1000;

//////const char *pStrPipeName = "\\\\.\\pipe\\Name_pipe_demon_get";
//////const int BUFFER_MAX_LEN = 1024;
//////char buf[BUFFER_MAX_LEN];
// random number in [0, 1)
map<int, float> index_mScore;
bool cmp(const pair<int, float> &p1, const pair<int, float> &p2)//要用常数，不然编译错误   
{
	return p1.second<p2.second;
}
//vector<pair<int, int>> vtMap;
double  duration;
extern vector<vector<float>> mdata_g;
extern vector<vector<float>> mdata_score;
extern vector<float> mdata_time;
extern int MFlage;
extern int m_gen;
extern double  durationMy;
extern bool myexit;
extern std::mutex mtx;
void convert_data_two(vector<float> tmp, float* data)
{
	for (int i = 0; i < tmp.size(); i++)
	{
		data[i] = tmp[i];
	}
}

float float_random(M_args_Bound Parameter_temp) {
	return Parameter_temp.g[0] + (Parameter_temp.g[1] - Parameter_temp.g[0])*((static_cast<float>(rand())) / (U_RAND_MAX));
}

struct ScoreWithId {
	float score;
	int id;
};

//struct HH_Param{
//	float Na_0 = 0;
//	float K_0 = 0;
//	float Ca_0 = 0;
//
//};
struct M_args_Tset
{
	float *spike_TestData;
	int length;
};



//function HH_return for test
void HH_return(float *List_param, int VAR_NUMBER, double Mtime, double tempVB, double TimeStep, vector<float> m_I, int FlagParameter[], float gL, float C, float *Array_Data, int POPULATION_SIZE)
{
	//float *Array_Data = new float[int(Mtime / TimeStep)];
	//HMODULE hDll = LoadLibrary("CHH_dll.dll");
	
	if (myexit) return;
	if (hDll != NULL)
	{
		AddFunc HH_Entrance = (AddFunc)GetProcAddress(hDll, "HH_Entrance");
		if (HH_Entrance != NULL)
		{
			//printf("aa %f %f\n", List_param[0], List_param[1]);
			//HH_Entrance(120, 36, 100, -69, 0.01, 10);
			float *pI = new float[m_I.size()];
			for (int i = 0; i < m_I.size(); i++)
				pI[i] = m_I[i];
			HH_Entrance(List_param, Mtime, tempVB, TimeStep, pI, FlagParameter, gL, C, Array_Data, POPULATION_SIZE);
			//Array_Data = (float*)HH_Entrance(120, 36, Mtime, tempVB, TimeStep, m_I);
			delete pI;
			
		}
		else
		{
			{
				printf("not find dll");
			}
		}
		
		FreeLibrary(hDll);
	}
	else
	{
		{
			printf("not find dll");
		}
	}
	for (int i = 0; i < int(Mtime) * 100; i++)
	{
		printf("%f\n", Array_Data[i]);
	}
	////////////////
}



vector<float> HH_SpikeTime(float *data, double Mtime, double TimeStep, int &length,int spike_length)
{
	int size = int(Mtime / TimeStep);
	int k = 0;
	vector<float>Mdata_copy;
	//float *Mdata_copy = new float[spike_length];
	for (int i = 1; i < size - 1; i++)
	{
		if (data[i]>data[i - 1] && data[i] > data[i + 1] && data[i] >= 0)
		{
			Mdata_copy.push_back(i*TimeStep);
			k++;
		}
	}
	
	if (k == 0)
	{
		k = 1;
		Mdata_copy.push_back(Mtime / TimeStep);
	}
	length = k;
	//for (int j = 0; j < spike_length; j++)
	//{
	//	Mdata_copy[j] = j*10000;
	//}
	//length=spike_length;
	return Mdata_copy;
}

float fitness_B(M_args deviceParameter, M_args_Tset deviceParameter_Tset, float tau) {
	float result = 0;
	//printf("%d_a",deviceParameter.spike_data_num);
	//printf("%d_b", deviceParameter_Tset.length);
	//#pragma omp parallel for reduction(+:result)
	for (int i = 0; i<deviceParameter_Tset.length; i++)
	for (int j = 0; j<deviceParameter_Tset.length; j++)
	{
		result += exp(-abs(deviceParameter_Tset.spike_TestData[i] - deviceParameter_Tset.spike_TestData[j])*1.0 / tau);
		//printf("%f_3 ", result);
		// ++curPos;
	}
	//#pragma omp parallel for reduction(+:result)
	for (int i = 0; i<deviceParameter.spike_data_num; i++)
	for (int j = 0; j<deviceParameter.spike_data_num; j++)
	{
		result += exp(-abs(deviceParameter.spike_data[i] - deviceParameter.spike_data[j])*1.0 / tau);
		//printf("%f_2 ", result);
		// ++curPos;
	}
	//#pragma omp parallel for reduction(-:result)
	for (int i = 0; i<deviceParameter.spike_data_num; i++)
	for (int j = 0; j<deviceParameter_Tset.length; j++)
	{
		//printf("%f_c ", deviceParameter.spike_data[i]);
		//printf("%f_d ", deviceParameter_Tset.spike_TestData[j]);
		result -= 2 * exp(-abs(deviceParameter.spike_data[i] - deviceParameter_Tset.spike_TestData[j])*1.0 / tau);
		//printf("%f_1 ", result);
		// ++curPos;
	}
	//printf("%f_4 ", result);
	//cout << result << endl;
	return result;
}


float FinallResult(int Generation, float tau, float *population, double Mtime, double tempVB, double TimeStep, vector<float> m_I, M_args &Parameter_, int FlagParameter[], float gL, float C, M_args_Bound Parameter_Bound[], const int POPULATION_SIZE, stringstream &strResult){
	ScoreWithId *score = new ScoreWithId[POPULATION_SIZE];
	float *population_tmp = new float[POPULATION_SIZE*VAR_NUMBER];
	M_args_Tset *Parameter_Tmp = new M_args_Tset[POPULATION_SIZE];
	float score_min = 1000;
	float score_max = 0;
	int index = 0;
	int max_index = 0;
	//mdata_score.clear();
	mdata_g.clear();
	vector<float>tmp_score;
	int size = int(Mtime / TimeStep);
	
	clock_t start, finish;
	string name1= "str_g_" + to_string(Generation)+".txt";
	string name2 = "str_score_" + to_string(Generation) + ".txt";
	string name3 = "str_5dim_" + to_string(Generation) + ".txt";
	const char * mystr = name1.c_str();
	const char * mystr_score = name2.c_str();
	const char * mystr_5dim = name3.c_str();
	ofstream outf;
	ofstream outf2;
	ofstream outf3;
	outf.open(mystr);
	outf2.open(mystr_score);
	outf3.open(mystr_5dim);
	//std::cout << Parameter_Bound[0].g[0] << ' ' << Parameter_Bound[0].g[1] << ' ' << Parameter_Bound[1].g[0] << ' ' << Parameter_Bound[1].g[1] << endl;
	vector<int> tmp_param;
	vector<int> max_id;
	vector<pair<int, float> > arr_pair;
	///#pragma omp parallel for num_threads(4)
	////for (int j = 0; j < POPULATION_SIZE; j++)
	////{
	////	vector<float> temp_mdata_g;
	////	vector<float> temp_mdata_param;
	////	for (int i = 0; i < VAR_NUMBER; i++)
	////	{

	////		if (!(population[j*VAR_NUMBER + i] >= Parameter_Bound[i].g[0] && population[j*VAR_NUMBER + i] <= Parameter_Bound[i].g[1]))
	////		{
	////			outf <<"Y"<<i<<": "<< population[j*VAR_NUMBER + i] << " ";
	////			temp_mdata_param.push_back(i);
	////			float mtp = population[index*VAR_NUMBER + i] + (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.3*((static_cast<float>(rand())) / (U_RAND_MAX)-0.5);
	////			mtp = mtp < Parameter_Bound[i].g[0] ? Parameter_Bound[i].g[0] + (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
	////			mtp = mtp > Parameter_Bound[i].g[1] ? Parameter_Bound[i].g[1] - (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
	////			population[j*VAR_NUMBER + i] = mtp;
	////		}
	////		if (i < 2)
	////		{
	////			outf << population[j*VAR_NUMBER + i] << " ";
	////			temp_mdata_g.push_back(population[j*VAR_NUMBER + i]);
	////		}
	////	}
	////	tmp_param.push_back(temp_mdata_param);
	////	outf << endl;
	////	mdata_g.push_back(temp_mdata_g);
	////	}
		//float *temp_spike_TestData;

		//Parameter_Tset[j].spike_TestData = HH_SpikeTime(temp_data, Mtime, TimeStep, Parameter_Tset[j].length, Parameter_.spike_data_num);

	if (best_score < 1)
	{
		int rd_id = POPULATION_SIZE*(static_cast<float>(rand())) / (U_RAND_MAX);
		for (int i = 0; i < VAR_NUMBER; i++)
			best_population[i] = population[rd_id*VAR_NUMBER + i];
	}

	for (int j = 0; j < POPULATION_SIZE; j++)
		{
			vector<float> temp_mdata_g;
			vector<int> temp_mdata_param;
			for (int i = 0; i < VAR_NUMBER; i++)
			{
				if (!(population[j*VAR_NUMBER + i] >= Parameter_Bound[i].g[0] && population[j*VAR_NUMBER + i] <= Parameter_Bound[i].g[1]))
				{
					outf << "Y" << i << ": " << population[j*VAR_NUMBER + i] << " ";
					//temp_mdata_param.push_back(i);
					//float mtp = population[index*VAR_NUMBER + i] + (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.3*((static_cast<float>(rand())) / (U_RAND_MAX)-0.5);
					//mtp = mtp < Parameter_Bound[i].g[0] ? Parameter_Bound[i].g[0] + (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
					//mtp = mtp > Parameter_Bound[i].g[1] ? Parameter_Bound[i].g[1] - (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
					float mtp = Parameter_Bound[i].g[0] + (Parameter_Bound[i].g[1] - Parameter_Bound[i].g[0])*((static_cast<float>(rand())) / (U_RAND_MAX));
					population[j*VAR_NUMBER + i] = mtp;
					tmp_param.push_back(j);
				}
				if (i < 2)
				{
					outf << population[j*VAR_NUMBER + i] << " ";
					temp_mdata_g.push_back(population[j*VAR_NUMBER + i]);
				}
			}
			//tmp_param.push_back(temp_mdata_param);
			outf << endl;
			mdata_g.push_back(temp_mdata_g);
		}
	    
		start = clock();
		hDll = LoadLibrary("CHdll.dll");
		int spike_length = int(Mtime / TimeStep);
		float *temp_data = new float[spike_length*POPULATION_SIZE];
		HH_return(population, VAR_NUMBER, Mtime, tempVB, TimeStep, m_I, FlagParameter, gL, C, temp_data, POPULATION_SIZE);
		finish = clock();
		FreeLibrary(hDll);
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "time4 " << duration << endl;
		for (int j = 0; j < POPULATION_SIZE; j++)
		{
			vector<float>mdata_tmp = HH_SpikeTime(&temp_data[j*spike_length], Mtime, TimeStep, Parameter_Tmp[j].length, Parameter_.spike_data_num);
			Parameter_Tmp[j].spike_TestData = new float[mdata_tmp.size()];
			convert_data_two(mdata_tmp, Parameter_Tmp[j].spike_TestData);
		////for (int i = 0; i < Parameter_.spike_data_num; i++)
		////	cout << Parameter_.spike_data[i] << " ";

			score[j].score = fitness_B(Parameter_, Parameter_Tmp[j], tau);
			outf2 << score[j].score << endl;
			tmp_score.push_back(score[j].score);
			index_mScore[j] = score[j].score;
			if (score[j].score < score_min)
			{
				score_min = score[j].score;
				index = j;
			}
			if (score[j].score > score_max)
			{
				score_max = score[j].score;
				max_index = j;
				max_id.push_back(max_index);
			}
		}
	mdata_score.push_back(tmp_score);
	delete temp_data;
	for (int j = 0; j < POPULATION_SIZE; j++)
	{
		int indexm;
		int ind_score = 1000;
		for (int i = 0; i <POPULATION_SIZE*0.002; i++)
		{
			int first = POPULATION_SIZE*(static_cast<float>(rand())) / (U_RAND_MAX);
			if (score[first].score < ind_score)
			{
				indexm = first;
				ind_score=	score[first].score ;
			}
		}
		
		for (int i = 0; i < VAR_NUMBER; i++)
			population_tmp[j*VAR_NUMBER + i] = population[indexm*VAR_NUMBER + i];
	}
	outf.close();
	outf2.close();
	strResult<<"Generation No:" <<to_string(Generation)<<"\r\n";
	strResult << "best result id is " << to_string(index) << "\r\n";
	cout<< "Generation No:" << to_string(Generation) << endl;
	cout<< "best result id is " << to_string(index) << endl;
	///printf("best result id is %d\n", index);
	////printf("Result: ");
	cout << "========================" << endl;
	cout << "Score(1),Parameter(2-)"<< endl;
	cout << "========================" << endl;
	cout << to_string(score[index].score) << " ";
	strResult << "========================" << "\r\n";
	strResult << "Score(1),Parameter(2-)" << "\r\n";
	strResult << "========================" << "\r\n";
	cout << "The best individual" << endl;
	strResult << "The best individual" << "\r\n";
	strResult << to_string(score[index].score) << " ";
	////for (int i = 0; i < tmp_param.size();i++)
	////	for (int temp_i = 0; temp_i < VAR_NUMBER; temp_i++)
	////	{
	////		float mtp = population[index*VAR_NUMBER + temp_i] + (Parameter_Bound[temp_i].g[1] - Parameter_Bound[temp_i].g[0])*0.3*((static_cast<float>(rand())) / (U_RAND_MAX)-0.5);
	////		mtp = mtp < Parameter_Bound[temp_i].g[0] ? Parameter_Bound[temp_i].g[0] + (Parameter_Bound[temp_i].g[1] - Parameter_Bound[temp_i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
	////		mtp = mtp > Parameter_Bound[temp_i].g[1] ? Parameter_Bound[temp_i].g[1] - (Parameter_Bound[temp_i].g[1] - Parameter_Bound[temp_i].g[0])*0.1*(static_cast<float>(rand())) / (U_RAND_MAX) : mtp;
	////		population[i*VAR_NUMBER + temp_i] = mtp;
	////}

	for (int j = 0; j < VAR_NUMBER; j++)
	{
			strResult << to_string(population[index*VAR_NUMBER + j]) << " ";
			cout << to_string(population[index*VAR_NUMBER + j]) << " ";
		}
	cout << endl;
	strResult << "\r\n"; 
	if (best_score>score[index].score)
	{
		best_score = score[index].score;
		for (int j = 0; j < VAR_NUMBER; j++)
			best_population[j] = population[index*VAR_NUMBER + j];
	}
	for (map<int, float>::iterator it = index_mScore.begin(); it != index_mScore.end(); ++it)
	{
		arr_pair.push_back(make_pair(it->first, it->second));
	}
	sort(arr_pair.begin(), arr_pair.end(), cmp);
	////int kl = 0;
	////for (vector<pair<int, float> >::iterator it = arr_pair.begin(), ite = arr_pair.end()-1; it != arr_pair.end(); ++it,--ite)
	////{
	////	for (int j = 0; j < VAR_NUMBER; j++)
	////	{
	////		population[ite->first * VAR_NUMBER + j] = population[it->first * VAR_NUMBER + j];
	////	}
	////	kl++;
	////	if (kl>POPULATION_SIZE*0.1) break;
	////}
	cout << endl;
	//cout << "All individuals" << endl;
	strResult << "All individuals" << "\r\n";
	#pragma omp parallel for
	for (int i = 0; i < THREADS_PER_BLOCK; i++)
	{
		//cout << to_string(score[i].score) << " ";
		strResult << to_string(score[i].score) << " ";
		for (int j = 0; j < VAR_NUMBER; j++)
		{
			//cout << to_string(population[i*VAR_NUMBER + j]) << " ";
			strResult << to_string(population[i*VAR_NUMBER + j]) << " ";
		}
		//cout << endl;
		strResult << "\r\n";
	}
	//cout << endl;
	strResult << "\r\n";

		////	printf("%f ", population[index*VAR_NUMBER+ j]);
	////printf("\n spike_time\n");
	////for (int j = 0; j < Parameter_Tmp[index].length; j++)
	////	printf("   %f \n", Parameter_Tmp[index].spike_TestData[j]);
	////float *temp_data = new float[int(Mtime / TimeStep)];
	////HH_return(&population[index], VAR_NUMBER, Mtime, tempVB, TimeStep, m_I, FlagParameter, gL, C, temp_data);
	//float aa[2] = { 129.229584, 32.623901 };
	//temp_data = HH_return(aa, VAR_NUMBER, Mtime, tempVB, TimeStep, m_I);
	////for (int j = 0; j < size; j++)
	////	printf("   %f \n", temp_data[j]);
	for (int i = 0; i < POPULATION_SIZE; i++)
	{
		for (int j = 0; j < VAR_NUMBER; j++)
		{
			population[i*VAR_NUMBER + j] = population_tmp[i*VAR_NUMBER + j];
			outf3 << population[i*VAR_NUMBER + j] << " ";
		}
		outf3 << endl;
	}
	outf3.close();
	delete population_tmp;
	for (int jk = 0; jk < POPULATION_SIZE; jk++)
		delete Parameter_Tmp[jk].spike_TestData;
	//delete temp_data;
	delete Parameter_Tmp;
	float finscore = score[index].score;
	delete score;
	
	return finscore;
}