#include<vector>
using namespace std;
struct M_args
{
	float *spike_data;
	//float *current_data;
	//float *spike_TestData;
	float dt;
	//int length;
	int spike_data_num;
	//int current_data_num;
};


struct M_args_Bound
{
	float g[2];
};

struct threadInfo
{
	struct M_args_Bound Parameter_Bound[5];
	double Mtime;//ms
	double tempVB;//��ֵ
	vector<float> m_I;//����
	float gL;
	float C;//����
	int MaxGeneration;
	float crossver;
	float mutations;
	double TimeStep;
	struct M_args Parameter_;
	int POPULATION_SIZE;
	int FlagParameter[5];//Na,K,KM,Kv,Ca
};

