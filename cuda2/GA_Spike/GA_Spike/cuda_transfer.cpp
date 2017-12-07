#include "stdafx.h"
#include "cuda_transfer.h"
#include <fstream> 
#include <sstream>

extern "C" float solveGPU(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult);
//模板函数：将string类型变量转换为常用的数值类型 
template <class Type>
Type stringToNum(const string& str){
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}

template <class T>
int getArrayLen(T& array)
{
	return (sizeof(array) / sizeof(array[0]) - 1);
}
vector<float> Read_Txt(string filename, int &num)
{
	//float *Mdata = new float[100000];
	vector<float> Mdata;
	ifstream in(filename);
	string line;
	int i = 0;
	if (in) // 有该文件  
	{
		while (getline(in, line)) // line中不包括每行的换行符  
		{
			//cout << stringToNum<float>(line)+0.015 << endl;
			Mdata.push_back(stringToNum<float>(line));
			i++;
		}
	}
	else // 没有该文件  
	{
		cout << "no such file" << endl;
		return Mdata;
	}
	num = i;
	return Mdata;
}

void convert_data(vector<float> tmp, float* data)
{
	for (int i = 0; i < tmp.size(); i++)
		data[i] = tmp[i];
}

float solveGPU_cpp(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult)
{
	return solveGPU(Parameter_, Mtime, tempVB, TimeStep, m_I, FlagParameter, Parameter_Bound, MaxGeneration, gL, C, POPULATION_SIZE, crossver, mutations, strResult);
}
