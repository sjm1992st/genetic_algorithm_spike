#include "math.h"
#include <iostream>
#include "struct_cu.h"
#include<vector>
#include <sstream>
#include<queue> 
#include<vector> 
using namespace std;
vector<float> Read_Txt(string filename);
void convert_data(vector<float> tmp, float* data);
float solveGPU_cpp(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, vector<float> m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult);
