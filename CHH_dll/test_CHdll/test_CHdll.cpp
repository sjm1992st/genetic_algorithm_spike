// test_CHdll.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <windows.h>
typedef void(*AddFunc)(float gNa, float gK, float gKM, float gKv, float gCa, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA);

int main(int argc, char *argv[])
{
	
	double Mtime = 200;
	double TimeStep = 0.01;
	int FlagParameter[5] = {1,1,1,1,1};
	float gL = 0.3;
	float C = 1;
	float *pArrayA = new float[int(Mtime / TimeStep)];
	HMODULE hDll = LoadLibrary("CHH_dll.dll");
	if (hDll != NULL)
	{
		AddFunc HH_Entrance = (AddFunc)GetProcAddress(hDll, "HH_Entrance");
		if (HH_Entrance != NULL)
		{
			HH_Entrance(122.699539, 4.984550, 3.272842, 1579.503906, 28.251514, Mtime, -69, TimeStep, 10, FlagParameter, gL, C, pArrayA);
		}
		FreeLibrary(hDll);
	}
	else
	{
		{
			printf("not find dll");
		}
	}
	freopen("V2_output.txt", "w", stdout);
	for (int i = 0; i < int(Mtime) * 100; i++)
	{
		printf("%f\n", pArrayA[i]);
	}
}
