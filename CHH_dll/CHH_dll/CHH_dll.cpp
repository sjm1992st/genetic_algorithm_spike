// CHH_dll.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Cell.h"

/////////////////////////
#ifndef _MYCODE_H_
#define _MYCODE_H_
#ifdef DLLDEMO1_EXPORTS
#define EXPORTS_DEMO _declspec( dllexport )
#else
#define EXPORTS_DEMO _declspec(dllimport)
#endif
extern "C" EXPORTS_DEMO void HH_Entrance(float gNa, float gK, float gKM, float gKv, float gCa, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA);
#endif
/////////////////////////////////////////////////////
void Calculation_I(CCell Cell_B, double Mtime, int FlagParameter[], float gL, float C, float *pArrayA)
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
			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I)/C;

		cell_B_K_I = Cell_B.CK_a.Compute_KI(y_preB, Cell_B.GetTimeStep());

		cell_B_Na_I = Cell_B.CNa_a.Compute_NaI(y_preB, Cell_B.GetTimeStep());
		cell_B_Ca_I = Cell_B.CCa_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
		cell_B_KM_I = Cell_B.CKM_a.Compute_KI(y_preB, Cell_B.GetTimeStep());
		cell_B_Kv_I = Cell_B.CKv_a.Compute_KI(y_preB, Cell_B.GetTimeStep());

		y_nextB = tempVB + Cell_B.GetTimeStep()*(gL*(Cell_B.GetStart() - y_preB) + Cell_B.GetCurrent()
			- cell_B_Na_I - cell_B_K_I - cell_B_Ca_I - cell_B_KM_I - cell_B_Kv_I)/C;

		temp = (y_preB + y_nextB)*1.0 / 2;
		pArrayA[k] = temp;
		tempVB = temp;
	}
}


void HH_Entrance(float gNa, float gK, float gKM, float gKv, float gCa, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], float gL, float C, float *pArrayA)
{
	//freopen("V_output.txt", "w", stdout);
	//double Mtime = 100;
	//double tempVB = -69;
	//double TimeStep = 0.01;
	CK K_a(gK, 0.318);//36
	CNa Na_a(gNa, 0.053, 0.596);///(120, 0.053, 0.596)
	CKM KM_a(gKM, 0);///(120, 0.053, 0.596)
	CKv Kv_a(gKv, 0);///(120, 0.053, 0.596)
	CCa Ca_a(gCa, 0);///(120, 0.053, 0.596)
	//double m_I = 10;
	CCell Cell_D(tempVB, m_I, TimeStep, Na_a, K_a, KM_a, Kv_a, Ca_a);
	Calculation_I(Cell_D, Mtime, FlagParameter, gL, C, pArrayA);
	//for (int i = 0; i < int(Mtime) * 100; i++)
	//{
	//	printf("%f\n", pArrayA[i]);
	//}
	////////////////
	//delete pArrayA;

}



