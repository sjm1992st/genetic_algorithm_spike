#pragma once
#include "KM.h"
class CKM
{
protected:
	double g_KM;
	double m_z;
public:
	CKM(double gk, double n);
	double Compute_n(double tempV, double TimeStep);
	double Compute_KI(double tempV, double TimeStep);
	double Iterative(double tmp_n, double an, double bn, double TimeStep);
	CKM();
	~CKM();
};
