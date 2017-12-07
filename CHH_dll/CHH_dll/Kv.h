#pragma once
#include "Kv.h"
class CKv
{
protected:
	double g_Kv;
	double m_p;
public:
	CKv(double gk, double n);
	double Compute_n(double tempV, double TimeStep);
	double Compute_KI(double tempV, double TimeStep);
	double Iterative(double tmp_n, double an, double bn, double TimeStep);
	CKv();
	~CKv();
};
