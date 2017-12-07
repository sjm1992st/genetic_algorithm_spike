#pragma once
#include "K.h"
class CK
{
protected:
	double g_K;
	double m_Kn;
public:
CK(double gk,double n);
double Compute_n(double tempV,double TimeStep);
double Compute_KI(double tempV,double TimeStep);
double Iterative(double tmp_n, double an, double bn, double TimeStep);
CK();
	~CK();
};

