#pragma once
#include "Ca.h"
class CCa
{
protected:
	double g_Ca;
	double m_c;
public:
	CCa(double gca, double n);
	double Compute_n(double tempV, double TimeStep);
	double Compute_KI(double tempV, double TimeStep);
	double Iterative(double tmp_n, double an, double bn, double TimeStep);
	CCa();
	~CCa();
};



















