#pragma once
#include "Na.h"


class CNa
{
protected:
	double g_Na;
	double m_Nam;
	double m_Nah;
public:
	CNa(double gna, double m, double h);
	double Compute_m(double tempV,double TimeStep);
	double Compute_h(double tempV,double TimeStep);
	double Compute_NaI(double tempV,double TimeStep);
	double Iterative(double tmp_n, double an, double bn, double TimeStep);
	CNa();
	~CNa();
};

