#include "stdafx.h"
#include "K.h"
#include "Na.h"
CK::CK(double gk,double n)
{
	g_K = gk;
	m_Kn=n;
}
double CK::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double CK::Compute_n(double tempV,double TimeStep)
{
		   double an = 0.01*(tempV + 55) / (1 - exp(-0.1*(tempV + 55)));
		   double bn = 0.125*exp(-(tempV + 65) *1.0 / 80);
		   m_Kn = Iterative(m_Kn, an, bn,TimeStep);
		   return m_Kn;
}
double CK::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV,TimeStep);
	double K_I = g_K * pow(tmp_n, 4)*(tempV+77);
	return K_I;
}
CK::CK()
{}
CK::~CK()
{
}
