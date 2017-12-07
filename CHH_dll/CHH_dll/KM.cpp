#include "stdafx.h"
#include "KM.h"
CKM::CKM(double gk, double n)
{
	g_KM = gk;
	m_z = n;
}
double CKM::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double CKM::Compute_n(double tempV, double TimeStep)
{
	double an = 1.0 / (1 + exp(-0.2*(tempV + 39)))/75;
	double bn = 1.0/75-an;
	m_z= Iterative(m_z, an, bn, TimeStep);
	return m_z;
}
double CKM::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double K_I = g_KM * pow(tmp_n, 1)*(tempV + 77);
	return K_I;
}
CKM::CKM()
{}
CKM::~CKM()
{
}
