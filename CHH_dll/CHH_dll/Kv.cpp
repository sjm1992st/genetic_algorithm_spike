#include "stdafx.h"
#include "Kv.h"
CKv::CKv(double gk, double n)
{
	g_Kv = gk;
	m_p = n;
}
double CKv::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double CKv::Compute_n(double tempV, double TimeStep)
{
	double an = (tempV - 95)*1.0 / (1 - exp(-(tempV - 95) / 11.8));
	double bn = 0.025*exp(-tempV / 22.222);
	m_p = Iterative(m_p, an, bn, TimeStep);
	return m_p;
}
double CKv::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double K_I = g_Kv * pow(tmp_n, 2)*(tempV + 77);
	return K_I;
}
CKv::CKv()
{}
CKv::~CKv()
{
}
