#include "stdafx.h"
#include "stdafx.h"
#include "Ca.h"
CCa::CCa(double gca, double n)
{
	g_Ca= gca;
	m_c = n;
}
double CCa::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double CCa::Compute_n(double tempV, double TimeStep)
{
	double an = 0.3*(tempV + 13) / (1 - exp(-(tempV + 13)*1.0 / 10));
	double bn = 10.0 * exp(-(tempV + 38) *1.0/ 18);
	m_c = Iterative(m_c, an, bn, TimeStep);
	return m_c;
}
double CCa::Compute_KI(double tempV, double TimeStep)
{
	double tmp_n = Compute_n(tempV, TimeStep);
	double Ca_I = g_Ca * pow(tmp_n, 3)*(tempV -120);
	return Ca_I;
}
CCa::CCa()
{}
CCa::~CCa()
{
}
