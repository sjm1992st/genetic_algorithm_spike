#include "stdafx.h"
#include "K.h"
#include "Na.h"
CNa::CNa(double gna, double m,double h)
{
	g_Na = gna;
	m_Nam = m;
	m_Nah = h;
}
double CNa::Compute_m(double tempV, double TimeStep)
{
	double am = 0.1*(tempV + 40) / (1 - exp(-0.1*(tempV + 40)));
	double bm = 4 * exp(-(tempV + 65) *1.0 / 18);//////////18
	m_Nam =Iterative(m_Nam, am, bm,TimeStep);
	return m_Nam;
}
double CNa::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double CNa::Compute_h(double tempV,double TimeStep)
{
	double ah = 0.07*exp(-(tempV + 65)*1.0 / 20);//////////20
	double bh = 1.0 / (exp(-0.1*(tempV + 35)) + 1);
	m_Nah = Iterative(m_Nah, ah, bh,TimeStep);
	return m_Nah;
}
double CNa::Compute_NaI(double tempV,double TimeStep)
{
	double tmp_m = Compute_m(tempV,TimeStep);
	double tmp_h = Compute_h(tempV,TimeStep);
	double Na_I = g_Na * pow(tmp_m, 3)*tmp_h*(tempV - 55);
	return Na_I;
}
CNa::CNa()
{}
CNa::~CNa()
{}
