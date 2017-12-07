#include "stdafx.h"
#include "Cell.h"
#include "K.h"
#include "Na.h"

CCell::CCell(double start, double Current, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a)
{
	m_start = start;
	CNa_a = Na_a;
	CK_a = K_a;
	CKM_a = KM_a;
	CKv_a = Kv_a;
	CCa_a = Ca_a;
	//m_pArray_V = pArray_V;
	m_Current = Current;
	m_TimeStep = TimeStep;
}
double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep)
{
	double y_pre = tmp_n + (an*(1 - tmp_n) - bn*tmp_n)*TimeStep;
	double y_next = tmp_n + (an*(1 - y_pre) - bn*y_pre)*TimeStep;;
	double tmp_avgn = (y_pre + y_next)*1.0 / 2;
	return tmp_avgn;
}
double  CCell::GetStart()
{
	double start = m_start;
	return start;
}

double  CCell::GetTimeStep()
{
	double 	 TimeStep = m_TimeStep;
	return TimeStep;
}
double  CCell::GetCurrent()
{
	double 	 Current = m_Current;
	return Current;
}


//CSynapse  CCell::GetSyn_a()
//{
//	CSynapse syn_a =Csyn_a;
//	return syn_a;
//}
//CSynapse  CCell::GetSyn_b()
//{
//	CSynapse syn_b = Csyn_b;
//	return syn_b;
//}
//CSynapse  CCell::GetSyn_c()
//{
//	CSynapse syn_c = Csyn_c;
//	return syn_c;
//}
//CNMDA  CCell::GetNmda_a()
//{
//	CNMDA Nmda_a =CNmda_a;
//	return Nmda_a;
//}
//CStringArray CCell::GetPArray_V()
//{
//	CStringArray pArray_V= m_pArray_V;
//	return pArray_V;
//}

CCell::~CCell()
{
}
