#pragma once

#include "K.h"
#include "Na.h"
#include "Ca.h"
#include "KM.h"
#include "Kv.h"
class CCell
{
protected:
	double m_start;
	double m_Current;
	double m_TimeStep;
public:
	CNa CNa_a;
	CK CK_a;
	CKM CKM_a;
	CKv CKv_a;
	CCa CCa_a;
	float *m_pArray_V;

public:
	CCell(double start, double TimeStep, CNa Na_a, CK K_a, CKM KM_a, CKv Kv_a, CCa Ca_a);
	double CCell::Iterative(double tmp_n, double an, double bn, double TimeStep);

	double  CCell::GetStart();
	double  CCell::GetTimeStep();
	//CStringArray CCell::GetPArray_V();

	~CCell();
};

