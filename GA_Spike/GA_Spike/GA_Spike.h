
// GA_Spike.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CGA_SpikeApp: 
// �йش����ʵ�֣������ GA_Spike.cpp
//

class CGA_SpikeApp : public CWinApp
{
public:
	CGA_SpikeApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CGA_SpikeApp theApp;