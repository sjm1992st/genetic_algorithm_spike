
// GA_SpikeDlg.h : ͷ�ļ�
//

#pragma once
#include "afxwin.h"
#include "cuda_transfer.h"
#include <thread>
#include <sstream>
#include <stdio.h>
#include <conio.h>  
#include <tchar.h>
#include <Windows.h>
#include <process.h>
#include <stdlib.h>
#include <memory>
#include <atomic> 

// CGA_SpikeDlg �Ի���



UINT ThreadFunc(LPVOID lpParam);
class CGA_SpikeDlg : public CDialogEx
{
// ����
public:
	CGA_SpikeDlg(CWnd* pParent = NULL);	// ��׼���캯��
// �Ի�������
	virtual ~CGA_SpikeDlg();
	enum { IDD = IDD_GA_SPIKE_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��


// ʵ��
protected:
	HICON m_hIcon;
	CWinThread* pThread;
	void CGA_SpikeDlg::DrawPictureB(vector<vector<float>> mdata_g, vector<float>mdata_score);
	void CGA_SpikeDlg::DrawPictureA(vector<float>mdata_score);
	void CGA_SpikeDlg::DrawPictureC(double duration);
	// ���ɵ���Ϣӳ�亯��
	
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeEdit4();
	float m_Na_low;
	float m_Na_up;
	float m_K_low;
	float m_K_up;
	float m_Ca_low;
	float m_Ca_up;
	float m_KM_low;
	float m_KM_up;
	float m_Kv_low;
	float m_Kv_up;
	int m_MaxGeneration;
	float m_crossover;
	int m_POPULATION_SIZE;
	float m_mutations;
	double m_V;
	float m_C;
	double m_time;
	float m_gL;
	
	afx_msg void OnBnClickedButtonRun();
	////CEdit m_result;
	void RunTimer();
	void threadNew();//***�����½���
	std::shared_ptr<std::thread>mThread;
	void threadNew2(threadInfo &Info);
	afx_msg void OnBnClickedCancel();

	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnBnClickedButtonOpen();
//	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
};
