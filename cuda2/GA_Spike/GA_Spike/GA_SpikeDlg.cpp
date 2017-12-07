
// GA_SpikeDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "GA_Spike.h"
#include "GA_SpikeDlg.h"
#include "afxdialogex.h"
#include <mutex> 
using namespace std;
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
std::mutex mtx;
std::string s;
std::stringstream strResult(s);
int m_length = 0;
std::atomic<bool> update = FALSE;
threadInfo Info;
M_args Parameter_;
vector<float> mdata_score;
vector<float> mdata_time;
vector<vector<float>> mdata_g;
int MFlage = -1;
int m_gen = 0;
double  duration;
bool myexit;
CString FilePathName;
//////typedef struct SthData
//////{
//////	M_args_Bound Parameter_Bound[VAR_NUMBER];
//////	double Mtime;//ms
//////	double tempVB;//初值
//////	double m_I;//电流
//////	float gL;
//////	float C;//电容
//////	int MaxGeneration;
//////	float crossver;
//////	float mutations;
//////	double TimeStep;
//////	M_args Parameter_;
//////	const int POPULATION_SIZE;
//////	int FlagParameter[5];//Na,K,KM,Kv,Ca
//////
//////}*pSthData;
//CFileDialog fileDlg(TRUE);
CGA_SpikeDlg::~CGA_SpikeDlg()
{
	if (!myexit)
		mThread->join();
}

CGA_SpikeDlg::CGA_SpikeDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CGA_SpikeDlg::IDD, pParent)
	, m_Na_low(0)
	, m_Na_up(0)
	, m_K_low(0)
	, m_K_up(0)
	, m_Ca_low(0)
	, m_Ca_up(0)
	, m_KM_low(0)
	, m_KM_up(0)
	, m_Kv_low(0)
	, m_Kv_up(0)
	, m_MaxGeneration(0)
	, m_crossover(0)
	, m_POPULATION_SIZE(0)
	, m_mutations(0)
	, m_V(0)
	, m_C(0)
	, m_time(0)
	, m_gL(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CGA_SpikeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_Na_Low, m_Na_low);
	DDX_Text(pDX, IDC_EDIT_Na_Up, m_Na_up);
	DDX_Text(pDX, IDC_EDIT_K_Low, m_K_low);
	DDX_Text(pDX, IDC_EDIT_K_Up, m_K_up);
	DDX_Text(pDX, IDC_EDIT_Ca_Low, m_Ca_low);
	DDX_Text(pDX, IDC_EDIT_Ca_Up, m_Ca_up);
	DDX_Text(pDX, IDC_EDIT_KM_Low, m_KM_low);
	DDX_Text(pDX, IDC_EDIT_KM_Up, m_KM_up);
	DDX_Text(pDX, IDC_EDIT_Kv_Low, m_Kv_low);
	DDX_Text(pDX, IDC_EDIT_Kv_Up, m_Kv_up);
	DDX_Text(pDX, IDC_EDIT_MaxGeneration, m_MaxGeneration);
	DDX_Text(pDX, IDC_EDIT_Crossover, m_crossover);
	DDX_Text(pDX, IDC_EDIT_Population_size, m_POPULATION_SIZE);
	DDX_Text(pDX, IDC_EDIT_Mutations, m_mutations);
	DDX_Text(pDX, IDC_EDIT_V, m_V);
	DDX_Text(pDX, IDC_EDIT_Capacitance, m_C);
	DDX_Text(pDX, IDC_EDIT_Time, m_time);
	DDX_Text(pDX, IDC_EDIT_gL, m_gL);
	////DDX_Control(pDX, IDC_EDIT_Result, m_result);
}

BEGIN_MESSAGE_MAP(CGA_SpikeDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//ON_EN_CHANGE(IDC_EDIT4, &CGA_SpikeDlg::OnEnChangeEdit4)
	ON_BN_CLICKED(IDC_BUTTON_RUN, &CGA_SpikeDlg::OnBnClickedButtonRun)
	ON_BN_CLICKED(IDCANCEL, &CGA_SpikeDlg::OnBnClickedCancel)
	ON_WM_TIMER()
	ON_BN_CLICKED(IDC_BUTTON_OPEN, &CGA_SpikeDlg::OnBnClickedButtonOpen)
//	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()


// CGA_SpikeDlg 消息处理程序

BOOL CGA_SpikeDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();
	CString   str_NaLow = "100";
	CString   str_NaUp = "150";
	CString   str_Up = "10";
	CString   str_KMLow = "1550";
	CString   str_KMUp = "1600";
	CString   str_KvLow = "20";
	CString   str_KvUp = "50";
	CString   str_MaxGeneration = "3";
	CString   str_V = "-69";
	CString   str_gL = "0.3";
	CString   str_Capacitance = "1";
	CString   str_Population_size = "2048";
	CString   str_Crossover = "0.5";
	CString   str_Mutations = "0.5";
	CString   str_Time = "100";
	//GetDlgItem(IDC_EDIT_Inputcurrent)->SetWindowText(str_Current);
	GetDlgItem(IDC_EDIT_Na_Low)->SetWindowText(str_NaLow);
	GetDlgItem(IDC_EDIT_Na_Up)->SetWindowText(str_NaUp);
	GetDlgItem(IDC_EDIT_K_Up)->SetWindowText(str_Up);
	GetDlgItem(IDC_EDIT_Ca_Up)->SetWindowText(str_Up);
	GetDlgItem(IDC_EDIT_KM_Low)->SetWindowText(str_KMLow);
	GetDlgItem(IDC_EDIT_KM_Up)->SetWindowText(str_KMUp);
	GetDlgItem(IDC_EDIT_Kv_Low)->SetWindowText(str_KvLow);
	GetDlgItem(IDC_EDIT_Kv_Up)->SetWindowText(str_KvUp);
	GetDlgItem(IDC_EDIT_Crossover)->SetWindowText(str_Crossover);
	GetDlgItem(IDC_EDIT_Population_size)->SetWindowText(str_Population_size);
	GetDlgItem(IDC_EDIT_Mutations)->SetWindowText(str_Mutations);
	GetDlgItem(IDC_EDIT_V)->SetWindowText(str_V);
	GetDlgItem(IDC_EDIT_Capacitance)->SetWindowText(str_Capacitance);
	GetDlgItem(IDC_EDIT_gL)->SetWindowText(str_gL);
	GetDlgItem(IDC_EDIT_Time)->SetWindowText(str_Time);
	GetDlgItem(IDC_EDIT_MaxGeneration)->SetWindowText(str_MaxGeneration);
	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO:  在此添加额外的初始化代码
	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CGA_SpikeDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
		/////////////////////////
		/////////////////////////
		////CDC *pDc = GetDlgItem(IDC_STATIC_PictureB)->GetDC();
		////GetDlgItem(IDC_STATIC_PictureB)->GetClientRect(&rect);
		////CDC memDC;
		////int w = rect.right - rect.left;
		////int h = rect.bottom - rect.top;
		////memDC.CreateCompatibleDC(pDc);//创建与目标DC相兼容的内存DC，
		////memDC.SelectObject(&memBitmapB);
		////pDc->BitBlt(0, 0, w, h, &memDC, 0, 0, SRCCOPY);//将内存DC上的图象拷贝到前台
		////memDC.DeleteDC();

		///////////////////////////
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CGA_SpikeDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

UINT postThread(LPVOID pParam)
{
	CGA_SpikeDlg*p = (CGA_SpikeDlg*)pParam;
	p->RunTimer();
	return 0;
}

void CGA_SpikeDlg::threadNew()
{
	AfxBeginThread(postThread, this);
}

////void CGA_SpikeDlg::threadNew2(threadInfo &Info)
////{
////	AfxBeginThread(ThreadFunc, &Info);
////}
////UINT ThreadFunc(LPVOID lpParam)
////{
////	threadInfo* pInfo = (threadInfo*)lpParam;
////	solveGPU_cpp(pInfo->Parameter_, pInfo->Mtime, pInfo->tempVB, pInfo->TimeStep, pInfo->m_I, pInfo->FlagParameter, pInfo->Parameter_Bound, pInfo->MaxGeneration, pInfo->gL, pInfo->C, pInfo->POPULATION_SIZE, pInfo->crossver, pInfo->mutations, strResult);
////
////}

void CGA_SpikeDlg::RunTimer()
{
	SetTimer(1, 100, NULL);
}

void task1(threadInfo &Info)
{
	solveGPU_cpp(Info.Parameter_, Info.Mtime, Info.tempVB, Info.TimeStep, Info.m_I, Info.FlagParameter, Info.Parameter_Bound, Info.MaxGeneration, Info.gL, Info.C, Info.POPULATION_SIZE, Info.crossver, Info.mutations, strResult);
	
	cout << "task1 says: " << endl;
}
////void task1(M_args &Parameter_, double Mtime, double tempVB, double TimeStep, double m_I, int FlagParameter[], M_args_Bound Parameter_Bound[], int MaxGeneration, float gL, float C, const int POPULATION_SIZE, float crossver, float mutations, stringstream &strResult)
////{
////	solveGPU_cpp(Parameter_, Mtime, tempVB, TimeStep, m_I, FlagParameter, Parameter_Bound, MaxGeneration, gL, C, POPULATION_SIZE, crossver, mutations, strResult);
////
////	//cout << "task1 says: " << endl;
////}

void CGA_SpikeDlg::OnBnClickedButtonRun()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(true);
	
	threadNew();
	
	/////////////////////////////
	//SetTimer(1, 1000, NULL);
	//////const int VAR_NUMBER = 5;
	const int VAR_NUMBER = 5;
	M_args_Bound Parameter_Bound[VAR_NUMBER];
	Parameter_Bound[0].g[0] = m_Na_low;
	Parameter_Bound[0].g[1] = m_Na_up;

	Parameter_Bound[1].g[0] = m_K_low;
	Parameter_Bound[1].g[1] = m_K_up;

	Parameter_Bound[2].g[0] = m_Ca_low;
	Parameter_Bound[2].g[1] = m_Ca_up;

	Parameter_Bound[3].g[0] = m_KM_low;
	Parameter_Bound[3].g[1] = m_KM_up;

	Parameter_Bound[4].g[0] = m_Kv_low;
	Parameter_Bound[4].g[1] = m_Kv_up;
	

	const int POPULATION_SIZE = m_POPULATION_SIZE;
	
	freopen("output.txt", "w", stdout);
	srand(1900);
	srand(static_cast<unsigned>(time(0)));
	//float *spike_data, *current_data, *spike_TestData;
	//float gNa = 120;
	//float gK = 36;
	double Mtime = m_time;//ms
	double tempVB = m_V;//初值
	double TimeStep = 0.01;
	double m_I = 10;//电流
	float gL = m_gL;
	float C = m_C;//电容
	int FlagParameter[5] = { 1, 1, 1, 1, 1 };//Na,K,KM,Kv,Ca

	//Parameter_Bound.length = 2;
	int MaxGeneration = m_MaxGeneration;
	float crossver = m_crossover;
	float mutations = m_mutations;
	std::string strFilePathName;
	strFilePathName = FilePathName.GetBuffer(0);
	vector<float>Mtemp_data = Read_Txt(strFilePathName, Parameter_.spike_data_num);
	Parameter_.spike_data = new float[Mtemp_data.size()];
	convert_data(Mtemp_data, Parameter_.spike_data);
	//float ans = solveGPU_cpp(Parameter_, Mtime, tempVB, TimeStep, m_I, FlagParameter, Parameter_Bound, MaxGeneration, gL, C, POPULATION_SIZE, crossver, mutations, strResult);
	//////////////////////////////////////////
	
	//////////////////////////////////////////
	//////pSthData pDataValue1 = new SthData();
	//////memset(pDataValue1, 0x00, sizeof(SthData));
	Info.Parameter_ = Parameter_;
	Info.Mtime = m_time;
	Info.tempVB = m_V;
	Info.m_I = 10;
	Info.FlagParameter[0] = FlagParameter[0];
	Info.FlagParameter[1] = FlagParameter[1];
	Info.FlagParameter[2] = FlagParameter[2];
	Info.FlagParameter[3] = FlagParameter[3];
	Info.FlagParameter[4] = FlagParameter[4];
	Info.Parameter_Bound[0] = Parameter_Bound[0];
	Info.Parameter_Bound[1] = Parameter_Bound[1];
	Info.Parameter_Bound[2] = Parameter_Bound[2];
	Info.Parameter_Bound[3] = Parameter_Bound[3];
	Info.Parameter_Bound[4] = Parameter_Bound[4];
	Info.MaxGeneration = m_MaxGeneration;
	Info.gL = m_gL;
	Info.C = m_C;
	Info.POPULATION_SIZE = m_POPULATION_SIZE;
	Info.crossver = m_crossover;
	Info.mutations = m_mutations;
	Info.TimeStep = 0.01;
	//////threadNew2(Info);
	update = TRUE;
	/////////////////////////////////
	///mThread = std::make_shared<std::thread>(task1, Parameter_, Mtime, tempVB, TimeStep, m_I, FlagParameter, Parameter_Bound, MaxGeneration, gL, C, POPULATION_SIZE, crossver, mutations, std::ref(strResult));
	mThread = std::make_shared<std::thread>(task1, Info);
	////thread t1(task1, Parameter_, m_time, m_V, TimeStep, m_I, FlagParameter, Parameter_Bound, m_MaxGeneration, m_gL, m_C, m_POPULATION_SIZE, m_crossover, m_mutations, std::ref(strResult));
	//t1.join();
	//Sleep(1000);
	////t1.join();
	//mThread->detach();
	//mThread->join();
	//delete Parameter_.spike_data;
	//std::cout << "GPU answer = " << ans << std::endl;
	
	//printf("%f seconds\n", duration);
	//KillTimer(1);
}


void CGA_SpikeDlg::OnBnClickedCancel()
{
	// TODO:  在此添加控件通知处理程序代码
	delete Parameter_.spike_data;
	KillTimer(1);
	myexit = true;
	if (update)
		mThread->join();
	CDialogEx::OnCancel();
	
}

void CGA_SpikeDlg::DrawPictureB(vector<vector<float>> mdata_g, vector<float>mdata_score)
{   ///////////////gNa,gK
	//Invalidate();
	CPaintDC dc(this);
	RECT rect;
	CDC memDC;
	CDC *pDc = GetDlgItem(IDC_STATIC_PictureB)->GetDC();
	GetDlgItem(IDC_STATIC_PictureB)->GetClientRect(&rect);
	InvalidateRect(&rect);
	memDC.CreateCompatibleDC(pDc);//创建与目标DC相兼容的内存DC，
	CBitmap memBitmapB;
	int w = rect.right - rect.left;
	int h = rect.bottom - rect.top;
	//创建一个内存中的图像
	///CBitmap memBitmapB;
	memBitmapB.CreateCompatibleBitmap(pDc, w, h);//根据目标DC创建位图
	memDC.SelectObject(&memBitmapB);//把位图选入内存DC
	memDC.FillSolidRect(rect.left,rect.top,rect.right,rect.bottom, pDc->GetBkColor());//按原来背景填充客户区，不然会是黑色
	memDC.MoveTo(30, h - 30);
	memDC.LineTo(w - 30, h - 30);
	memDC.MoveTo(30, h - 30);
	memDC.LineTo(30, 30);
	memDC.MoveTo(w - 35, h - 25);
	memDC.LineTo(w - 30, h - 30);
	memDC.LineTo(w - 35, h - 35);
	memDC.MoveTo(25, 35);
	memDC.LineTo(30, 30);
	memDC.LineTo(35, 35);
	//////x轴刻度
	//////////y轴刻度
	double dividx = (m_Na_up - m_Na_low) / 10;
	double dividy = (m_K_up - m_K_low) / 10;

	for (int j = 1; j <= 10; j++)
	{
		CString strY;
		CString strX;
		strY.Format("%g", m_K_low + j * dividy);
		strX.Format("%g", m_Na_low + j * dividx);
		memDC.TextOut(23 + j*1.0 * (w - 60) / 10, h - 24, strX);
		memDC.TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
	}
	//pDc->TextOut(5, h - 24, 0);
	//////////////////
	CBrush *oldbrush;
	CBrush brush;
	brush.CreateSolidBrush(RGB(255, 0, 0));
	//brush.CreateSolidBrush(RGB(255, 255, 0));
	oldbrush = memDC.SelectObject(&brush);
	for (int j = 0; j < mdata_g.size(); j++)////描点
	{
		float nscore = (mdata_score[j] * 200 / 20>255) ? 255 : (mdata_score[j] * 200 / 20);
		brush.DeleteObject();
		brush.CreateSolidBrush(RGB(nscore, nscore ,nscore));
		oldbrush = memDC.SelectObject(&brush);
		float m_x = mdata_g[j][0];
		float m_y = mdata_g[j][1];
		double x = 30 + (m_x - m_Na_low)*(w - 60) / (m_Na_up - m_Na_low);
		double y = h - 30 - (m_y - m_K_low)*(h - 60) / (m_K_up - m_K_low);
		memDC.Ellipse(x - 2, y - 2, x + 2, y + 2);
	}
	memDC.SelectObject(oldbrush);
	pDc->BitBlt(0, 0, w, h, &memDC, 0, 0, SRCCOPY);//将内存DC上的图象拷贝到前台
	memDC.DeleteDC();                                       //删除DC
	///memBitmapB.DeleteObject();                                        //删除位图
	ReleaseDC(pDc);
}
void CGA_SpikeDlg::DrawPictureA(vector<float>mdata_score)
{
	CPaintDC dc(this);
	RECT rect;
	///////////////score
	CDC *pDcB = GetDlgItem(IDC_STATIC_PictureA)->GetDC();
	GetDlgItem(IDC_STATIC_PictureA)->GetClientRect(&rect);
	int w = rect.right - rect.left;
	int h = rect.bottom - rect.top;
	pDcB->MoveTo(30, h - 30);
	pDcB->LineTo(w - 30, h - 30);
	pDcB->MoveTo(30, h - 30);
	pDcB->LineTo(30, 30);
	pDcB->MoveTo(w - 35, h - 25);
	pDcB->LineTo(w - 30, h - 30);
	pDcB->LineTo(w - 35, h - 35);
	pDcB->MoveTo(25, 35);
	pDcB->LineTo(30, 30);
	pDcB->LineTo(35, 35);
	//////
	//////////y轴刻度
	double dividx = (m_MaxGeneration - 0) / 10;
	double dividy = (50 - 0) / 10;

	for (int j = 1; j <= 10; j++)
	{
		CString strY;
		strY.Format("%g", 0 + j * dividy);
		pDcB->TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
	}
	////
	//x轴刻度
	for (int j = 1; j <= m_MaxGeneration; j++)
	{
		CString strX;
		strX.Format("%d", j);
		pDcB->TextOut(23 + j*1.0 * (w - 60) / m_MaxGeneration, h - 24, strX);
	}
	//pDc->TextOut(5, h - 24, 0);
	//////////////////
	for (int j = 0; j <mdata_score.size(); j++)////描点
	{
		float m_x = m_gen;
		float m_y = mdata_score[j];
		double x = 30 + (m_x - 0)*(w - 60) / (m_MaxGeneration - 0);
		double y = h - 30 - (m_y - 0)*(h - 60) / (50 - 0);
		pDcB->Ellipse(x - 2, y - 2, x + 2, y + 2);
	}
	//////////////////
	CBrush *oldbrush;
	CBrush brush;
	brush.CreateSolidBrush(RGB(255, 0, 0));
	oldbrush = pDcB->SelectObject(&brush);
	ReleaseDC(pDcB);
}

void CGA_SpikeDlg::DrawPictureC(double duration)
{
	CPaintDC dc(this);
	RECT rect;
	///////////////score
	CDC *pDcB = GetDlgItem(IDC_STATIC_PictureC)->GetDC();
	GetDlgItem(IDC_STATIC_PictureC)->GetClientRect(&rect);
	int w = rect.right - rect.left;
	int h = rect.bottom - rect.top;
	CBitmap memBitmap;
	pDcB->MoveTo(30, h - 30);
	pDcB->LineTo(w - 30, h - 30);
	pDcB->MoveTo(30, h - 30);
	pDcB->LineTo(30, 30);
	pDcB->MoveTo(w - 35, h - 25);
	pDcB->LineTo(w - 30, h - 30);
	pDcB->LineTo(w - 35, h - 35);
	pDcB->MoveTo(25, 35);
	pDcB->LineTo(30, 30);
	pDcB->LineTo(35, 35);
	//////
	//////////y轴刻度
	double dividx = (m_MaxGeneration - 0) / 10;
	double dividy = (500 - 0) / 10;

	for (int j = 1; j <= 10; j++)
	{
		CString strY;
		strY.Format("%g", 0 + j * dividy);
		pDcB->TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
	}
	////
	//x轴刻度
	for (int j = 1; j <= m_MaxGeneration; j++)
	{
		CString strX;
		strX.Format("%d", j);
		pDcB->TextOut(23 + j*1.0 * (w - 60) / m_MaxGeneration, h - 24, strX);
	}
	//pDc->TextOut(5, h - 24, 0);
	//////////////////
	float m_x = m_gen;
	float m_y = duration;
	double x = 30 + (m_x - 0)*(w - 60) / (m_MaxGeneration - 0);
	double y = h - 30 - (m_y - 0)*(h - 60) / (500 - 0);
	pDcB->Ellipse(x - 2, y - 2, x + 2, y + 2);
	//////////////////
	//CBrush *oldbrush;
	CBrush brush;
	brush.CreateSolidBrush(RGB(255, 0, 0));
	ReleaseDC(pDcB);
}

void CGA_SpikeDlg::OnTimer(UINT_PTR nIDEvent)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值
	switch (nIDEvent)
	{
		case 1:   //定时器1处理函数，定时发送数据进行更新
		{	
		mtx.lock();
		{
			//MFlage = 2;
			//cout << m_gen << endl;
			DrawPictureB(mdata_g, mdata_score);
			DrawPictureA(mdata_score);
			DrawPictureC(duration);
			//MFlage = 4;
		}
		mtx.unlock();
		}
	}
	CDialogEx::OnTimer(nIDEvent);
}


void CGA_SpikeDlg::OnBnClickedButtonOpen()
{
	// TODO:  在此添加控件通知处理程序代码
	
	CFileDialog fDlg(TRUE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT);
	fDlg.m_ofn.lpstrTitle = "Open file";
	fDlg.m_ofn.lpstrFilter = "Text Files(*.txt)\0*.txt\0All Files(*.*)\0*.*\0\0";
	if (fDlg.DoModal() == IDOK)
	{
		FilePathName = fDlg.GetPathName();//获取要打开文件的路径
		SetDlgItemText(IDC_EDIT_FilePath, FilePathName);
	}
	else
	{
		return;
	}

	////CFile f;
	////CString str_;
	////f.Open(FilePathName, CFile::modeCreate|CFile::modeReadWrite);
	////f.Read(str_.GetBuffer(f.GetLength()), f.GetLength());
	////f.Close();
	////CFileDialog dlg2(TRUE, //TRUE为OPEN对话框，FALSE为SAVE AS对话框
	////	NULL,
	////	FilePathName,
	////	OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
	////	(LPCTSTR)_TEXT("文本文件(*.txt,*.ini,*.log)|*.txt;*.ini;*.log|全部文件(*.*)|*.*||"),
	////	NULL);
}


//BOOL CGA_SpikeDlg::OnEraseBkgnd(CDC* pDC)
//{
//	// TODO:  在此添加消息处理程序代码和/或调用默认值
//
//	//return CDialogEx::OnEraseBkgnd(pDC);
//	return TRUE;
//}
