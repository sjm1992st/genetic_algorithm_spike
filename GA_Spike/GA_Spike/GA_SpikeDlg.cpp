
// GA_SpikeDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "GA_Spike.h"
#include "GA_SpikeDlg.h"
#include "afxdialogex.h"
#include <mutex> 
#include <cuda_runtime.h>
#include <time.h> 

using namespace std;

typedef void(*AddFunc)(float gNa, float gK, float gKM, float gKv, float gCa, double Mtime, double tempVB, double TimeStep, vector<float> m_I, int FlagParameter[], float gL, float C, float *pArrayA);

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
vector<vector<float>> mdata_score;
vector<float> mdata_time;
vector<vector<float>> mdata_g;
int MFlage = -1;
int m_gen = 0;
int mtmp_gen = 0;
double  durationMy=0;
bool myexit = false;
CString FilePathName="";
CString FilePathName_I = "";
clock_t Mstart = clock();
int pos = 0;
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
	, m_current_start(0)
	, m_current_duration(0)
	, m_current_val(0)
	, m_realNa(0)
	, m_realK(0)
	, m_realCa(0)
	, m_realKM(0)
	, m_realKv(0)

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
	//	DDX_Control(pDX, IDC_SPLIT_Style, m_splitbtn);
	//DDX_Control(pDX, IDC_SPLIT_Style, m_splitbtn);
	DDX_Text(pDX, IDC_EDIT_Start, m_current_start);
	DDX_Text(pDX, IDC_EDIT_duration, m_current_duration);
	DDX_Text(pDX, IDC_EDIT_VALUE, m_current_val);
	//  DDX_Control(pDX, IDC_LIST_INPUT, m_List);
	DDX_Control(pDX, IDC_LIST_INPUT, m_wndList);
	DDX_Control(pDX, IDC_CHECK_Na, m_checkNa);
	DDX_Control(pDX, IDC_CHECK_K, m_checkK);
	//  DDX_Control(pDX, IDC_CHECK_Ca, m_checkKM);
	DDX_Control(pDX, IDC_CHECK_Ca, m_checkCa);
	DDX_Control(pDX, IDC_CHECK_KM, m_checkKM);
	DDX_Control(pDX, IDC_CHECK_Kv, m_checkKv);
	//  DDX_Control(pDX, IDC_EDIT_GNa, m_realNa);
	//  DDX_Control(pDX, IDC_EDIT_GK, m_realK);
	//  DDX_Control(pDX, IDC_EDIT_GCa, m_realCa);
	//  DDX_Control(pDX, IDC_EDIT_GKM, m_realKM);
	//  DDX_Control(pDX, IDC_EDIT_GKv, m_realKv);
	DDX_Text(pDX, IDC_EDIT_GNa, m_realNa);
	DDX_Text(pDX, IDC_EDIT_GK, m_realK);
	DDX_Text(pDX, IDC_EDIT_GCa, m_realCa);
	DDX_Text(pDX, IDC_EDIT_GKM, m_realKM);
	DDX_Text(pDX, IDC_EDIT_GKv, m_realKv);
	DDX_Control(pDX, IDC_CHECK_Na2, m_checkNa2);
	DDX_Control(pDX, IDC_CHECK_K2, m_checkK2);
	DDX_Control(pDX, IDC_CHECK_Ca2, m_checkCa2);
	DDX_Control(pDX, IDC_CHECK_KM2, m_checkKM2);
	DDX_Control(pDX, IDC_CHECK_Kv2, m_checkKv2);
	DDX_Control(pDX, IDC_CHARTCTRL, m_ChartCtrl1);
	//DDX_Control(pDX, IDC_CHARTCTRL, m_ChartCtrl2);
	DDX_Control(pDX, IDC_CUSTOM_ChartCtrl2, m_ChartCtrl2);
	DDX_Control(pDX, IDC_CUSTOM_ChartCtrl3, m_ChartCtrl3);
	DDX_Control(pDX, IDC_CUSTOM_ChartCtrl4, m_ChartCtrl4);
	DDX_Control(pDX, IDC_CUSTOM_ChartCtrl5, m_ChartCtrl5);
	DDX_Control(pDX, IDC_PROGRESS_Step, m_progress);
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
//ON_EN_CHANGE(IDC_MFCEDITBROWSE_Filepath, &CGA_SpikeDlg::OnEnChangeMfceditbrowseFilepath)
ON_BN_CLICKED(IDC_BUTTON_ADD, &CGA_SpikeDlg::OnBnClickedButtonAdd)
ON_BN_CLICKED(IDC_BUTTON_DEL, &CGA_SpikeDlg::OnBnClickedButtonDel)
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
	CString   str_MaxGeneration = "5";
	CString   str_V = "-69";
	CString   str_gL = "0.3";
	CString   str_Capacitance = "1";
	CString   str_Population_size = "2048";
	CString   str_Crossover = "0.5";
	CString   str_Mutations = "0.001";
	CString   str_Time = "100";
//	m_splitbtn.SetDropDownMenu(IDR_MENU1, 0);
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
	m_wndList.SetHeadings(_T("start, 120; duration, 120; value, 120"));
	m_wndList.SetGridLines(TRUE);
	//m_List.InsertColumn(0, "start", LVCFMT_LEFT, 100);
	//m_List.InsertColumn(1, "duration", LVCFMT_LEFT, 100);
	//m_List.InsertColumn(2, "value", LVCFMT_LEFT, 100);
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	m_progress.SetRange(0, 1000);
	
	CButton* pBtn = (CButton*)GetDlgItem(IDC_CHECK_GPU);
	int count;
	// get the cuda device count
	cudaGetDeviceCount(&count);
	if (count)
		pBtn->SetCheck(1);
	///////////////////////////////
	//////////////////////////////
	CButton* pBtn1 = (CButton*)GetDlgItem(IDC_CHECK_Na);
	pBtn1->SetCheck(1);
	CButton* pBtn2 = (CButton*)GetDlgItem(IDC_CHECK_K);
	pBtn2->SetCheck(1);
	CButton* pBtn3 = (CButton*)GetDlgItem(IDC_CHECK_Ca);
	pBtn3->SetCheck(1);
	CButton* pBtn4 = (CButton*)GetDlgItem(IDC_CHECK_KM);
	pBtn4->SetCheck(1);
	CButton* pBtn5 = (CButton*)GetDlgItem(IDC_CHECK_Kv);
	pBtn5->SetCheck(1);

	CButton* pBtn21 = (CButton*)GetDlgItem(IDC_CHECK_Na2);
	pBtn21->SetCheck(1);
	CButton* pBtn22 = (CButton*)GetDlgItem(IDC_CHECK_K2);
	pBtn22->SetCheck(1);
	CButton* pBtn23 = (CButton*)GetDlgItem(IDC_CHECK_Ca2);
	pBtn23->SetCheck(1);
	CButton* pBtn24 = (CButton*)GetDlgItem(IDC_CHECK_KM2);
	pBtn24->SetCheck(1);
	CButton* pBtn25 = (CButton*)GetDlgItem(IDC_CHECK_Kv2);
	pBtn25->SetCheck(1);
	///////////////////////////////////
	//SetCheck(1)表示设置复选框为“选中”状态；
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
UINT postThread2(LPVOID pParam)
{
	CGA_SpikeDlg*p2 = (CGA_SpikeDlg*)pParam;
	p2->RunTimer2();
	return 0;
}
void CGA_SpikeDlg::threadNew()
{
	AfxBeginThread(postThread, this);
}
void CGA_SpikeDlg::threadNew2()
{
	AfxBeginThread(postThread2, this);
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
void CGA_SpikeDlg::RunTimer2()
{
	SetTimer(2, 100, NULL);
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
vector<float> MHH_SpikeTime(float *data, double Mtime, double TimeStep, int &length)
{
	int size = int(Mtime / TimeStep);
	int k = 0;
	vector<float>Mdata_copy;
	//float *Mdata_copy = new float[spike_length];
	for (int i = 1; i < size - 1; i++)
	{
		if (data[i]>data[i - 1] && data[i] > data[i + 1] && data[i] >= 0)
		{
			Mdata_copy.push_back(i*TimeStep);
			k++;
		}
	}

	if (k == 0)
	{
		k = 1;
		Mdata_copy.push_back(Mtime / TimeStep);
	}
	length = k;
	//for (int j = 0; j < spike_length; j++)
	//{
	//	Mdata_copy[j] = j*10000;
	//}
	//length=spike_length;
	return Mdata_copy;
}

void CGA_SpikeDlg::MyCurvesHH(vector<float> m_I, vector<float> VrealSpike, double TimeStep)
{
	m_ChartCtrl1.EnableRefresh(false);
	float I_max = 0; float V_max = 0; float V_min = 100;
	for (int i = 0; i < m_I.size(); i++)
	{
		if (m_I[i]>I_max) I_max = m_I[i];
		if (VrealSpike.size() > 0)
		{
			if (VrealSpike[i] > V_max) V_max = VrealSpike[i];
			if (VrealSpike[i] < V_min) V_min = VrealSpike[i];
		}
	}
	CChartStandardAxis* pBottomAxis =
		m_ChartCtrl1.CreateStandardAxis(CChartCtrl::BottomAxis);
	pBottomAxis->SetMinMax(0, m_I.size()*TimeStep);
	CChartStandardAxis* pLeftAxis =
		m_ChartCtrl1.CreateStandardAxis(CChartCtrl::LeftAxis);
	pLeftAxis->SetMinMax(0, I_max);
	CChartStandardAxis* pRightAxis =
		m_ChartCtrl1.CreateStandardAxis(CChartCtrl::RightAxis);
	pRightAxis->SetMinMax(V_min, V_max);
	//pBottomAxis->SetTickIncrement(false, 1.0);
	//pBottomAxis->SetDiscrete(false);
	CChartLineSerie* pSeries = m_ChartCtrl1.CreateLineSerie();
	CChartLineSerie* pSeries2 = m_ChartCtrl1.CreateLineSerie(false, true);
	double *XVal = new double[m_I.size()];
	double *YVal = new double[m_I.size()];
	double *YVal2 = new double[m_I.size()];
	for (int i = 0; i<m_I.size(); i++)
	{
		XVal[i] =  i*TimeStep ;
		YVal[i] = m_I[i];
		if (VrealSpike.size()>0)
			YVal2[i] = VrealSpike[i];
		//cout << XVal[i]<<"//"<< YVal[i] << " ";
	}
	pSeries->SetPoints(XVal, YVal, m_I.size());
	pSeries2->SetPoints(XVal, YVal2, VrealSpike.size());
	m_ChartCtrl1.EnableRefresh(true);
	delete XVal;
	delete YVal;
	delete YVal2;
}
///////////////////////////
///////////////////////////////
void CGA_SpikeDlg::OnBnClickedButtonRun()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(true);
	if (myexit)
	{
		KillTimer(1);
	}
	int count;
	cudaGetDeviceCount(&count);
	if (count==0)
	{
		MessageBox("There is no device");
	}
	/////////////////////////
	threadNew();
	threadNew2();
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
	GetDlgItemText(IDC_MFCEDITBROWSE_Filepath, FilePathName);
	GetDlgItemText(IDC_MFCEDITBROWSE_Ifile, FilePathName_I);
	freopen("output.txt", "w", stdout);

	srand(1900);
	srand(static_cast<unsigned>(time(0)));
	//float *spike_data, *current_data, *spike_TestData;
	//float gNa = 120;
	//float gK = 36;
	double Mtime = m_time;//ms
	double tempVB = m_V;//初值
	double TimeStep = 0.01;
	//int I_length = int(Mtime / TimeStep);
	//double m_I = 10;//电流
	vector<float> m_I;
	vector<float>Mtemp_data;
	int FlagParameter[5] = { m_checkNa2.GetCheck(), m_checkK2.GetCheck(), m_checkCa2.GetCheck(), m_checkKM2.GetCheck(), m_checkKv2.GetCheck() };//Na,K,KM,Kv,Ca
	int FlagParameterReal[5] = { m_checkNa.GetCheck(), m_checkK.GetCheck(), m_checkCa.GetCheck(), m_checkKM.GetCheck(), m_checkKv.GetCheck() };//Na,K,KM,Kv,Ca
	vector<float>VrealSpike;
	//Parameter_Bound.length = 2;
	std::string strFilePathName;
	strFilePathName = FilePathName.GetBuffer(0);
	std::string strFilePathI;
	strFilePathI = FilePathName_I.GetBuffer(0);
	if (strFilePathI == "")
	{
		int nitem = m_wndList.GetItemCount();
		vector<vector<float>>Vn_I;
		for (int i = 0; i < nitem; i++)
		{
			vector<float> Vrow;
			for (int j = 0; j < 3; j++)
			{
				CString tempI = m_wndList.GetItemText(i, j);
				float tpI=(float)atof((char *)(LPTSTR)(LPCTSTR)tempI);
				Vrow.push_back(tpI);
			}
			Vn_I.push_back(Vrow);
		}
		int kcount = 0;
		for (int i = 0; i < Vn_I.size(); i++)
		{
			while (kcount*TimeStep < Vn_I[i][0]) { m_I.push_back(0); kcount++; }
			while (kcount*TimeStep < Vn_I[i][0] + Vn_I[i][1]) { m_I.push_back(Vn_I[i][2]); kcount++; }
		}
	}
	else
	{
		m_I = Read_Txt(strFilePathI);
	}

	if (strFilePathName != "")
	{
		Mtemp_data = Read_Txt(strFilePathName);
		Parameter_.spike_data_num = Mtemp_data.size();
	}
	else
	{
		Mtime = m_I.size()*TimeStep;
		float *pArrayA = new float[m_I.size()];
		HMODULE hDll = LoadLibrary("CHH_dll.dll");
		if (hDll != NULL)
		{
			AddFunc HH_Entrance = (AddFunc)GetProcAddress(hDll, "HH_Entrance");
			if (HH_Entrance != NULL)
			{
				HH_Entrance(m_realNa, m_realK, m_realCa, m_realKM, m_realKv, Mtime, tempVB, TimeStep, m_I, FlagParameterReal, m_gL, m_C, pArrayA);
			}
			FreeLibrary(hDll);
		}
		else
		{
			{
				printf("not find dll");
			}
		}
		
		for (int i = 0; i < int(Mtime / TimeStep); i++)
		{
			VrealSpike.push_back(pArrayA[i]);
		}
		//MyCurves(VrealSpike, TimeStep);
		Mtemp_data = MHH_SpikeTime(pArrayA, Mtime, TimeStep, Parameter_.spike_data_num);
		delete pArrayA;
	}
	Parameter_.spike_data = new float[Mtemp_data.size()];
	convert_data(Mtemp_data, Parameter_.spike_data);
	//float ans = solveGPU_cpp(Parameter_, Mtime, tempVB, TimeStep, m_I, FlagParameter, Parameter_Bound, MaxGeneration, gL, C, POPULATION_SIZE, crossver, mutations, strResult);
	//////////////////////////////////////////
	MyCurvesHH(m_I, VrealSpike, TimeStep);
	//////////////////////////////////////////
	//////pSthData pDataValue1 = new SthData();
	//////memset(pDataValue1, 0x00, sizeof(SthData));
	Info.Parameter_ = Parameter_;
	Info.Mtime = Mtime;
	Info.tempVB = tempVB;
	Info.m_I = m_I;
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
	Info.TimeStep = TimeStep;
	//////threadNew2(Info);

	/////////////////////////////////
	update = TRUE;
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

//void CGA_SpikeDlg::DrawPictureB(vector<vector<float>> mdata_g, vector<float>mdata_score)
//{   ///////////////gNa,gK
//	//Invalidate();
//	CPaintDC dc(this);
//	RECT rect;
//	CDC memDC;
//	CDC *pDc = GetDlgItem(IDC_STATIC_PictureB)->GetDC();
//	GetDlgItem(IDC_STATIC_PictureB)->GetClientRect(&rect);
//	InvalidateRect(&rect);
//	memDC.CreateCompatibleDC(pDc);//创建与目标DC相兼容的内存DC，
//	CBitmap memBitmapB;
//	int w = rect.right - rect.left;
//	int h = rect.bottom - rect.top;
//	//创建一个内存中的图像
//	///CBitmap memBitmapB;
//	memBitmapB.CreateCompatibleBitmap(pDc, w, h);//根据目标DC创建位图
//	memDC.SelectObject(&memBitmapB);//把位图选入内存DC
//	memDC.FillSolidRect(rect.left,rect.top,rect.right,rect.bottom, pDc->GetBkColor());//按原来背景填充客户区，不然会是黑色
//	memDC.MoveTo(30, h - 30);
//	memDC.LineTo(w - 30, h - 30);
//	memDC.MoveTo(30, h - 30);
//	memDC.LineTo(30, 30);
//	memDC.MoveTo(w - 35, h - 25);
//	memDC.LineTo(w - 30, h - 30);
//	memDC.LineTo(w - 35, h - 35);
//	memDC.MoveTo(25, 35);
//	memDC.LineTo(30, 30);
//	memDC.LineTo(35, 35);
//	//////x轴刻度
//	//////////y轴刻度
//	double dividx = (m_Na_up - m_Na_low) / 10;
//	double dividy = (m_K_up - m_K_low) / 10;
//
//	for (int j = 1; j <= 10; j++)
//	{
//		CString strY;
//		CString strX;
//		strY.Format("%g", m_K_low + j * dividy);
//		strX.Format("%g", m_Na_low + j * dividx);
//		memDC.TextOut(23 + j*1.0 * (w - 60) / 10, h - 24, strX);
//		memDC.TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
//	}
//	//pDc->TextOut(5, h - 24, 0);
//	//////////////////
//	CBrush *oldbrush;
//	CBrush brush;
//	brush.CreateSolidBrush(RGB(255, 0, 0));
//	//brush.CreateSolidBrush(RGB(255, 255, 0));
//	oldbrush = memDC.SelectObject(&brush);
//	for (int j = 0; j < mdata_g.size(); j++)////描点
//	{
//		float nscore = (mdata_score[j] * 200 / 20>255) ? 255 : (mdata_score[j] * 200 / 20);
//		brush.DeleteObject();
//		brush.CreateSolidBrush(RGB(nscore, nscore ,nscore));
//		oldbrush = memDC.SelectObject(&brush);
//		float m_x = mdata_g[j][0];
//		float m_y = mdata_g[j][1];
//		double x = 30 + (m_x - m_Na_low)*(w - 60) / (m_Na_up - m_Na_low);
//		double y = h - 30 - (m_y - m_K_low)*(h - 60) / (m_K_up - m_K_low);
//		memDC.Ellipse(x - 2, y - 2, x + 2, y + 2);
//	}
//	memDC.SelectObject(oldbrush);
//	pDc->BitBlt(0, 0, w, h, &memDC, 0, 0, SRCCOPY);//将内存DC上的图象拷贝到前台
//	memDC.DeleteDC();                                       //删除DC
//	///memBitmapB.DeleteObject();                                        //删除位图
//	ReleaseDC(pDc);
//}
void CGA_SpikeDlg::DrawPictureB(vector<vector<float>> mdata_g, vector<vector<float>>mdata_score, threadInfo Info)
{
	if (mdata_g.size() > 0)
	{
		m_ChartCtrl2.EnableRefresh(false);
		m_ChartCtrl2.RemoveAllSeries();//先清空 
		CChartStandardAxis* pBottomAxis =
			m_ChartCtrl2.CreateStandardAxis(CChartCtrl::BottomAxis);
		pBottomAxis->SetMinMax(Info.Parameter_Bound[0].g[0], Info.Parameter_Bound[0].g[1]);
		CChartStandardAxis* pLeftAxis =
			m_ChartCtrl2.CreateStandardAxis(CChartCtrl::LeftAxis);
		pLeftAxis->SetMinMax(Info.Parameter_Bound[1].g[0], Info.Parameter_Bound[1].g[1]);
		//pBottomAxis->SetDiscrete(false);
		//m_ChartCtrl2.RemoveAllSeries();//先清空 
		CChartPointsSerie* pSeries = m_ChartCtrl2.CreatePointsSerie();
		double *XVal = new double[mdata_g.size()];
		double *YVal = new double[mdata_g.size()];
		for (int i = 0; i < mdata_g.size(); i++)
		{
			XVal[i] = mdata_g[i][0];
			YVal[i] = mdata_g[i][1];

		}
		pSeries->SetPoints(XVal, YVal, mdata_g.size());
		m_ChartCtrl2.EnableRefresh(true);
		delete XVal;
		delete YVal;
	}

}
////////////////////////////////////
//void CGA_SpikeDlg::DrawPictureA(vector<float>mdata_score)
//{
//	CPaintDC dc(this);
//	RECT rect;
//	///////////////score
//	CDC *pDcB = GetDlgItem(IDC_STATIC_PictureA)->GetDC();
//	GetDlgItem(IDC_STATIC_PictureA)->GetClientRect(&rect);
//	int w = rect.right - rect.left;
//	int h = rect.bottom - rect.top;
//	pDcB->MoveTo(30, h - 30);
//	pDcB->LineTo(w - 30, h - 30);
//	pDcB->MoveTo(30, h - 30);
//	pDcB->LineTo(30, 30);
//	pDcB->MoveTo(w - 35, h - 25);
//	pDcB->LineTo(w - 30, h - 30);
//	pDcB->LineTo(w - 35, h - 35);
//	pDcB->MoveTo(25, 35);
//	pDcB->LineTo(30, 30);
//	pDcB->LineTo(35, 35);
//	//////
//	//////////y轴刻度
//	double dividx = (m_MaxGeneration - 0) / 10;
//	double dividy = (50 - 0) / 10;
//
//	for (int j = 1; j <= 10; j++)
//	{
//		CString strY;
//		strY.Format("%g", 0 + j * dividy);
//		pDcB->TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
//	}
//	////
//	//x轴刻度
//	for (int j = 1; j <= m_MaxGeneration; j++)
//	{
//		CString strX;
//		strX.Format("%d", j);
//		pDcB->TextOut(23 + j*1.0 * (w - 60) / m_MaxGeneration, h - 24, strX);
//	}
//	//pDc->TextOut(5, h - 24, 0);
//	//////////////////
//	for (int j = 0; j <mdata_score.size(); j++)////描点
//	{
//		float m_x = m_gen;
//		float m_y = mdata_score[j];
//		double x = 30 + (m_x - 0)*(w - 60) / (m_MaxGeneration - 0);
//		double y = h - 30 - (m_y - 0)*(h - 60) / (50 - 0);
//		pDcB->Ellipse(x - 2, y - 2, x + 2, y + 2);
//	}
//	//////////////////
//	CBrush *oldbrush;
//	CBrush brush;
//	brush.CreateSolidBrush(RGB(255, 0, 0));
//	oldbrush = pDcB->SelectObject(&brush);
//	ReleaseDC(pDcB);
//}
void CGA_SpikeDlg::DrawPictureA(vector<vector<float>>mdata_score)
{
	if (mdata_score.size() > 0)
	{
		m_ChartCtrl3.EnableRefresh(false);
		m_ChartCtrl3.RemoveAllSeries();//先清空 
		
		CChartStandardAxis* pBottomAxis =
			m_ChartCtrl3.CreateStandardAxis(CChartCtrl::BottomAxis);
		pBottomAxis->SetMinMax(0, m_MaxGeneration);
		pBottomAxis->SetTickIncrement(false, 1);
		CChartStandardAxis* pLeftAxis =
			m_ChartCtrl3.CreateStandardAxis(CChartCtrl::LeftAxis);
		pLeftAxis->SetMinMax(0,1);
		///////////////////////////
		//CChartStandardAxis* pBottomAxis6 =
		//	m_ChartCtrl6.CreateStandardAxis(CChartCtrl::BottomAxis);
		//pBottomAxis->SetMinMax(0, m_MaxGeneration);
		//CChartStandardAxis* pLeftAxis6 =
		//	m_ChartCtrl6.CreateStandardAxis(CChartCtrl::LeftAxis);
		//pLeftAxis6->SetMinMax(0, 2);
		//////////////////////////////////
		//CChartStandardAxis* pRightAxis =
		//	m_ChartCtrl1.CreateStandardAxis(CChartCtrl::RightAxis);
		//pRightAxis->SetMinMax(0, 2);
		//pBottomAxis->SetDiscrete(false);
		//m_ChartCtrl2.RemoveAllSeries();//先清空 
		CChartPointsSerie* pSeries = m_ChartCtrl3.CreatePointsSerie();
		CChartPointsSerie* pSeries6 = m_ChartCtrl3.CreatePointsSerie();
		//CChartPointsSerie* pSeries6 = m_ChartCtrl6.CreatePointsSerie();
		double *XVal = new double[mdata_score[0].size()];
		double *YVal = new double[mdata_score[0].size()];
		for (int j = 0; j < mdata_score.size(); j++)
		{
			double Ymin = 100;
			for (int i = 0; i < mdata_score[0].size(); i++)
			{
				XVal[i] = j+1;
				YVal[i] = 1 - exp(-0.1*mdata_score[j][i]);
				if (YVal[i] < Ymin) Ymin = YVal[i];
			}
			pSeries->SetColor(RGB(210, 210, 210));
			pSeries->AddPoints(XVal, YVal, mdata_score[0].size());
			pSeries6->SetColor(RGB(250, 20, 20));
			pSeries6->AddPoint(j+1, Ymin);
		}

		m_ChartCtrl3.EnableRefresh(true);
		//m_ChartCtrl6.EnableRefresh(true);
		delete XVal;
		delete YVal;
		//delete YVal6;
	}
}
/////////////////////////
//void CGA_SpikeDlg::DrawPictureC(double duration)
//{
//	CPaintDC dc(this);
//	RECT rect;
//	///////////////score
//	CDC *pDcB = GetDlgItem(IDC_STATIC_PictureC)->GetDC();
//	GetDlgItem(IDC_STATIC_PictureC)->GetClientRect(&rect);
//	int w = rect.right - rect.left;
//	int h = rect.bottom - rect.top;
//	CBitmap memBitmap;
//	pDcB->MoveTo(30, h - 30);
//	pDcB->LineTo(w - 30, h - 30);
//	pDcB->MoveTo(30, h - 30);
//	pDcB->LineTo(30, 30);
//	pDcB->MoveTo(w - 35, h - 25);
//	pDcB->LineTo(w - 30, h - 30);
//	pDcB->LineTo(w - 35, h - 35);
//	pDcB->MoveTo(25, 35);
//	pDcB->LineTo(30, 30);
//	pDcB->LineTo(35, 35);
//	//////
//	//////////y轴刻度
//	double dividx = (m_MaxGeneration - 0) / 10;
//	double dividy = (500 - 0) / 10;
//
//	for (int j = 1; j <= 10; j++)
//	{
//		CString strY;
//		strY.Format("%g", 0 + j * dividy);
//		pDcB->TextOut(5, h - 30 - j*1.0 * (h - 60) / 10, strY);
//	}
//	////
//	//x轴刻度
//	for (int j = 1; j <= m_MaxGeneration; j++)
//	{
//		CString strX;
//		strX.Format("%d", j);
//		pDcB->TextOut(23 + j*1.0 * (w - 60) / m_MaxGeneration, h - 24, strX);
//	}
//	//pDc->TextOut(5, h - 24, 0);
//	//////////////////
//	float m_x = m_gen;
//	float m_y = duration;
//	double x = 30 + (m_x - 0)*(w - 60) / (m_MaxGeneration - 0);
//	double y = h - 30 - (m_y - 0)*(h - 60) / (500 - 0);
//	pDcB->Ellipse(x - 2, y - 2, x + 2, y + 2);
//	//////////////////
//	CBrush *oldbrush;
//	CBrush brush;
//	brush.CreateSolidBrush(RGB(255, 0, 0));
//	ReleaseDC(pDcB);
//}
void CGA_SpikeDlg::DrawPictureC(double duration)
{
}
//////////////////
void CGA_SpikeDlg::OnTimer(UINT_PTR nIDEvent)
{
	// TODO:  在此添加消息处理程序代码和/或调用默认值
	switch (nIDEvent)
	{
		case 1:   //定时器1处理函数，定时发送数据进行更新
		{	

						  //cout << m_gen << endl;
			mtx.lock();
			//cout << m_gen << endl;
			if (mdata_g.size() > 0 && mdata_score.size()>0)
			{
				DrawPictureB(mdata_g, mdata_score, Info);
				DrawPictureA(mdata_score);
				DrawPictureC(durationMy);
				mdata_g.clear();
			}
			mtx.unlock();
			break;
		}
		case 2:   //定时器2处理函数
		{
					  
			//if (durationMy)
			//{
				clock_t Mend= clock();
				if (m_gen==1)
					pos += 5;
				//int pos = 1000 * durationMy / (CLOCKS_PER_SEC*durationMy*(m_MaxGeneration));
				else if (m_gen>1)
					pos += 1000 / (10*CLOCKS_PER_SEC*durationMy*(m_MaxGeneration-m_gen));
				
				m_progress.SetPos(pos);
				break;
			//}

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


}





void CGA_SpikeDlg::OnBnClickedCheck1()
{
	// TODO:  在此添加控件通知处理程序代码
}


void CGA_SpikeDlg::OnBnClickedButtonAdd()
{
	UpdateData(TRUE);

	
	//AR: Add the CSomeClass object to the List
	//int nPos = m_List.GetItemCount();
	////m_List.InsertItem(nPos, m_sStrValue);

	CString sIntValue1,sIntValue2, sIntValue3;
	sIntValue1.Format("%g", m_current_start);
	sIntValue2.Format("%g", m_current_duration);
	sIntValue3.Format("%g", m_current_val);
	//m_List.InsertItem(nPos, sIntValue1);
	//m_List.InsertItem(nPos, sIntValue2);
	//m_List.InsertItem(nPos, sIntValue3);
	////m_List.SetItemText(nPos, 0, sIntValue1);
	////m_List.SetItemText(nPos, 1, sIntValue2);
	////m_List.SetItemText(nPos, 2, sIntValue3);
	m_wndList.InsertItem(INT_MAX, sIntValue1, sIntValue2, sIntValue3);

	
}


void CGA_SpikeDlg::OnBnClickedButtonDel()
{
	m_wndList.DeleteItem(m_wndList.GetFirstSelectedItem(), TRUE, ItemdataProc, (LPARAM)this);// TODO:  在此添加控件通知处理程序代码
}
BOOL CGA_SpikeDlg::ItemdataProc(DWORD dwData, LPARAM lParam)
{
	// TODO: Process your item data here

	// Please return TRUE to proceed the deletion, return FALSE to abort.
	return TRUE;
}


