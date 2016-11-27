
#pragma once
using namespace std;

#include <stdio.h>
#include <tchar.h>
#include <time.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <Windows.h>
#include <cstring>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <math.h>


#define DIMENSION 3
#define A_X 0
#define A_Y 1
#define A_Z 2
#define G 9.8
#define ON 1
#define OFF 0

void gauss(double *matrix,double *B,int N);
int newton();
int main_hy0_dpdq();
int main_hy0();
int main_n1();
int q_newton();
int SQP_method();
int SUMT_IPM();
int MM_method();
int MM_method2();
int QP_with_2_newton_way();
int QP_with_2_newton_way_SDM();
int h_p();
int h_q();

void PM_AVS(int t_max,int t, double *qn);
void MM_AVS(int t_max, int t,double Dt, double mi, double *qn,double *pn);

int main()	
{
	system("mkdir E_QN");
	system("mkdir p_QN");
	system("mkdir q_QN");
	system("mkdir lam_QN");
	system("mkdir mu_QN");
	system("mkdir Tr_QN");
	system("mkdir h_QN");
	system("mkdir dhdt_QN");

	system("mkdir E_QP");
	system("mkdir count_QP");
	system("mkdir h_QP");
	system("mkdir dhdt_QP");
	system("mkdir p_QP");
	system("mkdir q_QP");
	system("mkdir lam_QP");
	system("mkdir mu_QP");
	system("mkdir Tr_QP");


	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=0.1;
	double nG[DIMENSION]={0,0,1};
	double n[DIMENSION]={0,0,1};
	double a[DIMENSION]={0,0,0};

	double lam=0;
	double mu=0;

	double old_lam=0;
	double old_mu=0;

	double Ln=0;
	double L=0;
	double dL[2]={0,0};
	
	double h=0;		//h(q)での定義が正解
	double dh[2]={0,0};
	double dhdt=0;
	double d_dhdt[2]={0,0};

	double Tr=0;
	double dTr[2]={0,0};


	double r=0.1;
	double ep=1e-5;
	double ep_min=1e-5;
	double theta_h=0;
	double theta_dh=0;
	int c_max=5000;
	double E=1;
	double E_min=1;

	double *d=new double [Nx];	
	double *Nr=new double [Nx];
	double *B=new double [Nx*Nx];

	int t=0;
	int t_max=1000;

	double p0[DIMENSION]={0,0,0};
	double q0[DIMENSION]={0,0,5};
	double pn[DIMENSION]={0,0,0};
	double qn[DIMENSION]={0,0,0};
	double pn1[DIMENSION]={0,0,0};
	double qn1[DIMENSION]={0,0,0};
	double pn1_2[DIMENSION]={0,0,0};


	ofstream fs("intial_data.csv");
	fs<<"r, "<<r<<endl;
	fs<<"Dt, "<<Dt<<endl;
	fs<<"mi, "<<mi<<endl;
	fs<<"p0, "<<p0[A_X]<<","<<p0[A_Y]<<","<<p0[A_Z]<<endl;
	fs<<"q0, "<<q0[A_X]<<","<<q0[A_Y]<<","<<q0[A_Z]<<endl;
	fs<<"t_max, "<<t_max<<endl;
	fs<<"ep_min, "<<ep_min<<endl;
	fs<<"ep, "<<ep<<endl;
	fs.close();


	cout<<"Calculation Start"<<endl;


	pn[A_X]=p0[A_X];
	pn[A_Y]=p0[A_Y];
	pn[A_Z]=p0[A_Z];
	qn[A_X]=q0[A_X];
	qn[A_Y]=q0[A_Y];
	qn[A_Z]=q0[A_Z];

	double E0=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);


	while(t<t_max)
	{
		t++;
		clock_t t0=clock();

		double pn_1[DIMENSION]={pn[A_X], pn[A_Y], pn[A_Z]};
		double qn_1[DIMENSION]={qn[A_X], qn[A_Y], qn[A_Z]};

		pn1_2[A_X]=pn[A_X]-Dt*mi*G*nG[A_X];	pn1_2[A_Y]=pn[A_Y]-Dt*mi*G*nG[A_Y];	pn1_2[A_Z]=pn[A_Z]-Dt*mi*G*nG[A_Z];
		qn[A_X]+=Dt/mi*pn1_2[A_X];	qn[A_Y]+=Dt/mi*pn1_2[A_Y];	qn[A_Z]+=Dt/mi*pn1_2[A_Z];
		pn[A_X]+=-Dt*mi*G*nG[A_X];	pn[A_Y]+=-Dt*mi*G*nG[A_Y];	pn[A_Z]+=-Dt*mi*G*nG[A_Z];

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream fdh("dhdt.csv", ios::trunc);
			fdh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();
			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*( (qn[A_X]-a[A_X])*n[A_X] + (qn[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]-a[A_Z])*n[A_Z]);
		dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);
		if(h>0 || dhdt>0)
		{
			pn[A_X]=pn_1[A_X];	pn[A_Y]=pn_1[A_Y];	pn[A_Z]=pn_1[A_Z];
			qn[A_X]=qn_1[A_X];	qn[A_Y]=qn_1[A_Y];	qn[A_Z]=qn_1[A_Z];//

			cout<<"接触"<<t<<" ,QP_method"<<endl;
			
			theta_dh=0;
			theta_h=0;
			r=0.1;
			lam=1;
			mu=1;

			Ln=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])-mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);

			E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_QP/E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count_QP/count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./p_QP/p"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./q_QP/q"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
	
					stringstream se_lam;
					se_lam<<"./lam_QP/lam"<<t<<".csv";
					ofstream init_lam(se_lam.str(), ios::trunc);
					init_lam.close();

					stringstream se_mu;
					se_mu<<"./mu_QP/mu"<<t<<".csv";
					ofstream init_mu(se_mu.str(), ios::trunc);
					init_mu.close();

					stringstream se_Tr;
					se_Tr<<"./Tr_QP/Tr"<<t<<".csv";
					ofstream init_Tr(se_Tr.str(), ios::trunc);
					init_Tr.close();

					stringstream se_h;
					se_h<<"./h_QP/h"<<t<<".csv";
					ofstream init_h(se_h.str(), ios::trunc);
					init_h.close();

					stringstream se_dhdt;
					se_dhdt<<"./dhdt_QP/dhdt"<<t<<".csv";
					ofstream init_dhdt(se_dhdt.str(), ios::trunc);
					init_dhdt.close();
				}

				if(count_min>c_max && h<=0)	break;

				old_lam=lam;
				old_mu=mu;

				for(int i=0;i<Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<Nx;j++)
					{
						if(j==i)	B[i*Nx+j]=1;
						else
						{
							B[i*Nx+j]=0;
						}
					}
				}

				pn1[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam)+Dt*n[A_X]*mu;	
				pn1[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam)+Dt*n[A_Y]*mu;
				pn1[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam)+Dt*n[A_Z]*mu;

				qn1[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam));
				qn1[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam));
				qn1[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam));


				L=0.5/mi*(pn1[A_X]*pn1[A_X] + pn1[A_Y]*pn1[A_Y] + pn1[A_Z]*pn1[A_Z])-mi*G*(qn1[A_X]*nG[A_X] + qn1[A_Y]*nG[A_Y] + qn1[A_Z]*nG[A_Z]);

				Tr=(L-Ln)*(L-Ln);

				h=-1*( (qn1[A_X] - a[A_X])*n[A_X] + (qn1[A_Y] - a[A_Y])*n[A_Y] + (qn1[A_Z] - a[A_Z])*n[A_Z]);
				dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

				if(h+theta_h>0)	Tr+=0.5*r*(h+theta_h)*(h+theta_h);
				if(dhdt+theta_dh>0)	Tr+=0.5*r*(dhdt+theta_dh)*(dhdt+theta_dh);

				dL[0]=1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z])-G*Dt*Dt*(n[A_X]*nG[A_X]+n[A_Y]*nG[A_Y]+n[A_Z]*nG[A_Z]);
				dL[1]=1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

				dTr[0]=2*dL[0]*(L-Ln);
				dTr[1]=2*dL[1]*(L-Ln);

				dh[0]=-Dt*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
				dh[1]=0;
				d_dhdt[0]=-Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
				d_dhdt[1]=-Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);

				if(h+theta_h>0)	dTr[0]+=r*dh[0]*(h+theta_h);	
				if(dhdt+theta_dh>0)
				{
					dTr[0]+=r*d_dhdt[0]*(dhdt+theta_dh);
					dTr[1]+=r*d_dhdt[1]*(dhdt+theta_dh);
				}

				E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);

				int count=0;
 				if(E<ep)	break;
					while(E>ep)
					{
						count++;
						if(count==1)
						{
							stringstream se;
							se<<"./E_QN/E_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream init_e(se.str(), ios::trunc);
							init_e.close();

							stringstream ssp;
							ssp<<"./p_QN/p_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream init_p(ssp.str(), ios::trunc);
							init_p.close();
	
							stringstream ssq;
							ssq<<"./q_QN/q_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream init_q(ssq.str(), ios::trunc);
							init_q.close();

							stringstream sslam;
							sslam<<"./lam_QN/lam_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream init_lam(sslam.str(), ios::trunc);
							init_lam.close();
	
							stringstream ssmu;
							ssmu<<"./mu_QN/mu_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream init_mu(ssmu.str(), ios::trunc);
							init_mu.close();

							stringstream se_Tr;
							se_Tr<<"./Tr_QN/Tr_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fTr(se_Tr.str(), ios::trunc);
							fTr.close();

							stringstream se_h;
							se_h<<"./h_QN/h_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fh(se_h.str(), ios::trunc);
							fh.close();

							stringstream se_dhdt;
							se_dhdt<<"./dhdt_QN/dhdt_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fdhdt(se_dhdt.str(), ios::trunc);
							fdhdt.close();

						}
						if(count>c_max)	break;

						double lam_k=lam;
						double mu_k=mu;
						double dTr_k[2]={dTr[0],dTr[1]};

						Nr[0]=dTr[0];	Nr[1]=dTr[1];
						if(count>1)
						{
	/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
							cout<<"B=";
							for(int i=0;i<3*Nx;i++)
							{
								for(int j=0;j<3*Nx;j++)
								{
									cout<<B[i*3*Nx+j]<<", ";
								}
								cout<<endl;
							}*/

							gauss(B,Nr,Nx);
						}
						d[0]=-1*Nr[0];	d[1]=-1*Nr[1];

						//cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
	//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

						double Tr_min=Tr;
						double a_min=0.1;//1e-4;
						double pn1_a[DIMENSION]={0,0,0};
						double qn1_a[DIMENSION]={0,0,0};
						double L_a=0,	h_a=0,	dhdt_a=0,	Tr_a=0;


						for(int i=0;i<10000;i++)
						{
							double alpha=(i+1)*1e-4;
							double lam_a=lam+d[0]*alpha;
							double mu_a=mu+d[1]*alpha;
						
							pn1_a[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam_a)+Dt*n[A_X]*mu_a;	
							pn1_a[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam_a)+Dt*n[A_Y]*mu_a;
							pn1_a[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam_a)+Dt*n[A_Z]*mu_a;

							qn1_a[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam_a));
							qn1_a[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam_a));
							qn1_a[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam_a));

							L_a=0.5/mi*(pn1_a[A_X]*pn1_a[A_X] + pn1_a[A_Y]*pn1_a[A_Y] + pn1_a[A_Z]*pn1_a[A_Z])-mi*G*(qn1_a[A_X]*nG[A_X] + qn1_a[A_Y]*nG[A_Y] + qn1_a[A_Z]*nG[A_Z]);
							Tr_a=(L_a-Ln)*(L_a-Ln);

							h_a=-1*( (qn1_a[A_X] - a[A_X])*n[A_X] + (qn1_a[A_Y] - a[A_Y])*n[A_Y] + (qn1_a[A_Z] - a[A_Z])*n[A_Z]);
							dhdt_a=-1/mi*(n[A_X]*pn1_a[A_X]+n[A_Y]*pn1_a[A_Y]+n[A_Z]*pn1_a[A_Z]);

							if(h_a+theta_h>0)	Tr_a+=0.5*r*(h_a+theta_h)*(h_a+theta_h);

							if(dhdt_a+theta_dh>0)	Tr_a+=0.5*r*(dhdt_a+theta_dh)*(dhdt_a+theta_dh);

							if(Tr_a<Tr_min)
							{
								Tr_min=Tr_a;
								a_min=alpha;
							}
						}//*/
						//cout<<"Tr"<<count<<"="<<Tr<<", alpha="<<a_min<<endl;
				
						lam+=d[0]*a_min;
						mu+=d[1]*a_min;

						//cout<<"lam"<<count<<"="<<lam<<" ,mu"<<count<<"="<<mu<<endl;

						pn1[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam)+Dt*n[A_X]*mu;	
						pn1[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam)+Dt*n[A_Y]*mu;
						pn1[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam)+Dt*n[A_Z]*mu;

						qn1[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]-n[A_X]*lam));
						qn1[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]-n[A_Y]*lam));
						qn1[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]-n[A_Z]*lam));


						L=0.5/mi*(pn1[A_X]*pn1[A_X] + pn1[A_Y]*pn1[A_Y] + pn1[A_Z]*pn1[A_Z])-mi*G*(qn1[A_X]*nG[A_X] + qn1[A_Y]*nG[A_Y] + qn1[A_Z]*nG[A_Z]);

						Tr=(L-Ln)*(L-Ln);

						h=-1*( (qn1[A_X] - a[A_X])*n[A_X] + (qn1[A_Y] - a[A_Y])*n[A_Y] + (qn1[A_Z] - a[A_Z])*n[A_Z]);
						dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

						if(h+theta_h>0)		Tr+=0.5*r*(h+theta_h)*(h+theta_h);
						if(dhdt+theta_dh>0)	Tr+=0.5*r*(dhdt+theta_dh)*(dhdt+theta_dh);

						dL[0]=1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z])-G*Dt*Dt*(n[A_X]*nG[A_X]+n[A_Y]*nG[A_Y]+n[A_Z]*nG[A_Z]);
						dL[1]=1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

						dTr[0]=2*dL[0]*(L-Ln);
						dTr[1]=2*dL[1]*(L-Ln);

						dh[0]=-Dt*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
						dh[1]=0;
						d_dhdt[0]=-Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
						d_dhdt[1]=-Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);


						if(h+theta_h>0)	dTr[0]+=r*dh[0]*(h+theta_h);
						if(dhdt+theta_dh>0)
						{
							dTr[0]+=r*d_dhdt[0]*(dhdt+theta_dh);
							dTr[1]+=r*d_dhdt[1]*(dhdt+theta_dh);
						}

						E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);

						//cout<<"Tr"<<count<<"="<<Tr<<endl;
						//cout<<"h="<<h<<", dhdt="<<dhdt<<endl;
						//cout<<"theta_h="<<theta_h<<", theta_dhdt="<<theta_dh<<endl;
						//cout<<"pn1="<<pn1[0]<<", "<<pn1[1]<<", "<<pn1[2]<<endl;
						//cout<<"qn1="<<qn1[0]<<", "<<qn1[1]<<", "<<qn1[2]<<endl;
						if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				
						//if(count%100==0)
						{
							stringstream se;
							se<<"./E_QN/E_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fe(se.str(), ios::app);
							fe<<count<<","<<E<<endl;
							fe.close();

							stringstream ssp;
							ssp<<"./p_QN/p_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fp(ssp.str(), ios::app);
							fp<<count<<","<<pn1[A_X]<<","<<pn1[A_Y]<<","<<pn1[A_Z]<<endl;
							fp.close();
	
							stringstream ssq;
							ssq<<"./q_QN/q_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fq(ssq.str(), ios::app);
							fq<<count<<","<<qn1[A_X]<<","<<qn1[A_Y]<<","<<qn1[A_Z]<<endl;
							fq.close();

							stringstream sslam;
							sslam<<"./lam_QN/lam_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream flam(sslam.str(), ios::app);
							flam<<count<<","<<lam<<endl;
							flam.close();
	
							stringstream ssmu;
							ssmu<<"./mu_QN/mu_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fmu(ssmu.str(), ios::app);
							fmu<<count<<","<<mu<<endl;
							fmu.close();

							stringstream se_Tr;
							se_Tr<<"./Tr_QN/Tr_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fTr(se_Tr.str(), ios::app);
							fTr<<Tr<<endl;	
							fTr.close();

							stringstream se_h;
							se_h<<"./h_QN/h_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fh(se_h.str(), ios::app);
							fh<<h<<endl;	
							fh.close();

							stringstream se_dhdt;
							se_dhdt<<"./dhdt_QN/dhdt_t"<<t<<"cmin"<<count_min<<".csv";
							ofstream fdhdt(se_dhdt.str(), ios::app);
							fdhdt<<dhdt<<endl;	
							fdhdt.close();
						}

						if(E<ep && h<=0)	break;
						/*
						if(E>old_E && h<=0)
						{
							dp[A_X]=dp_k[A_X];	dp[A_Y]=dp_k[A_Y];	dp[A_X]=dp_k[A_Z];
							dq[A_X]=dq_k[A_X];	dq[A_Y]=dq_k[A_Y];	dq[A_X]=dq_k[A_Z];
							break;
						}//*/
						double old_E=E;

						double s[2]={lam-lam_k,	mu-mu_k};
						double y[2]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1]};

						double beta=y[0]*s[0]+y[1]*s[1];

						if(beta>0)
						{

							double sigma=(s[0]*B[0]+s[1]*B[2])*s[0]+(s[0]*B[2]+s[1]*B[3])*s[1];
							double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
							double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

							B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
							B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
							B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
							B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
							//cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;					
						}
					}
				double old_E_min=E_min;
				E_min=sqrt( (old_lam-lam)*(old_lam-lam) + (old_mu-mu)*(old_mu-mu)	);	
				if(E_min<ep_min*1000)	r*=4;
		
				if(h+theta_h>0)	theta_h+=h;
				if(dhdt+theta_dh>0)	theta_dh+=dhdt;
	
				cout<<"E_min"<<count_min<<"="<<E_min<<endl;
					
				//if(count_min%100==0)
				{
					stringstream se_min;
					se_min<<"./E_QP/E_min"<<t<<".csv";
					ofstream f_emin(se_min.str(), ios::app);
					f_emin<<E_min<<endl;
					f_emin.close();

					stringstream sh;
					sh<<"./h_QP/h"<<t<<".csv";
					ofstream fh(sh.str(), ios::app);
					fh<<h<<endl;
					fh.close();

					stringstream se_c;
					se_c<<"./count_QP/count"<<t<<".csv";
					ofstream fc(se_c.str(), ios::app);
					fc<<count_min<<","<<count<<endl;
					fc.close();
		
					stringstream se_p;
					se_p<<"./p_QP/p"<<t<<".csv";
					ofstream fp(se_p.str(), ios::app);
					fp<<count_min<<","<<pn1[A_X]<<","<<pn1[A_Y]<<","<<pn1[A_Z]<<endl;
					fp.close();

					stringstream se_q;
					se_q<<"./q_QP/q"<<t<<".csv";
					ofstream fq(se_q.str(), ios::app);
					fq<<count_min<<","<<qn1[A_X]<<","<<qn1[A_Y]<<","<<qn1[A_Z]<<endl;
					fq.close();

					stringstream se_lam;
					se_lam<<"./lam_QP/lam"<<t<<".csv";
					ofstream flam(se_lam.str(), ios::app);
					flam<<count_min<<","<<lam<<endl;
					flam.close();

					stringstream se_mu;
					se_mu<<"./mu_QP/mu"<<t<<".csv";
					ofstream fmu(se_mu.str(), ios::app);
					fmu<<count_min<<","<<mu<<endl;
					fmu.close();

					stringstream se_Tr;
					se_Tr<<"./Tr_QP/Tr"<<t<<".csv";
					ofstream fTr(se_Tr.str(), ios::app);
					fTr<<Tr<<endl;	
					fTr.close();

					stringstream sse_h;
					sse_h<<"./h_QP/h"<<t<<".csv";
					ofstream fsh(sse_h.str(), ios::app);
					fsh<<h<<endl;	
					fsh.close();

					stringstream se_dhdt;
					se_dhdt<<"./dhdt_QP/dhdt"<<t<<".csv";
					ofstream fdhdt(se_dhdt.str(), ios::app);
					fdhdt<<dhdt<<endl;	
					fdhdt.close();
				}
/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
				if(count_min>c_max)	break;
			}
			pn[A_X]=pn1[A_X];	pn[A_Y]=pn1[A_Y];	pn[A_Z]=pn1[A_Z];
			qn[A_X]=qn1[A_X];	qn[A_Y]=qn1[A_Y];	qn[A_Z]=qn1[A_Z];


//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
			cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

			cout<<"------------------------------OK"<<endl;
			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn[A_X]<<","<<pn[A_Y]<<","<<pn[A_Z]<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn[A_X]<<","<<qn[A_Y]<<","<<qn[A_Z]<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fdh("dhdt.csv", ios::app);
		fdh<<dhdt<<endl;
		fdh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


		double *qm=new double [DIMENSION];
		double *pm=new double [DIMENSION];

		pm[A_X]=pn[A_X];	pm[A_Y]=pn[A_Y];	pm[A_Z]=pn[A_Z];
		qm[A_X]=qn[A_X];	qm[A_Y]=qn[A_Y];	qm[A_Z]=qn[A_Z];

		PM_AVS(t_max,t,qm);
		MM_AVS(t_max,t,Dt,mi,qm, pm);
		
		delete[]	qm;
		delete[]	pm;

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<t<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();
	}
	
	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}


int main_hy0()	
{

	system("mkdir E");
	system("mkdir p");
	system("mkdir q");
	system("mkdir Tr");

	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=0.1;
	double nG[DIMENSION]={0,0,1};
	double n[DIMENSION]={0,0,1};
	double a[DIMENSION]={0,0,0};

	double p[DIMENSION]={1,1,1};
	double q[DIMENSION]={1,1,1};

	double old_p[DIMENSION]={0,0,0};
	double old_q[DIMENSION]={0,0,0};

	double Ln=0;
	double L=0;
	double dpL[DIMENSION]={0,0,0};
	double dqL[DIMENSION]={0,0,0};
	
	double Phi=0;	
	double dpPhi[DIMENSION]={0,0,0};
	double dqPhi[DIMENSION]={0,0,0};

	double h=0;		//h(q)での定義が正解
	double dhdt=0;		//h(q)での定義が正解
	double Tr=0;
	double dTr[6]={0,0,0,0,0,0};

	double r=0.01;
	double ep=1e-5;
	double ep_min=1e-5;
	double theta=0;
	double theta_dhdt=0;
	int c_max=2000;
	double E=1;
	double E_min=1;

	double *d=new double [3*Nx];	
	double *Nr=new double [3*Nx];
	double *B=new double [9*Nx*Nx];

	int t=0;
	int t_max=10000;

	double p0[DIMENSION]={0,0,0};
	double q0[DIMENSION]={0,0,5};
	double pn[DIMENSION]={0,0,0};
	double qn[DIMENSION]={0,0,0};


	ofstream fs("intial_data.csv");
	fs<<"r, "<<r<<endl;
	fs<<"Dt, "<<Dt<<endl;
	fs<<"mi, "<<mi<<endl;
	fs<<"p0, "<<p0[A_X]<<","<<p0[A_Y]<<","<<p0[A_Z]<<endl;
	fs<<"q0, "<<q0[A_X]<<","<<q0[A_Y]<<","<<q0[A_Z]<<endl;
	fs<<"t_max, "<<t_max<<endl;
	fs<<"ep_min, "<<ep_min<<endl;
	fs<<"ep, "<<ep<<endl;
	fs.close();



	cout<<"Calculation Start"<<endl;


	pn[A_X]=p0[A_X];
	pn[A_Y]=p0[A_Y];
	pn[A_Z]=p0[A_Z];
	qn[A_X]=q0[A_X];
	qn[A_Y]=q0[A_Y];
	qn[A_Z]=q0[A_Z];

	double E0=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);


	while(t<t_max)
	{
		t++;
		clock_t t0=clock();

		double pn_1[DIMENSION]={pn[A_X], pn[A_Y], pn[A_Z]};
		double qn_1[DIMENSION]={qn[A_X], qn[A_Y], qn[A_Z]};
		Ln=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])-mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);


		qn[A_X]+=Dt/mi*pn_1[A_X];	qn[A_Y]+=Dt/mi*pn_1[A_Y];	qn[A_Z]+=Dt/mi*pn_1[A_Z];
		pn[A_X]+=-Dt*mi*G*nG[A_X];	pn[A_Y]+=-Dt*mi*G*nG[A_Y];	pn[A_Z]+=-Dt*mi*G*nG[A_Z];

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();
			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*( (qn[A_X]-a[A_X])*n[A_X] + (qn[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]-a[A_Z])*n[A_Z]);

		if(h>0)
		{
			p[A_X]=pn[A_X];	p[A_Y]=pn[A_Y];	p[A_Z]=pn[A_Z];
			q[A_X]=qn[A_X];	q[A_Y]=qn[A_Y];	q[A_Z]=qn[A_Z];

			cout<<"接触"<<t<<" ,QP_method"<<endl;


			E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./p"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./q"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
	
					stringstream se_Tr;
					se_Tr<<"./Tr"<<t<<".csv";
					ofstream init_Tr(se_Tr.str(), ios::trunc);
					init_Tr.close();
				}

				if(count_min>c_max && h<=0)	break;

				old_p[A_X]=p[A_X];	old_p[A_Y]=p[A_Y];	old_p[A_Z]=p[A_Z];
				old_q[A_X]=q[A_X];	old_q[A_Y]=q[A_Y];	old_q[A_Z]=q[A_Z];

				for(int i=0;i<3*Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<3*Nx;j++)
					{
						if(j==i)	B[i*3*Nx+j]=1;
						else
						{
							B[i*3*Nx+j]=0;
						}
					}
				}


				L=0.5/mi*(p[A_X]*p[A_X] + p[A_Y]*p[A_Y] + p[A_Z]*p[A_Z])-mi*G*(q[A_X]*nG[A_X] + q[A_Y]*nG[A_Y] + q[A_Z]*nG[A_Z]);


				h=-1*( (q[A_X] - a[A_X])*n[A_X] + (q[A_Y] - a[A_Y])*n[A_Y] + (q[A_Z] - a[A_Z])*n[A_Z]);
				dhdt=-1*(p[A_X]*n[A_X] +p[A_Y]*n[A_Y]+p[A_Z]*n[A_Z]);
				Phi=0.5*r*(h+theta)*(h+theta);

				Tr=(L-Ln)*(L-Ln);
				if(h+theta>0)	Tr+=Phi;
				if(dhdt+theta_dhdt>0)	Tr+=0.5*r*(dhdt+theta_dhdt)*(dhdt+theta_dhdt);
				dpL[A_X]=1/mi*p[A_X];	dpL[A_Y]=1/mi*p[A_Y];	dpL[A_Z]=1/mi*p[A_Z];
				dqL[A_X]=-mi*G*nG[A_X];				dqL[A_Y]=-mi*G*nG[A_Y];				dqL[A_Z]=-mi*G*nG[A_Z];

				dTr[0]=2*dpL[A_X]*(L-Ln);
				dTr[1]=2*dpL[A_Y]*(L-Ln);
				dTr[2]=2*dpL[A_Z]*(L-Ln);
				dTr[3]=2*dqL[A_X]*(L-Ln);
				dTr[4]=2*dqL[A_Y]*(L-Ln);
				dTr[5]=2*dqL[A_Z]*(L-Ln);

				if(h+theta>0)
				{
					dqPhi[A_X]=-r*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*n[A_Z]*(h+theta);		

					dTr[3]+=dqPhi[A_X];
					dTr[4]+=dqPhi[A_Y];
					dTr[5]+=dqPhi[A_Z];
				}
				if(dhdt+theta_dhdt>0)
				{
					dTr[0]+=-r*n[A_X]*(dhdt+theta_dhdt);
					dTr[1]+=-r*n[A_X]*(dhdt+theta_dhdt);
					dTr[2]+=-r*n[A_X]*(dhdt+theta_dhdt);
				}
				E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);
				int count=0;
 				if(E<ep)	break;
				while(E>ep)
				{
					count++;
					if(count==1)
					{
						stringstream se;
						se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_e(se.str(), ios::trunc);
						init_e.close();

						stringstream ssp;
						ssp<<"./p/p_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_p(ssp.str(), ios::trunc);
						init_p.close();
	
						stringstream ssq;
						ssq<<"./q/q_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_q(ssq.str(), ios::trunc);
						init_q.close();
						stringstream se_Tr;
						se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream fTr(se_Tr.str(), ios::trunc);
						fTr.close();
					}
					if(count>c_max)	break;

					double p_k[DIMENSION]={p[A_X],p[A_Y], p[A_Z]};
					double q_k[DIMENSION]={q[A_X],q[A_Y], q[A_Z]};
					double dTr_k[6]={dTr[0], dTr[1], dTr[2], dTr[3], dTr[4], dTr[5]};

					for(int i=0;i<3*Nx;i++)	Nr[i]=dTr[i];
					if(count>1)
					{
/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
						cout<<"B=";
						for(int i=0;i<3*Nx;i++)
						{
							for(int j=0;j<3*Nx;j++)
							{
								cout<<B[i*3*Nx+j]<<", ";
							}
							cout<<endl;
						}*/

						gauss(B,Nr,3*Nx);
					}

					for(int i=0;i<3*Nx;i++)	d[i]=-1*Nr[i];

	//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

					double Tr_min=Tr;
					double a_min=0.1;//1e-4;

					for(int i=0;i<10000;i++)
					{
						double alpha=(i+1)*1e-4;
						double p_a[DIMENSION]={p[A_X]+d[0]*alpha, p[A_Y]+d[1]*alpha, p[A_Z]+d[2]*alpha};
						double q_a[DIMENSION]={q[A_X]+d[3]*alpha, q[A_Y]+d[4]*alpha, q[A_Z]+d[5]*alpha};
	
						double L_a=0.5/mi*(p_a[A_X]*p_a[A_X] + p_a[A_Y]*p_a[A_Y] + p_a[A_Z]*p_a[A_Z])-mi*G*(q_a[A_X]*nG[A_X] + q_a[A_Y]*nG[A_Y] + q_a[A_Z]*nG[A_Z]);

						double h_a=-1*( (q_a[A_X] - a[A_X])*n[A_X] + (q_a[A_Y] - a[A_Y])*n[A_Y] + (q_a[A_Z] - a[A_Z])*n[A_Z]);
						double dhdt_a=-1*( p_a[A_X]*n[A_X]+p_a[A_Y]*n[A_Y]+p_a[A_Z]*n[A_Z]);
						double Phi_a=0.5*r*(h_a+theta)*(h_a+theta);

						double Tr_a=(L_a-Ln)*(L_a-Ln);
						if(h_a+theta>0)	Tr_a+=Phi_a;
						if(dhdt_a+theta_dhdt>0)	Tr_a+=r*0.5*(dhdt_a+theta_dhdt)*(dhdt_a+theta_dhdt);

						if(Tr_a<Tr_min)
						{
							Tr_min=Tr_a;
							a_min=alpha;
						}
					}//*/
//					cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
					p[A_X]+=d[0]*a_min;
					p[A_Y]+=d[1]*a_min;
					p[A_Z]+=d[2]*a_min;
					q[A_X]+=d[3]*a_min;
					q[A_Y]+=d[4]*a_min;
					q[A_Z]+=d[5]*a_min;
//					cout<<"dp"<<count<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//					cout<<"dq"<<count<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;



					L=0.5/mi*(p[A_X]*p[A_X] + p[A_Y]*p[A_Y] + p[A_Z]*p[A_Z])-mi*G*(q[A_X]*nG[A_X] + q[A_Y]*nG[A_Y] + q[A_Z]*nG[A_Z]);


					h=-1*( (q[A_X] - a[A_X])*n[A_X] + (q[A_Y] - a[A_Y])*n[A_Y] + (q[A_Z] - a[A_Z])*n[A_Z]);
					dhdt=-1*(p[A_X]*n[A_X] +p[A_Y]*n[A_Y]+p[A_Z]*n[A_Z]);
					Phi=0.5*r*(h+theta)*(h+theta);

					Tr=(L-Ln)*(L-Ln);
					if(h+theta>0)	Tr+=Phi;
					if(dhdt+theta_dhdt>0)	Tr+=r*0.5*(dhdt+theta_dhdt)*(dhdt+theta_dhdt);

					dpL[A_X]=1/mi*p[A_X];	dpL[A_Y]=1/mi*p[A_Y];	dpL[A_Z]=1/mi*p[A_Z];
					dqL[A_X]=-mi*G*nG[A_X];				dqL[A_Y]=-mi*G*nG[A_Y];				dqL[A_Z]=-mi*G*nG[A_Z];

					dTr[0]=2*dpL[A_X]*(L-Ln);
					dTr[1]=2*dpL[A_Y]*(L-Ln);
					dTr[2]=2*dpL[A_Z]*(L-Ln);
					dTr[3]=2*dqL[A_X]*(L-Ln);
					dTr[4]=2*dqL[A_Y]*(L-Ln);
					dTr[5]=2*dqL[A_Z]*(L-Ln);

					if(h+theta>0)
					{
						dqPhi[A_X]=-r*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*n[A_Z]*(h+theta);		

						dTr[3]+=dqPhi[A_X];
						dTr[4]+=dqPhi[A_Y];
						dTr[5]+=dqPhi[A_Z];
					}
					if(dhdt+theta_dhdt>0)
					{
						dTr[0]+=-r*n[A_X]*(dhdt+theta_dhdt);
						dTr[1]+=-r*n[A_X]*(dhdt+theta_dhdt);
						dTr[2]+=-r*n[A_X]*(dhdt+theta_dhdt);
					}

					E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);



					double old_E=E;
					if(count%500==0)cout<<"E"<<count<<"="<<E<<endl;
						
					stringstream se;
					se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fe(se.str(), ios::app);
					fe<<count<<","<<E<<endl;
					fe.close();

					stringstream ssp;
					ssp<<"./p/p_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fp(ssp.str(), ios::app);
					fp<<count<<","<<p[A_X]<<","<<p[A_Y]<<","<<p[A_Z]<<endl;
					fp.close();
	
					stringstream ssq;
					ssq<<"./q/q_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fq(ssq.str(), ios::app);
					fq<<count<<","<<q[A_X]<<","<<q[A_Y]<<","<<q[A_Z]<<endl;
					fq.close();

					stringstream se_Tr;
					se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fTr(se_Tr.str(), ios::app);
					fTr<<Tr<<endl;	
					fTr.close();

					if(E<ep && h<=0)	break;
					
					
					old_E=E;

					double s[6]={p[A_X]-p_k[A_X],p[A_Y]-p_k[A_Y], p[A_Z]-p_k[A_Z], q[A_X]-q_k[A_X], q[A_Y]-q_k[A_Y], q[A_Z]-q_k[A_Z]};
					double y[6]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1], dTr[2]-dTr_k[2], dTr[3]-dTr_k[3], dTr[4]-dTr_k[4], dTr[5]-dTr_k[5]};



					double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3]+y[4]*s[4]+y[5]*s[5];
					double sigma=0;
					for(int j=0;j<6;j++)
					{
						double sB_j=0;
						for(int k=0;k<6;k++)
						{
							sB_j+=s[k]*B[k*6+j];
						}
						sigma+=sB_j*s[j];
					}


					//if(count%100==0)	cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

					if(beta>0)//(beta>=0.2*sigma)
					{
						double bs[6]={0,0,0,0,0,0};
						for(int j=0;j<6;j++)
						{
							double sB_j=0;
							for(int k=0;k<6;k++)
							{
								sB_j+=B[j*6+k]*s[k];
							}
							bs[j]=sB_j;
						}						

						double bss[36]={bs[0]*s[0],bs[0]*s[1],bs[0]*s[2],bs[0]*s[3],bs[0]*s[4],bs[0]*s[5],
										bs[1]*s[0],bs[1]*s[1],bs[1]*s[2],bs[1]*s[3],bs[1]*s[4],bs[1]*s[5],
										bs[2]*s[0],bs[2]*s[1],bs[2]*s[2],bs[2]*s[3],bs[2]*s[4],bs[2]*s[5],
										bs[3]*s[0],bs[3]*s[1],bs[3]*s[2],bs[3]*s[3],bs[3]*s[4],bs[3]*s[5],
										bs[4]*s[0],bs[4]*s[1],bs[4]*s[2],bs[4]*s[3],bs[4]*s[4],bs[4]*s[5],
										bs[5]*s[0],bs[5]*s[1],bs[5]*s[2],bs[5]*s[3],bs[5]*s[4],bs[5]*s[5]};

						double bssb[36];
						for(int j=0;j<6;j++)
						{
							for(int k=0;k<6;k++)
							{
								bssb[j*6+k]=0;
								for(int i=0;i<6;i++)
								{
									bssb[j*6+k]+=bss[j*6+i]*B[i*6+k];
								}
							}
						}
								
						for(int j=0;j<6;j++)
						{
							for(int k=0;k<6;k++)	B[j*6+k]+=1/beta*y[j]*y[k]-1/sigma*bssb[j*6+k];
						}
					}
				}
				double old_E_min=E_min;
				E_min=sqrt( (old_p[A_X]-p[A_X])*(old_p[A_X]-p[A_X]) + (old_p[A_Y]-p[A_Y])*(old_p[A_Y]-p[A_Y]) + (old_p[A_Z]-p[A_Z])*(old_p[A_Z]-p[A_Z]) +
							(old_q[A_X]-q[A_X])*(old_q[A_X]-q[A_X]) + (old_q[A_Y]-q[A_Y])*(old_q[A_Y]-q[A_Y]) + (old_q[A_Z]-q[A_Z])*(old_q[A_Z]-q[A_Z]));	
				if(E_min<ep_min*1000)	r*=4;
				if(h+theta>0)	theta+=h;
				if(dhdt+theta_dhdt>0)	theta_dhdt+=dhdt;
	
				if(count_min%100==0)	cout<<"E_min"<<count_min<<"="<<E_min<<endl;
					
				
				stringstream se_min;
				se_min<<"./E_min"<<t<<".csv";
				ofstream f_emin(se_min.str(), ios::app);
				f_emin<<E_min<<endl;
				f_emin.close();

				stringstream sh;
				sh<<"./h"<<t<<".csv";
				ofstream fh(sh.str(), ios::app);
				fh<<h<<endl;
				fh.close();

				stringstream se_c;
				se_c<<"./count"<<t<<".csv";
				ofstream fc(se_c.str(), ios::app);
				fc<<count_min<<","<<count<<endl;
				fc.close();

		
				stringstream se_p;
				se_p<<"./p"<<t<<".csv";
				ofstream fp(se_p.str(), ios::app);
				fp<<count_min<<","<<p[A_X]<<","<<p[A_Y]<<","<<p[A_Z]<<endl;
				fp.close();

				stringstream se_q;
				se_q<<"./q"<<t<<".csv";
				ofstream fq(se_q.str(), ios::app);
				fq<<count_min<<","<<q[A_X]<<","<<q[A_Y]<<","<<q[A_Z]<<endl;
				fq.close();

				stringstream se_Tr;
				se_Tr<<"./Tr"<<t<<".csv";
				ofstream fTr(se_Tr.str(), ios::app);
				fTr<<Tr<<endl;	
				fTr.close();

/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
				if(count_min>c_max)	break;
			}

			double mu=0;
			double fc=h;

			//v0　計算
			//int nv=0;
			//cout<<"fc="<<h<<endl;
			//if(h>-ep && h<ep)
			//{
			//	mu=-n[A_Z]/dqL[A_Z];
			//}
			//cout<<"fc="<<fc<<", mu="<<mu<<endl;

			qn[A_X]=q[A_X];		qn[A_Y]=q[A_Y];		qn[A_Z]=q[A_Z];
			pn[A_X]=p[A_X];		pn[A_Y]=p[A_Y];		pn[A_Z]=p[A_Z];
	
			//pn[A_X]-=-mu*(-1/mi*(p[A_X]*n[A_X]+p[A_Y]*n[A_Y]+p[A_Z]*n[A_Z]));
			//pn[A_Y]-=-mu*(-1/mi*(p[A_X]*n[A_X]+p[A_Y]*n[A_Y]+p[A_Z]*n[A_Z]));
			//pn[A_Z]-=-mu*(-1/mi*(p[A_X]*n[A_X]+p[A_Y]*n[A_Y]+p[A_Z]*n[A_Z]));



//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
			cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

			cout<<"------------------------------OK"<<endl;
			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn[A_X]<<","<<pn[A_Y]<<","<<pn[A_Z]<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn[A_X]<<","<<qn[A_Y]<<","<<qn[A_Z]<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


		double *qm=new double [DIMENSION];
		double *pm=new double [DIMENSION];

		pm[A_X]=pn[A_X];	pm[A_Y]=pn[A_Y];	pm[A_Z]=pn[A_Z];
		qm[A_X]=qn[A_X];	qm[A_Y]=qn[A_Y];	qm[A_Z]=qn[A_Z];

		PM_AVS(t_max,t,qm);
		MM_AVS(t_max,t,Dt,mi,qm, pm);
		
		delete[]	qm;
		delete[]	pm;

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<t<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();
	}
	
	//v0　計算
//	int nv=0;
//	cout<<"fc="<<c0<<", "<<c1<<endl;
	

	//double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);
	//double v1=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	//cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	//cout<<"c="<<c0<<endl;
	//cout<<"v="<<v0<<endl;
	//cout<<"fx="<<fx<<endl;
	//cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/

	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}





int main_minus()	
{

	system("mkdir E_QN");
	system("mkdir p_QN");
	system("mkdir q_QN");
	system("mkdir lam_QN");
	system("mkdir mu_QN");
	system("mkdir Tr_QN");

	system("mkdir E_QP");
	system("mkdir count_QP");
	system("mkdir h_QP");
	system("mkdir p_QP");
	system("mkdir q_QP");
	system("mkdir lam_QP");
	system("mkdir mu_QP");
	system("mkdir Tr_QP");


	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=1e-3;
	double nG[DIMENSION]={0,0,1};
	double n[DIMENSION]={0,0,1};
	double a[DIMENSION]={0,0,0};

	double lam=0;
	double mu=0;

	double old_lam=0;
	double old_mu=0;

	double Ln=0;
	double L=0;
	double dL[2]={0,0};
	
	double h=0;		//h(q)での定義が正解
	double dh[2]={0,0};
	double dhdt=0;
	double d_dhdt[2]={0,0};

	double Tr=0;
	double dTr[2]={0,0};


	double r=0.01;
	double ep=1e-10;
	double ep_min=1e-10;
	double theta_h=0;
	double theta_dh=0;
	int c_max=2000;
	double E=1;
	double E_min=1;

	double *d=new double [Nx];	
	double *Nr=new double [Nx];
	double *B=new double [Nx*Nx];

	int t=0;
	int t_max=10000;

	double p0[DIMENSION]={0,0,0};
	double q0[DIMENSION]={0,0,3};
	double pn[DIMENSION]={0,0,0};
	double qn[DIMENSION]={0,0,0};
	double pn1[DIMENSION]={0,0,0};
	double qn1[DIMENSION]={0,0,0};
	double pn1_2[DIMENSION]={0,0,0};


	ofstream fs("intial_data.csv");
	fs<<"r, "<<r<<endl;
	fs<<"Dt, "<<Dt<<endl;
	fs<<"mi, "<<mi<<endl;
	fs<<"p0, "<<p0[A_X]<<","<<p0[A_Y]<<","<<p0[A_Z]<<endl;
	fs<<"q0, "<<q0[A_X]<<","<<q0[A_Y]<<","<<q0[A_Z]<<endl;
	fs<<"t_max, "<<t_max<<endl;
	fs<<"ep_min, "<<ep_min<<endl;
	fs<<"ep, "<<ep<<endl;
	fs.close();


	cout<<"Calculation Start"<<endl;


	pn[A_X]=p0[A_X];
	pn[A_Y]=p0[A_Y];
	pn[A_Z]=p0[A_Z];
	qn[A_X]=q0[A_X];
	qn[A_Y]=q0[A_Y];
	qn[A_Z]=q0[A_Z];

	double E0=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);


	while(t<t_max)
	{
		t++;
		clock_t t0=clock();

		double pn_1[DIMENSION]={pn[A_X], pn[A_Y], pn[A_Z]};
		double qn_1[DIMENSION]={qn[A_X], qn[A_Y], qn[A_Z]};

		pn1_2[A_X]=pn[A_X]-Dt*mi*G*nG[A_X];	pn1_2[A_Y]=pn[A_Y]-Dt*mi*G*nG[A_Y];	pn1_2[A_Z]=pn[A_Z]-Dt*mi*G*nG[A_Z];
		qn[A_X]+=Dt/mi*pn1_2[A_X];	qn[A_Y]+=Dt/mi*pn1_2[A_Y];	qn[A_Z]+=Dt/mi*pn1_2[A_Z];
		pn[A_X]+=-Dt*mi*G*nG[A_X];	pn[A_Y]+=-Dt*mi*G*nG[A_Y];	pn[A_Z]+=-Dt*mi*G*nG[A_Z];

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();
			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*( (qn[A_X]-a[A_X])*n[A_X] + (qn[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]-a[A_Z])*n[A_Z]);
		if(h>0)
		{
			pn[A_X]=pn_1[A_X];	pn[A_Y]=pn_1[A_Y];	pn[A_Z]=pn_1[A_Z];
			qn[A_X]=qn_1[A_X];	qn[A_Y]=qn_1[A_Y];	qn[A_Z]=qn_1[A_Z];//

			cout<<"接触"<<t<<" ,QP_method"<<endl;
			
			theta_dh=0;
			theta_h=0;
			r=0.1;
			lam=1;
			mu=1;

			Ln=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])-mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);

			E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_QP/E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count_QP/count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./p_QP/p"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./q_QP/q"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
	
					stringstream se_lam;
					se_lam<<"./lam_QP/lam"<<t<<".csv";
					ofstream init_lam(se_lam.str(), ios::trunc);
					init_lam.close();

					stringstream se_mu;
					se_mu<<"./mu_QP/mu"<<t<<".csv";
					ofstream init_mu(se_mu.str(), ios::trunc);
					init_mu.close();

					stringstream se_Tr;
					se_Tr<<"./Tr_QP/Tr"<<t<<".csv";
					ofstream init_Tr(se_Tr.str(), ios::trunc);
					init_Tr.close();
				}

				if(count_min>c_max && h<=0)	break;

				old_lam=lam;
				old_mu=mu;

				for(int i=0;i<Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<Nx;j++)
					{
						if(j==i)	B[i*Nx+j]=1;
						else
						{
							B[i*Nx+j]=0;
						}
					}
				}

				pn1[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam)-Dt*n[A_X]*mu;	
				pn1[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam)-Dt*n[A_Y]*mu;
				pn1[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam)-Dt*n[A_Z]*mu;

				qn1[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam));
				qn1[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam));
				qn1[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam));


				L=0.5/mi*(pn1[A_X]*pn1[A_X] + pn1[A_Y]*pn1[A_Y] + pn1[A_Z]*pn1[A_Z])-mi*G*(qn1[A_X]*nG[A_X] + qn1[A_Y]*nG[A_Y] + qn1[A_Z]*nG[A_Z]);

				Tr=(L-Ln)*(L-Ln);

				h=-1*( (qn1[A_X] - a[A_X])*n[A_X] + (qn1[A_Y] - a[A_Y])*n[A_Y] + (qn1[A_Z] - a[A_Z])*n[A_Z]);
				dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

				if(h+theta_h>0)	Tr+=0.5*r*(h+theta_h)*(h+theta_h);
				if(dhdt+theta_dh>0)	Tr+=0.5*r*(dhdt+theta_dh)*(dhdt+theta_dh);

				dL[0]=-1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z])+G*Dt*Dt*(n[A_X]*nG[A_X]+n[A_Y]*nG[A_Y]+n[A_Z]*nG[A_Z]);
				dL[1]=-1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

				dTr[0]=2*dL[0]*(L-Ln);
				dTr[1]=2*dL[1]*(L-Ln);

				dh[0]=Dt*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
				dh[1]=0;
				d_dhdt[0]=Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
				d_dhdt[1]=Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);

				if(h+theta_h>0)	dTr[0]+=r*dh[0]*(h+theta_h);	
				if(dhdt+theta_dh>0)
				{
					dTr[0]+=r*d_dhdt[0]*(dhdt+theta_dh);
					dTr[1]+=r*d_dhdt[1]*(dhdt+theta_dh);
				}

				E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);

				int count=0;
 				if(E<ep)	break;
				while(E>ep)
				{
					count++;
					if(count==1)
					{
						stringstream se;
						se<<"./E_QN/E_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_e(se.str(), ios::trunc);
						init_e.close();

						stringstream ssp;
						ssp<<"./p_QN/p_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_p(ssp.str(), ios::trunc);
						init_p.close();
	
						stringstream ssq;
						ssq<<"./q_QN/q_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_q(ssq.str(), ios::trunc);
						init_q.close();

						stringstream sslam;
						sslam<<"./lam_QN/lam_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_lam(sslam.str(), ios::trunc);
						init_lam.close();
	
						stringstream ssmu;
						ssmu<<"./mu_QN/mu_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_mu(ssmu.str(), ios::trunc);
						init_mu.close();

						stringstream se_Tr;
						se_Tr<<"./Tr_QN/Tr_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream fTr(se_Tr.str(), ios::trunc);
						fTr.close();
					}
					if(count>c_max)	break;

					double lam_k=lam;
					double mu_k=mu;
					double dTr_k[2]={dTr[0],dTr[1]};

					Nr[0]=dTr[0];	Nr[1]=dTr[1];
					if(count>1)
					{
/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
						cout<<"B=";
						for(int i=0;i<3*Nx;i++)
						{
							for(int j=0;j<3*Nx;j++)
							{
								cout<<B[i*3*Nx+j]<<", ";
							}
							cout<<endl;
						}*/

						gauss(B,Nr,Nx);
					}
					d[0]=-1*Nr[0];	d[1]=-1*Nr[1];

	//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

					double Tr_min=Tr;
					double a_min=0.1;//1e-4;
					double pn1_a[DIMENSION]={0,0,0};
					double qn1_a[DIMENSION]={0,0,0};
					double L_a=0,	h_a=0,	dhdt_a=0,	Tr_a=0;


					for(int i=0;i<10000;i++)
					{
						double alpha=(i+1)*1e-4;
						double lam_a=lam+d[0]*alpha;
						double mu_a=mu+d[1]*alpha;
						
						pn1_a[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam_a)-Dt*n[A_X]*mu_a;	
						pn1_a[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam_a)-Dt*n[A_Y]*mu_a;
						pn1_a[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam_a)-Dt*n[A_Z]*mu_a;

						qn1_a[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam_a));
						qn1_a[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam_a));
						qn1_a[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam_a));

						L_a=0.5/mi*(pn1_a[A_X]*pn1_a[A_X] + pn1_a[A_Y]*pn1_a[A_Y] + pn1_a[A_Z]*pn1_a[A_Z])-mi*G*(qn1_a[A_X]*nG[A_X] + qn1_a[A_Y]*nG[A_Y] + qn1_a[A_Z]*nG[A_Z]);
						Tr_a=(L_a-Ln)*(L_a-Ln);

						h_a=-1*( (qn1_a[A_X] - a[A_X])*n[A_X] + (qn1_a[A_Y] - a[A_Y])*n[A_Y] + (qn1_a[A_Z] - a[A_Z])*n[A_Z]);
						dhdt_a=-1/mi*(n[A_X]*pn1_a[A_X]+n[A_Y]*pn1_a[A_Y]+n[A_Z]*pn1_a[A_Z]);

						if(h_a+theta_h>0)	Tr_a+=0.5*r*(h_a+theta_h)*(h_a+theta_h);

						if(dhdt_a+theta_dh>0)	Tr_a+=0.5*r*(dhdt_a+theta_dh)*(dhdt_a+theta_dh);

						if(Tr_a<Tr_min)
						{
							Tr_min=Tr_a;
							a_min=alpha;
						}
					}//*/
//					cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
					lam+=d[0]*a_min;
					mu+=d[1]*a_min;

//					cout<<"dp"<<count<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//					cout<<"dq"<<count<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
	
					pn1[A_X]=pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam)-Dt*n[A_X]*mu;	
					pn1[A_Y]=pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam)-Dt*n[A_Y]*mu;
					pn1[A_Z]=pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam)-Dt*n[A_Z]*mu;

					qn1[A_X]=qn[A_X]+Dt/mi*(pn[A_X]-Dt*(mi*G*nG[A_X]+n[A_X]*lam));
					qn1[A_Y]=qn[A_Y]+Dt/mi*(pn[A_Y]-Dt*(mi*G*nG[A_Y]+n[A_Y]*lam));
					qn1[A_Z]=qn[A_Z]+Dt/mi*(pn[A_Z]-Dt*(mi*G*nG[A_Z]+n[A_Z]*lam));


					L=0.5/mi*(pn1[A_X]*pn1[A_X] + pn1[A_Y]*pn1[A_Y] + pn1[A_Z]*pn1[A_Z])-mi*G*(qn1[A_X]*nG[A_X] + qn1[A_Y]*nG[A_Y] + qn1[A_Z]*nG[A_Z]);

					Tr=(L-Ln)*(L-Ln);

					h=-1*( (qn1[A_X] - a[A_X])*n[A_X] + (qn1[A_Y] - a[A_Y])*n[A_Y] + (qn1[A_Z] - a[A_Z])*n[A_Z]);
					dhdt=-1/mi*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

					if(h+theta_h>0)	Tr+=0.5*r*(h+theta_h)*(h+theta_h);
					if(dhdt+theta_dh>0)	Tr+=0.5*r*(dhdt+theta_dh)*(dhdt+theta_dh);

					dL[0]=-1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z])+G*Dt*Dt*(n[A_X]*nG[A_X]+n[A_Y]*nG[A_Y]+n[A_Z]*nG[A_Z]);
					dL[1]=-1/mi*Dt*(n[A_X]*pn1[A_X]+n[A_Y]*pn1[A_Y]+n[A_Z]*pn1[A_Z]);

					dTr[0]=2*dL[0]*(L-Ln);
					dTr[1]=2*dL[1]*(L-Ln);

					dh[0]=Dt*Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
					dh[1]=0;
					d_dhdt[0]=Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);
					d_dhdt[1]=Dt/mi*(n[A_X]*n[A_X]+n[A_Y]*n[A_Y]+n[A_Z]*n[A_Z]);

					if(h+theta_h>0)	dTr[0]+=r*dh[0]*(h+theta_h);	
					if(dhdt+theta_dh>0)
					{
						dTr[0]+=r*d_dhdt[0]*(dhdt+theta_dh);
						dTr[1]+=r*d_dhdt[1]*(dhdt+theta_dh);
					}


					E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);

					if(count%500==0)	cout<<"E"<<count<<"="<<E<<endl;
				

					stringstream se;
					se<<"./E_QN/E_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fe(se.str(), ios::app);
					fe<<count<<","<<E<<endl;
					fe.close();

					stringstream ssp;
					ssp<<"./p_QN/p_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fp(ssp.str(), ios::app);
					fp<<count<<","<<pn1[A_X]<<","<<pn1[A_Y]<<","<<pn1[A_Z]<<endl;
					fp.close();
	
					stringstream ssq;
					ssq<<"./q_QN/q_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fq(ssq.str(), ios::app);
					fq<<count<<","<<qn1[A_X]<<","<<qn1[A_Y]<<","<<qn1[A_Z]<<endl;
					fq.close();

					stringstream sslam;
					sslam<<"./lam_QN/lam_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream flam(sslam.str(), ios::app);
					flam<<count<<","<<lam<<endl;
					flam.close();
	
					stringstream ssmu;
					ssmu<<"./mu_QN/mu_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fmu(ssmu.str(), ios::app);
					fmu<<count<<","<<mu<<endl;
					fmu.close();

					stringstream se_Tr;
					se_Tr<<"./Tr_QN/Tr_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fTr(se_Tr.str(), ios::app);
					fTr<<Tr<<endl;	
					fTr.close();

					if(E<ep && h<=0)	break;
					
					/*
					if(E>old_E && h<=0)
					{
						dp[A_X]=dp_k[A_X];	dp[A_Y]=dp_k[A_Y];	dp[A_X]=dp_k[A_Z];
						dq[A_X]=dq_k[A_X];	dq[A_Y]=dq_k[A_Y];	dq[A_X]=dq_k[A_Z];
						break;
					}//*/
					double old_E_min=E_min;

					double s[2]={lam-lam_k,	mu-mu_k};
					double y[2]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1]};

					double beta=y[0]*s[0]+y[1]*s[1];

					if(beta>0)
					{

						double sigma=(s[0]*B[0]+s[1]*B[2])*s[0]+(s[0]*B[2]+s[1]*B[3])*s[1];
						double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
						double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

						B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
						B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
						B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
						B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
						//if(count%1000==0)	cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;					
					}
					
				}
				double old_E_min=E_min;
				E_min=sqrt( (old_lam-lam)*(old_lam-lam) + (old_mu-mu)*(old_mu-mu)	);	
				if(E_min<ep_min*1000)	r*=4;
		
				if(h+theta_h>0)	theta_h+=h;
				if(dhdt+theta_dh>0)	theta_dh+=dhdt;
	
				if(count_min%500==0)	cout<<"E_min"<<count_min<<"="<<E_min<<endl;
					
				
				stringstream se_min;
				se_min<<"./E_QP/E_min"<<t<<".csv";
				ofstream f_emin(se_min.str(), ios::app);
				f_emin<<E_min<<endl;
				f_emin.close();

				stringstream sh;
				sh<<"./h_QP/h"<<t<<".csv";
				ofstream fh(sh.str(), ios::app);
				fh<<h<<endl;
				fh.close();

				stringstream se_c;
				se_c<<"./count_QP/count"<<t<<".csv";
				ofstream fc(se_c.str(), ios::app);
				fc<<count_min<<","<<count<<endl;
				fc.close();

		
				stringstream se_p;
				se_p<<"./p_QP/p"<<t<<".csv";
				ofstream fp(se_p.str(), ios::app);
				fp<<count_min<<","<<pn1[A_X]<<","<<pn1[A_Y]<<","<<pn1[A_Z]<<endl;
				fp.close();

				stringstream se_q;
				se_q<<"./q_QP/q"<<t<<".csv";
				ofstream fq(se_q.str(), ios::app);
				fq<<count_min<<","<<qn1[A_X]<<","<<qn1[A_Y]<<","<<qn1[A_Z]<<endl;
				fq.close();

				stringstream se_lam;
				se_lam<<"./lam_QP/lam"<<t<<".csv";
				ofstream flam(se_lam.str(), ios::app);
				flam<<count_min<<","<<lam<<endl;
				flam.close();

				stringstream se_mu;
				se_mu<<"./mu_QP/mu"<<t<<".csv";
				ofstream fmu(se_mu.str(), ios::app);
				fmu<<count_min<<","<<mu<<endl;
				fmu.close();

				stringstream se_Tr;
				se_Tr<<"./Tr_QP/Tr"<<t<<".csv";
				ofstream fTr(se_Tr.str(), ios::app);
				fTr<<Tr<<endl;	
				fTr.close();

/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
			}
			pn[A_X]=pn1[A_X];	pn[A_Y]=pn1[A_Y];	pn[A_Z]=pn1[A_Z];
			qn[A_X]=qn1[A_X];	qn[A_Y]=qn1[A_Y];	qn[A_Z]=qn1[A_Z];


//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
			cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

			cout<<"------------------------------OK"<<endl;
			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn[A_X]<<","<<pn[A_Y]<<","<<pn[A_Z]<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn[A_X]<<","<<qn[A_Y]<<","<<qn[A_Z]<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


		double *qm=new double [DIMENSION];
		double *pm=new double [DIMENSION];

		pm[A_X]=pn[A_X];	pm[A_Y]=pn[A_Y];	pm[A_Z]=pn[A_Z];
		qm[A_X]=qn[A_X];	qm[A_Y]=qn[A_Y];	qm[A_Z]=qn[A_Z];

		PM_AVS(t_max,t,qm);
		MM_AVS(t_max,t,Dt,mi,qm, pm);
		
		delete[]	qm;
		delete[]	pm;

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<t<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();
	}
	
	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}








int main_hy0_dpdq()	
{
	system("mkdir E");
	system("mkdir p");
	system("mkdir q");
	system("mkdir Tr");

	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=1e-3;
	double nG[DIMENSION]={0,0,1};
	double n[DIMENSION]={0,0,1};
	double a[DIMENSION]={0,0,0};

	double dp[DIMENSION]={1,1,1};
	double dq[DIMENSION]={1,1,1};

	double old_dp[DIMENSION]={0,0,0};
	double old_dq[DIMENSION]={0,0,0};

	double Ln=0;
	double L=0;
	double dpL[DIMENSION]={0,0,0};
	double dqL[DIMENSION]={0,0,0};
	
	double Phi=0;	
	double dpPhi[DIMENSION]={0,0,0};
	double dqPhi[DIMENSION]={0,0,0};

	double h=0;
	double Tr=0;
	double dTr[6]={0,0,0,0,0,0};

	double r=0.01;
	double ep=1e-10;
	double ep_min=1e-10;
	double theta=0;
	int c_max=5000;
	double E=1;
	double E_min=1;

	double *d=new double [3*Nx];	
	double *Nr=new double [3*Nx];
	double *B=new double [9*Nx*Nx];

	int t=0;
	int t_max=10000;

	double p0[DIMENSION]={0,0,0};
	double q0[DIMENSION]={0,0,3};
	double pn[DIMENSION]={0,0,0};
	double qn[DIMENSION]={0,0,0};

	ofstream fs("intial_data.csv");
	fs<<"r, "<<r<<endl;
	fs<<"Dt, "<<Dt<<endl;
	fs<<"mi, "<<mi<<endl;
	fs<<"p0, "<<p0[A_X]<<","<<p0[A_Y]<<","<<p0[A_Z]<<endl;
	fs<<"q0, "<<q0[A_X]<<","<<q0[A_Y]<<","<<q0[A_Z]<<endl;
	fs<<"t_max, "<<t_max<<endl;
	fs<<"ep_min, "<<ep_min<<endl;
	fs<<"ep, "<<ep<<endl;
	fs.close();



	cout<<"Calculation Start"<<endl;


	pn[A_X]=p0[A_X];
	pn[A_Y]=p0[A_Y];
	pn[A_Z]=p0[A_Z];
	qn[A_X]=q0[A_X];
	qn[A_Y]=q0[A_Y];
	qn[A_Z]=q0[A_Z];

	double E0=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);


	while(t<t_max)
	{
		t++;
		clock_t t0=clock();

		double pn_1[DIMENSION]={pn[A_X], pn[A_Y], pn[A_Z]};
		double qn_1[DIMENSION]={qn[A_X], qn[A_Y], qn[A_Z]};


		qn[A_X]+=Dt/mi*pn_1[A_X];	qn[A_Y]+=Dt/mi*pn_1[A_Y];	qn[A_Z]+=Dt/mi*pn_1[A_Z];
		pn[A_X]+=-Dt*mi*G*nG[A_X];	pn[A_Y]+=-Dt*mi*G*nG[A_Y];	pn[A_Z]+=-Dt*mi*G*nG[A_Z];

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();
			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*( (qn[A_X]-a[A_X])*n[A_X] + (qn[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]-a[A_Z])*n[A_Z]);
		if(h>0)
		{
			pn[A_X]=pn_1[A_X];	pn[A_Y]=pn_1[A_Y];	pn[A_Z]=pn_1[A_Z];
			qn[A_X]=qn_1[A_X];	qn[A_Y]=qn_1[A_Y];	qn[A_Z]=qn_1[A_Z];

			cout<<"接触"<<t<<" ,QP_method"<<endl;
			dp[A_X]=pn[A_X];	dp[A_Y]=pn[A_Y];	dp[A_Z]=pn[A_Z];
			dq[A_X]=qn[A_X];	dq[A_Y]=qn[A_Y];	dq[A_Z]=qn[A_Z];

			Ln=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])-mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);

			E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./p"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./q"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
	
					stringstream se_Tr;
					se_Tr<<"./Tr"<<t<<".csv";
					ofstream init_Tr(se_Tr.str(), ios::trunc);
					init_Tr.close();
				}

				if(count_min>c_max && h<=0)	break;

				old_dp[A_X]=dp[A_X];	old_dp[A_Y]=dp[A_Y];	old_dp[A_Z]=dp[A_Z];
				old_dq[A_X]=dq[A_X];	old_dq[A_Y]=dq[A_Y];	old_dq[A_Z]=dq[A_Z];

				for(int i=0;i<3*Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<3*Nx;j++)
					{
						if(j==i)	B[i*3*Nx+j]=1;
						else
						{
							B[i*3*Nx+j]=0;
						}
					}
				}


				L=0.5/mi*((pn[A_X]+Dt*dp[A_X])*(pn[A_X]+Dt*dp[A_X]) + (pn[A_Y]+Dt*dp[A_Y])*(pn[A_Y]+Dt*dp[A_Y]) + (pn[A_Z]+Dt*dp[A_Z])*(pn[A_Z]+Dt*dp[A_Z]))
					-mi*G*((qn[A_X]+Dt*dq[A_X])*nG[A_X] + (qn[A_Y]+Dt*dq[A_Y])*nG[A_Y] + (qn[A_Z]+Dt*dq[A_Z])*nG[A_Z]);


				h=-1*( (qn[A_X]+Dt*dq[A_X] - a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq[A_Y] - a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq[A_Z] - a[A_Z])*n[A_Z]);
				Phi=0.5*r*(h+theta)*(h+theta);

				Tr=(L-Ln)*(L-Ln);
				if(h+theta>0)	Tr+=Phi;

				dpL[A_X]=1/mi*Dt*(pn[A_X]+Dt*dp[A_X]);	dpL[A_Y]=1/mi*Dt*(pn[A_Y]+Dt*dp[A_Y]);	dpL[A_Z]=1/mi*Dt*(pn[A_Z]+Dt*dp[A_Z]);
				dqL[A_X]=-mi*G*Dt*nG[A_X];				dqL[A_Y]=-mi*G*Dt*nG[A_Y];				dqL[A_Z]=-mi*G*Dt*nG[A_Z];

				dTr[0]=2*dpL[A_X]*(L-Ln);
				dTr[1]=2*dpL[A_Y]*(L-Ln);
				dTr[2]=2*dpL[A_Z]*(L-Ln);
				dTr[3]=2*dqL[A_X]*(L-Ln);
				dTr[4]=2*dqL[A_Y]*(L-Ln);
				dTr[5]=2*dqL[A_Z]*(L-Ln);

				if(h>-ep && h<ep)
				{
					dqPhi[A_X]=-r*Dt*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*Dt*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*Dt*n[A_Z]*(h+theta);		

					if(dqPhi[A_X]>0)	dTr[3]=+dqPhi[A_X];
					if(dqPhi[A_Y]>0)	dTr[4]=+dqPhi[A_Y];
					if(dqPhi[A_Z]>0)	dTr[5]=+dqPhi[A_Z];
				}				
				else if(h+theta>0)
				{
					dqPhi[A_X]=-r*Dt*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*Dt*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*Dt*n[A_Z]*(h+theta);		

					dTr[3]=+dqPhi[A_X];
					dTr[4]=+dqPhi[A_Y];
					dTr[5]=+dqPhi[A_Z];
				}
				

				E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);

				int count=0;
 				if(E<ep)	break;
				while(E>ep)
				{
					count++;
					if(count==1)
					{
						stringstream se;
						se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_e(se.str(), ios::trunc);
						init_e.close();

						stringstream ssp;
						ssp<<"./p/p_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_p(ssp.str(), ios::trunc);
						init_p.close();
	
						stringstream ssq;
						ssq<<"./q/q_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_q(ssq.str(), ios::trunc);
						init_q.close();
						stringstream se_Tr;
						se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream fTr(se_Tr.str(), ios::trunc);
						fTr.close();
					}
					if(count>c_max)	break;

					double dp_k[DIMENSION]={dp[A_X],dp[A_Y], dp[A_Z]};
					double dq_k[DIMENSION]={dq[A_X],dq[A_Y], dq[A_Z]};
					double dTr_k[6]={dTr[0], dTr[1], dTr[2], dTr[3], dTr[4], dTr[5]};

					for(int i=0;i<3*Nx;i++)	Nr[i]=dTr[i];
					if(count>1)
					{
/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
						cout<<"B=";
						for(int i=0;i<3*Nx;i++)
						{
							for(int j=0;j<3*Nx;j++)
							{
								cout<<B[i*3*Nx+j]<<", ";
							}
							cout<<endl;
						}*/

						gauss(B,Nr,3*Nx);
					}

					for(int i=0;i<3*Nx;i++)	d[i]=-1*Nr[i];

	//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

					double Tr_min=Tr;
					double a_min=0.1;//1e-4;

					for(int i=0;i<10000;i++)
					{
						double alpha=(i+1)*1e-4;
						double dp_a[DIMENSION]={dp[A_X]+d[0]*alpha, dp[A_Y]+d[1]*alpha, dp[A_Z]+d[2]*alpha};
						double dq_a[DIMENSION]={dq[A_X]+d[3]*alpha, dq[A_Y]+d[4]*alpha, dq[A_Z]+d[5]*alpha};
	
						double	L_a=0.5/mi*((pn[A_X]+Dt*dp_a[A_X])*(pn[A_X]+Dt*dp_a[A_X]) + (pn[A_Y]+Dt*dp_a[A_Y])*(pn[A_Y]+Dt*dp_a[A_Y]) + (pn[A_Z]+Dt*dp_a[A_Z])*(pn[A_Z]+Dt*dp_a[A_Z]))
									-mi*G*((qn[A_X]+Dt*dq_a[A_X])*nG[A_X] + (qn[A_Y]+Dt*dq_a[A_Y])*nG[A_Y] + (qn[A_Z]+Dt*dq_a[A_Z])*nG[A_Z]);

						double h_a=-1*( (qn[A_X]+Dt*dq_a[A_X] - a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq_a[A_Y] - a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq_a[A_Z] - a[A_Z])*n[A_Z]);
						double Phi_a=0.5*r*(h_a+theta)*(h_a+theta);

						double Tr_a=(L_a-Ln)*(L_a-Ln);
						if(h_a+theta>0)	Tr_a+=Phi_a;


						if(Tr_a<Tr_min)
						{
							Tr_min=Tr_a;
							a_min=alpha;
						}
					}//*/
//					cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
					dp[A_X]+=d[0]*a_min;
					dp[A_Y]+=d[1]*a_min;
					dp[A_Z]+=d[2]*a_min;
					dq[A_X]+=d[3]*a_min;
					dq[A_Y]+=d[4]*a_min;
					dq[A_Z]+=d[5]*a_min;
//					cout<<"dp"<<count<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//					cout<<"dq"<<count<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;

					L=0.5/mi*((pn[A_X]+Dt*dp[A_X])*(pn[A_X]+Dt*dp[A_X]) + (pn[A_Y]+Dt*dp[A_Y])*(pn[A_Y]+Dt*dp[A_Y]) + (pn[A_Z]+Dt*dp[A_Z])*(pn[A_Z]+Dt*dp[A_Z]))
					-mi*G*((qn[A_X]+Dt*dq[A_X])*nG[A_X] + (qn[A_Y]+Dt*dq[A_Y])*nG[A_Y] + (qn[A_Z]+Dt*dq[A_Z])*nG[A_Z]);


					h=-1*( (qn[A_X]+Dt*dq[A_X] - a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq[A_Y] - a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq[A_Z] - a[A_Z])*n[A_Z]);
					Phi=0.5*r*(h+theta)*(h+theta);

					Tr=(L-Ln)*(L-Ln);
					if(h+theta>0)	Tr+=Phi;

					dpL[A_X]=1/mi*Dt*(pn[A_X]+Dt*dp[A_X]);	dpL[A_Y]=1/mi*Dt*(pn[A_Y]+Dt*dp[A_Y]);	dpL[A_Z]=1/mi*Dt*(pn[A_Z]+Dt*dp[A_Z]);
					dqL[A_X]=-mi*G*Dt*nG[A_X];				dqL[A_Y]=-mi*G*Dt*nG[A_Y];				dqL[A_Z]=-mi*G*Dt*nG[A_Z];

					dTr[0]=2*dpL[A_X]*(L-Ln);
					dTr[1]=2*dpL[A_Y]*(L-Ln);
					dTr[2]=2*dpL[A_Z]*(L-Ln);
					dTr[3]=2*dqL[A_X]*(L-Ln);
					dTr[4]=2*dqL[A_Y]*(L-Ln);
					dTr[5]=2*dqL[A_Z]*(L-Ln);

					if(h>-ep && h<ep)
					{
						dqPhi[A_X]=-r*Dt*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*Dt*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*Dt*n[A_Z]*(h+theta);		

						if(dqPhi[A_X]>0)	dTr[3]=+dqPhi[A_X];
						if(dqPhi[A_Y]>0)	dTr[4]=+dqPhi[A_Y];
						if(dqPhi[A_Z]>0)	dTr[5]=+dqPhi[A_Z];
					}
					else if(h+theta>0)
					{
						dqPhi[A_X]=-r*Dt*n[A_X]*(h+theta);		dqPhi[A_Y]=-r*Dt*n[A_Y]*(h+theta);		dqPhi[A_Z]=-r*Dt*n[A_Z]*(h+theta);		

						dTr[3]=+dqPhi[A_X];
						dTr[4]=+dqPhi[A_Y];
						dTr[5]=+dqPhi[A_Z];
					}

					E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);



					double old_E=E;

					E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);

					if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
						
					stringstream se;
					se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fe(se.str(), ios::app);
					fe<<count<<","<<E<<endl;
					fe.close();

					stringstream ssp;
					ssp<<"./p/p_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fp(ssp.str(), ios::app);
					fp<<count<<","<<dp[A_X]<<","<<dp[A_Y]<<","<<dp[A_Z]<<endl;
					fp.close();
	
					stringstream ssq;
					ssq<<"./q/q_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fq(ssq.str(), ios::app);
					fq<<count<<","<<dq[A_X]<<","<<dq[A_Y]<<","<<dq[A_Z]<<endl;
					fq.close();

					stringstream se_Tr;
					se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fTr(se_Tr.str(), ios::app);
					fTr<<Tr<<endl;	
					fTr.close();

					if(E<ep && h<=0)	break;
					
					/*
					if(E>old_E && h<=0)
					{
						dp[A_X]=dp_k[A_X];	dp[A_Y]=dp_k[A_Y];	dp[A_X]=dp_k[A_Z];
						dq[A_X]=dq_k[A_X];	dq[A_Y]=dq_k[A_Y];	dq[A_X]=dq_k[A_Z];
						break;
					}//*/
					double old_E_min=E_min;

					double s[6]={dp[A_X]-dp_k[A_X],dp[A_Y]-dp_k[A_Y], dp[A_Z]-dp_k[A_Z], dq[A_X]-dq_k[A_X], dq[A_Y]-dq_k[A_Y], dq[A_Z]-dq_k[A_Z]};
					double y[6]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1], dTr[2]-dTr_k[2], dTr[3]-dTr_k[3], dTr[4]-dTr_k[4], dTr[5]-dTr_k[5]};



					double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3]+y[4]*s[4]+y[5]*s[5];
					double sigma=0;
					for(int j=0;j<6;j++)
					{
						double sB_j=0;
						for(int k=0;k<6;k++)
						{
							sB_j+=s[k]*B[k*6+j];
						}
						sigma+=sB_j*s[j];
					}


					if(count%100==0)	cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

					if(beta>0)//(beta>=0.2*sigma)
					{
						double bs[6]={0,0,0,0,0,0};
						for(int j=0;j<6;j++)
						{
							double sB_j=0;
							for(int k=0;k<6;k++)
							{
								sB_j+=B[j*6+k]*s[k];
							}
							bs[j]=sB_j;
						}						

						double bss[36]={bs[0]*s[0],bs[0]*s[1],bs[0]*s[2],bs[0]*s[3],bs[0]*s[4],bs[0]*s[5],
										bs[1]*s[0],bs[1]*s[1],bs[1]*s[2],bs[1]*s[3],bs[1]*s[4],bs[1]*s[5],
										bs[2]*s[0],bs[2]*s[1],bs[2]*s[2],bs[2]*s[3],bs[2]*s[4],bs[2]*s[5],
										bs[3]*s[0],bs[3]*s[1],bs[3]*s[2],bs[3]*s[3],bs[3]*s[4],bs[3]*s[5],
										bs[4]*s[0],bs[4]*s[1],bs[4]*s[2],bs[4]*s[3],bs[4]*s[4],bs[4]*s[5],
										bs[5]*s[0],bs[5]*s[1],bs[5]*s[2],bs[5]*s[3],bs[5]*s[4],bs[5]*s[5]};

						double bssb[36];
						for(int j=0;j<6;j++)
						{
							for(int k=0;k<6;k++)
							{
								bssb[j*6+k]=0;
								for(int i=0;i<6;i++)
								{
									bssb[j*6+k]+=bss[j*6+i]*B[i*6+k];
								}
							}
						}
								
						for(int j=0;j<6;j++)
						{
							for(int k=0;k<6;k++)	B[j*6+k]+=1/beta*y[j]*y[k]-1/sigma*bssb[j*6+k];
						}
					}
				}
				double old_E_min=E_min;
				E_min=sqrt( (old_dp[A_X]-dp[A_X])*(old_dp[A_X]-dp[A_X]) + (old_dp[A_Y]-dp[A_Y])*(old_dp[A_Y]-dp[A_Y]) + (old_dp[A_Z]-dp[A_Z])*(old_dp[A_Z]-dp[A_Z]) +
							(old_dq[A_X]-dq[A_X])*(old_dq[A_X]-dq[A_X]) + (old_dq[A_Y]-dq[A_Y])*(old_dq[A_Y]-dq[A_Y]) + (old_dq[A_Z]-dq[A_Z])*(old_dq[A_Z]-dq[A_Z]));	
				if(E_min<ep_min*1000)	r*=4;
				if(h+theta>0)	theta+=h;
	
				if(count_min%100==0)	cout<<"E_min"<<count_min<<"="<<E_min<<endl;
					
				
				stringstream se_min;
				se_min<<"./E_min"<<t<<".csv";
				ofstream f_emin(se_min.str(), ios::app);
				f_emin<<E_min<<endl;
				f_emin.close();

				stringstream sh;
				sh<<"./h"<<t<<".csv";
				ofstream fh(sh.str(), ios::app);
				fh<<h<<endl;
				fh.close();

				stringstream se_c;
				se_c<<"./count"<<t<<".csv";
				ofstream fc(se_c.str(), ios::app);
				fc<<count_min<<","<<count<<endl;
				fc.close();

		
				stringstream se_p;
				se_p<<"./p"<<t<<".csv";
				ofstream fp(se_p.str(), ios::app);
				fp<<count_min<<","<<dp[A_X]<<","<<dp[A_Y]<<","<<dp[A_Z]<<endl;
				fp.close();

				stringstream se_q;
				se_q<<"./q"<<t<<".csv";
				ofstream fq(se_q.str(), ios::app);
				fq<<count_min<<","<<dq[A_X]<<","<<dq[A_Y]<<","<<dq[A_Z]<<endl;
				fq.close();

				stringstream se_Tr;
				se_Tr<<"./Tr"<<t<<".csv";
				ofstream fTr(se_Tr.str(), ios::app);
				fTr<<Tr<<endl;	
				fTr.close();

/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
			}

			double mu=0;
			double fc=h;

			//v0　計算
			int nv=0;
			cout<<"fc="<<h<<endl;
			if(h>-ep && h<ep)
			{
				double dqhz=-Dt*n[A_Z];
				mu=-dqhz/dqL[A_Z];
			}
			cout<<"fc="<<fc<<", mu="<<mu<<endl;

			qn[A_X]+=Dt*dq[A_X];		qn[A_Y]+=Dt*dq[A_Y];		qn[A_Z]+=Dt*dq[A_Z];
			pn[A_X]+=Dt*dp[A_X];		pn[A_Y]+=Dt*dp[A_Y];		pn[A_Z]+=Dt*dp[A_Z];

	/*		pn[A_X]-=-Dt/mi*(pn[A_X]+Dt*dp[A_X])*n[A_X]*mu;
			pn[A_Y]-=-Dt/mi*(pn[A_Y]+Dt*dp[A_Y])*n[A_Y]*mu;
			pn[A_Z]-=-Dt/mi*(pn[A_Z]+Dt*dp[A_Z])*n[A_Z]*mu;*/




//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
			cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

			cout<<"------------------------------OK"<<endl;
			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn[A_X]<<","<<pn[A_Y]<<","<<pn[A_Z]<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn[A_X]<<","<<qn[A_Y]<<","<<qn[A_Z]<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


		double *qm=new double [DIMENSION];
		double *pm=new double [DIMENSION];

		pm[A_X]=pn[A_X];	pm[A_Y]=pn[A_Y];	pm[A_Z]=pn[A_Z];
		qm[A_X]=qn[A_X];	qm[A_Y]=qn[A_Y];	qm[A_Z]=qn[A_Z];

		PM_AVS(t_max,t,qm);
		MM_AVS(t_max,t,Dt,mi,qm, pm);
		
		delete[]	qm;
		delete[]	pm;

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<t<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();
	}
	
	//v0　計算
//	int nv=0;
//	cout<<"fc="<<c0<<", "<<c1<<endl;
	

	//double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);
	//double v1=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	//cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	//cout<<"c="<<c0<<endl;
	//cout<<"v="<<v0<<endl;
	//cout<<"fx="<<fx<<endl;
	//cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/

	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}

int main_n1()	
{
	system("mkdir E");
	system("mkdir dp");
	system("mkdir dq");
	system("mkdir Tr");

	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=1e-3;

	double dp=0;
	double dq=0;

	double old_dp=0;
	double old_dq=0;

	double Ln=0;
	double L=0;
	double dpL=0;
	double dqL=0;
	
	double Phi=0;	
	double dpPhi=0;
	double dqPhi=0;

	double h=0;
	double Tr=0;
	double dTr[2]={0,0};

	double r=0.01;
	double ep=1e-10;
	double ep_min=1e-10;
	double theta=0;
	int c_max=5000;

	double *d=new double [Nx];	
	double *Nr=new double [Nx];
	double *B=new double [Nx*Nx];

	int t=0;
	int t_max=10000;

	double p0=0;
	double q0=5;
	double pn=0;
	double qn=0;



	ofstream fs("intial_data.csv");
	fs<<"r, "<<r<<endl;
	fs<<"Dt, "<<Dt<<endl;
	fs<<"mi, "<<mi<<endl;
	fs<<"p0, "<<p0<<endl;
	fs<<"q0, "<<q0<<endl;
	fs<<"t_max, "<<t_max<<endl;
	fs<<"ep_min, "<<ep_min<<endl;
	fs<<"ep, "<<ep<<endl;
	fs.close();



	cout<<"Calculation Start"<<endl;


	pn=p0;
	qn=q0;

	double E0=0.5/mi*pn*pn+mi*G*qn;


	while(t<t_max)
	{
		t++;
		clock_t t0=clock();

		double pn_1=pn;
		double qn_1=qn;

		qn+=Dt/mi*pn_1;
		pn+=-Dt*mi*G;	

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();
			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*qn;
		if(h>0)
		{
			pn=pn_1;
			qn=qn_1;

			cout<<"接触"<<t<<" ,QP_method"<<endl;
			dp=1;
			dq=1;

			double E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./dp"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./dq"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
	
					stringstream se_Tr;
					se_Tr<<"./Tr"<<t<<".csv";
					ofstream init_Tr(se_Tr.str(), ios::trunc);
					init_Tr.close();
				}

				if(count_min>c_max && h<=0)	break;

				old_dp=dp;
				old_dq=dq;

				for(int i=0;i<Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<Nx;j++)
					{
						if(j==i)	B[i*Nx+j]=1;
						else
						{
							B[i*Nx+j]=0;
						}
					}
				}

				Ln=0.5/mi*pn*pn-mi*G*qn;

				L=0.5/mi*(pn+Dt*dp)*(pn+Dt*dp)-mi*G*(qn+Dt*dq);


				h=-1*(qn+Dt*dq);
				Phi=0.5*r*(h+theta)*(h+theta);

				Tr=(L-Ln)*(L-Ln);
				if(h+theta>0)	Tr+=Phi;

				dpL=Dt/mi*(pn+Dt*dp);
				dqL=-mi*G*Dt;

				dTr[0]=2*dpL*(L-Ln);
				dTr[1]=2*dqL*(L-Ln);
			

				if(h+theta>0)
				{
					dqPhi=r*(-Dt)*(h+theta);

					dTr[0]=+dpPhi;
					dTr[1]=+dqPhi;
				}
				double E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);

				int count=0;
 				if(E<ep)	break;
				while(E>ep)
				{
					count++;
					if(count==1)
					{
						stringstream se;
						se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_e(se.str(), ios::trunc);
						init_e.close();

						stringstream ssp;
						ssp<<"./dp/dp_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_p(ssp.str(), ios::trunc);
						init_p.close();
	
						stringstream ssq;
						ssq<<"./dq/dq_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_q(ssq.str(), ios::trunc);
						init_q.close();
						stringstream se_Tr;
						se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream fTr(se_Tr.str(), ios::trunc);
						fTr.close();
					}
					if(count>c_max)	break;

					double dp_k=dp;
					double dq_k=dq;
					double dTr_k[2]={dTr[0], dTr[1]};

					for(int i=0;i<Nx;i++)	Nr[i]=dTr[i];
					if(count>1)
					{
/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
						cout<<"B=";
						for(int i=0;i<3*Nx;i++)
						{
							for(int j=0;j<3*Nx;j++)
							{
								cout<<B[i*3*Nx+j]<<", ";
							}
							cout<<endl;
						}*/

						gauss(B,Nr,Nx);
					}

					for(int i=0;i<Nx;i++)	d[i]=-1*Nr[i];

	//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

					double Tr_min=Tr;
					double a_min=0.1;//1e-4;

					for(int i=0;i<10000;i++)
					{
						double alpha=(i+1)*1e-4;
						double dp_a=dp+d[0]*alpha;
						double dq_a=dq+d[1]*alpha;
	
						double L_a=0.5/mi*(pn+Dt*dp_a)*(pn+Dt*dp_a)-mi*G*(qn+Dt*dq_a);

						double h_a=-1*( qn+Dt*dq_a);
						double Phi_a=0.5*r*(h_a+theta)*(h_a+theta);

						double Tr_a=(L_a-Ln)*(L_a-Ln);
						if(h_a+theta>0)	Tr_a+=Phi_a;


						if(Tr_a<Tr_min)
						{
							Tr_min=Tr_a;
							a_min=alpha;
						}
					}//*/
//					cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
					dp+=d[0]*a_min;
					dq+=d[1]*a_min;
//					cout<<"dp"<<count<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//					cout<<"dq"<<count<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;


					Ln=0.5/mi*pn*pn-mi*G*qn;
					L=0.5/mi*(pn+Dt*dp)*(pn+Dt*dp)-mi*G*(qn+Dt*dq);


					h=-1*(qn+Dt*dq);
					Phi=0.5*r*(h+theta)*(h+theta);

					Tr=(L-Ln)*(L-Ln);
					if(h+theta>0)	Tr+=Phi;

					dpL=Dt/mi*(pn+Dt*dp);
					dqL=-mi*G*Dt;

					dTr[0]=2*dpL*(L-Ln);
					dTr[1]=2*dqL*(L-Ln);
			

					if(h+theta>0)
					{
						dqPhi=r*(-Dt)*(h+theta);

						dTr[0]=+dpPhi;
						dTr[1]=+dqPhi;
					}
					double old_E=E;
					double E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]);



					if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
						
					stringstream se;
					se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fe(se.str(), ios::app);
					fe<<count<<","<<E<<endl;
					fe.close();

					stringstream ssp;
					ssp<<"./dp/dp_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fp(ssp.str(), ios::app);
					fp<<count<<","<<dp<<endl;
					fp.close();
	
					stringstream ssq;
					ssq<<"./dq/dq_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fq(ssq.str(), ios::app);
					fq<<count<<","<<dq<<endl;
					fq.close();

					stringstream se_Tr;
					se_Tr<<"./Tr/Tr_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fTr(se_Tr.str(), ios::app);
					fTr<<Tr<<endl;	
					fTr.close();

					if(E<ep)	break;
					
					/*
					if(E>old_E && h<=0)
					{
						dp[A_X]=dp_k[A_X];	dp[A_Y]=dp_k[A_Y];	dp[A_X]=dp_k[A_Z];
						dq[A_X]=dq_k[A_X];	dq[A_Y]=dq_k[A_Y];	dq[A_X]=dq_k[A_Z];
						break;
					}//*/
					double old_E_min=E_min;

					double s[2]={dp-dp_k,dq-dq_k};
					double y[2]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1]};



					double beta=y[0]*s[0]+y[1]*s[1];
					double sigma=0;
					for(int j=0;j<Nx;j++)
					{
						double sB_j=0;
						for(int k=0;k<Nx;k++)
						{
							sB_j+=s[k]*B[k*Nx+j];
						}
						sigma+=sB_j*s[j];
					}


					if(count%100==0)	cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

					if(beta>0)//(beta>=0.2*sigma)
					{
						double bs[6]={0,0,0,0,0,0};
						for(int j=0;j<Nx;j++)
						{
							double sB_j=0;
							for(int k=0;k<Nx;k++)
							{
								sB_j+=B[j*Nx+k]*s[k];
							}
							bs[j]=sB_j;
						}						

						double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

						double bssb[4];
						for(int j=0;j<Nx;j++)
						{
							for(int k=0;k<Nx;k++)
							{
								bssb[j*Nx+k]=0;
								for(int i=0;i<Nx;i++)
								{
									bssb[j*Nx+k]+=bss[j*Nx+i]*B[i*Nx+k];
								}
							}
						}
								
						for(int j=0;j<Nx;j++)
						{
							for(int k=0;k<Nx;k++)	B[j*Nx+k]+=1/beta*y[j]*y[k]-1/sigma*bssb[j*Nx+k];
						}
					}
				}
				double old_E_min=E_min;
				E_min=sqrt( (old_dp-dp)*(old_dp-dp)+(old_dp-dp)*(old_dq-dq));
				if(E_min<ep_min*1000)	r*=4;
				if(h+theta>0)	theta+=h;
	
				if(count_min%100==0)	cout<<"E_min"<<count_min<<"="<<E_min<<endl;
					
				
				stringstream se_min;
				se_min<<"./E_min"<<t<<".csv";
				ofstream f_emin(se_min.str(), ios::app);
				f_emin<<E_min<<endl;
				f_emin.close();

				stringstream sh;
				sh<<"./h"<<t<<".csv";
				ofstream fh(sh.str(), ios::app);
				fh<<h<<endl;
				fh.close();

				stringstream se_c;
				se_c<<"./count"<<t<<".csv";
				ofstream fc(se_c.str(), ios::app);
				fc<<count_min<<","<<count<<endl;
				fc.close();

		
				stringstream se_p;
				se_p<<"./dp"<<t<<".csv";
				ofstream fp(se_p.str(), ios::app);
				fp<<count_min<<","<<dp<<endl;
				fp.close();

				stringstream se_q;
				se_q<<"./dq"<<t<<".csv";
				ofstream fq(se_q.str(), ios::app);
				fq<<count_min<<","<<dq<<endl;
				fq.close();

				stringstream se_Tr;
				se_Tr<<"./Tr"<<t<<".csv";
				ofstream fTr(se_Tr.str(), ios::app);
				fTr<<Tr<<endl;	
				fTr.close();

/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
			}
//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			qn+=Dt*dq;
			pn+=Dt*dp;

			cout<<"pn"<<t<<"="<<pn<<endl;
			cout<<"qn"<<t<<"="<<qn<<endl;

			cout<<"------------------------------OK"<<endl;
			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*pn*pn+mi*G*qn;
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


/*		double *q;
		double *p;

		p=pn;
		q=qn;

		PM_AVS(t_max,t,q);
		MM_AVS(t_max,t,Dt,mi,q, p);
		
		delete[]	q;
		delete[]	p;*/

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<t<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();
	}
	
	//v0　計算
//	int nv=0;
//	cout<<"fc="<<c0<<", "<<c1<<endl;
	

	//double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);
	//double v1=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	//cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	//cout<<"c="<<c0<<endl;
	//cout<<"v="<<v0<<endl;
	//cout<<"fx="<<fx<<endl;
	//cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/

	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}


int h_q()	
{
	system("mkdir E");
	system("mkdir dp");
	system("mkdir dq");

	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double mi=1.0;
	double Dt=1e-3;
	double nG[DIMENSION]={0,0,1};
	double n[DIMENSION]={0,0,1};
	double a[DIMENSION]={0,0,0};

	double dp[DIMENSION]={1,1,1};
	double dq[DIMENSION]={1,1,1};

	double old_dp[DIMENSION]={0,0,0};
	double old_dq[DIMENSION]={0,0,0};

	double h=0;
	double Tr=0;
	double dTr[6]={0,0,0,0,0,0};

	double r=0.01;
	double ep=1e-10;
	double ep_min=1e-5;
	double theta=0;

	double *d=new double [3*Nx];	
	double *Nr=new double [3*Nx];
	double *B=new double [9*Nx*Nx];

	int t=0;
	int t_max=10000;

	double p0[DIMENSION]={0,0,0};
	double q0[DIMENSION]={0,0,5};
	double pn[DIMENSION]={0,0,0};
	double qn[DIMENSION]={0,0,0};

	pn[A_X]=p0[A_X];
	pn[A_Y]=p0[A_Y];
	pn[A_Z]=p0[A_Z];
	qn[A_X]=q0[A_X];
	qn[A_Y]=q0[A_Y];
	qn[A_Z]=q0[A_Z];

	double E0=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])+mi*G*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);

	while(t<t_max)
	{
		t++;
		clock_t t0=clock();


		double pn_1[DIMENSION]={pn[A_X], pn[A_Y], pn[A_Z]};
		double qn_1[DIMENSION]={qn[A_X], qn[A_Y], qn[A_Z]};

		qn[A_X]+=Dt/mi*pn_1[A_X];	qn[A_Y]+=Dt/mi*pn_1[A_Y];	qn[A_Z]+=Dt/mi*pn_1[A_Z];
		pn[A_X]+=-Dt*mi*G*nG[A_X];	pn[A_Y]+=-Dt*mi*G*nG[A_Y];	pn[A_Z]+=-Dt*mi*G*nG[A_Z];

		if(t==1)
		{
			ofstream init_p("p.csv", ios::trunc);
			init_p.close();
			ofstream init_q("q.csv", ios::trunc);
			init_q.close();
			ofstream fh("h.csv", ios::trunc);
			fh.close();
			ofstream init_t("time_log.csv", ios::trunc);
			init_t.close();
			ofstream init_En("Energy.csv", ios::trunc);
			init_En.close();

			ofstream init_cmin("cmin.csv", ios::trunc);
			init_cmin.close();
		}

//		cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
//		cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;

		h=-1*( (qn[A_X]-a[A_X])*n[A_X] + (qn[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]-a[A_Z])*n[A_Z]);
		if(h>0)
		{
			pn[A_X]=pn_1[A_X];	pn[A_Y]=pn_1[A_Y];	pn[A_Z]=pn_1[A_Z];
			qn[A_X]=qn_1[A_X];	qn[A_Y]=qn_1[A_Y];	qn[A_Z]=qn_1[A_Z];

			cout<<"QP_method"<<endl;
			dp[A_X]=1;	dp[A_Y]=1;	dp[A_Z]=1;
			dq[A_X]=1;	dq[A_Y]=1;	dq[A_Z]=1;
			double E_min=1;
			int count_min=0;
			while(E_min>ep_min)
			{
				count_min++;
				if(count_min==1)
				{
					stringstream se_min;
					se_min<<"./E_min"<<t<<".csv";
					ofstream init_emin(se_min.str(), ios::trunc);
					init_emin.close();
					
					stringstream se_c;
					se_c<<"./count"<<t<<".csv";
					ofstream init_c(se_c.str(), ios::trunc);
					init_c.close();
		
					stringstream se_p;
					se_p<<"./dp"<<t<<".csv";
					ofstream init_p(se_p.str(), ios::trunc);
					init_p.close();

					stringstream se_q;
					se_q<<"./dq"<<t<<".csv";
					ofstream init_q(se_q.str(), ios::trunc);
					init_q.close();
				}

				if(count_min>500 && h<=0)	break;

				old_dp[A_X]=dp[A_X];
				old_dp[A_Y]=dp[A_Y];
				old_dp[A_Z]=dp[A_Z];
				old_dq[A_X]=dq[A_X];
				old_dq[A_Y]=dq[A_Y];
				old_dq[A_Z]=dq[A_Z];

				for(int i=0;i<3*Nx;i++)
				{
					d[i]=0;
					Nr[i]=0;
					for(int j=0;j<3*Nx;j++)
					{
						if(j==i)	B[i*3*Nx+j]=1;
						else
						{
							B[i*3*Nx+j]=0;
						}
					}
				}

				Tr=0.5/mi*( 2*Dt*(pn[A_X]*dp[A_X]+pn[A_Y]*dp[A_Y]+pn[A_Z]*dp[A_Z]) + Dt*Dt*(dp[A_X]*dp[A_X]+dp[A_Y]*dp[A_Y]+dp[A_Z]*dp[A_Z]) ) -mi*G*Dt*(dq[A_X]*nG[A_X]+dq[A_Y]*nG[A_Y]+dq[A_Z]*nG[A_Z]);

				h=-1*( (qn[A_X]+Dt*dq[A_X]-a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq[A_Z]-a[A_Z])*n[A_Z]);
			
				if(h+theta>0)	Tr+=0.5*r*(h+theta)*(h+theta);

				dTr[0]=Dt/mi*(pn[A_X]+Dt*dp[A_X]);
				dTr[1]=Dt/mi*(pn[A_Y]+Dt*dp[A_Y]);
				dTr[2]=Dt/mi*(pn[A_Z]+Dt*dp[A_Z]);
				dTr[3]=-mi*G*Dt*nG[A_X];
				dTr[4]=-mi*G*Dt*nG[A_Y];
				dTr[5]=-mi*G*Dt*nG[A_Z];
				if(h+theta>0)
				{
					dTr[3]+=-r*Dt*n[A_X]*(h+theta);
					dTr[4]+=-r*Dt*n[A_Y]*(h+theta);
					dTr[5]+=-r*Dt*n[A_Z]*(h+theta);
				}
			
				double E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);
				cout<<"E0="<<E<<endl;

				int count=0;
 				if(E<ep)	break;
				while(E>ep)
				{
					count++;
					if(count==1)
					{
						stringstream se;
						se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_e(se.str(), ios::trunc);
						init_e.close();

						stringstream ssp;
						ssp<<"./dp/dp_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_p(ssp.str(), ios::trunc);
						init_p.close();
	
						stringstream ssq;
						ssq<<"./dq/dq_t"<<t<<"cmin"<<count_min<<".csv";
						ofstream init_q(ssq.str(), ios::trunc);
						init_q.close();
					}
					if(count>500)	break;
					double dp_k[DIMENSION]={dp[A_X],dp[A_Y], dp[A_Z]};
					double dq_k[DIMENSION]={dq[A_X],dq[A_Y], dq[A_Z]};
					double dTr_k[6]={dTr[0], dTr[1], dTr[2], dTr[3], dTr[4], dTr[5]};

					for(int i=0;i<3*Nx;i++)	Nr[i]=dTr[i];
					if(count>1)
					{
/*						cout<<"Nr="<<Nr[0]<<", "<<Nr[1]<<", "<<Nr[2]<<", "<<Nr[3]<<", "<<Nr[4]<<", "<<Nr[5]<<endl;
						cout<<"B=";
						for(int i=0;i<3*Nx;i++)
						{
							for(int j=0;j<3*Nx;j++)
							{
								cout<<B[i*3*Nx+j]<<", ";
							}
							cout<<endl;
						}*/

						gauss(B,Nr,3*Nx);
					}

					for(int i=0;i<3*Nx;i++)	d[i]=-1*Nr[i];

	//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;
//					cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<", "<<d[4]<<", "<<d[5]<<endl;

					double Tr_min=Tr;
					double a_min=0.1;//1e-4;

					for(int i=0;i<10000;i++)
					{
						double alpha=(i+1)*1e-4;
						double dp_a[DIMENSION]={dp[A_X]+d[0]*alpha, dp[A_Y]+d[1]*alpha, dp[A_Z]+d[2]*alpha};
						double dq_a[DIMENSION]={dq[A_X]+d[3]*alpha, dq[A_Y]+d[4]*alpha, dq[A_Z]+d[5]*alpha};

						double Tr_a=0.5/mi*( 2*Dt*(pn[A_X]*dp_a[A_X]+pn[A_Y]*dp_a[A_Y]+pn[A_Z]*dp_a[A_Z]) + Dt*Dt*(dp_a[A_X]*dp_a[A_X]+dp_a[A_Y]*dp_a[A_Y]+dp_a[A_Z]*dp_a[A_Z]) )
							-mi*G*Dt*( dq_a[A_X]*nG[A_X] + dq_a[A_Y]*nG[A_Y] + dq_a[A_Z]*nG[A_Z] );
						double h_a=-1*( (qn[A_X]+Dt*dq_a[A_X]-a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq_a[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq_a[A_Z]-a[A_Z])*n[A_Z]);
						if(h_a+theta>0)	Tr+=0.5*r*(h_a+theta)*(h_a+theta);

						if(Tr_a<Tr_min)
						{
							Tr_min=Tr_a;
							a_min=alpha;
						}
					}//*/
//					cout<<"Tr"<<count<<"="<<Tr_min<<", alpha="<<a_min<<endl;
				
					dp[A_X]+=d[0]*a_min;
					dp[A_Y]+=d[1]*a_min;
					dp[A_Z]+=d[2]*a_min;
					dq[A_X]+=d[3]*a_min;
					dq[A_Y]+=d[4]*a_min;
					dq[A_Z]+=d[5]*a_min;
//					cout<<"dp"<<count<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//					cout<<"dq"<<count<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
					Tr=0.5/mi*( 2*Dt*( pn[A_X]*dp[A_X] + pn[A_Y]*dp[A_Y] + pn[A_Z]*dp[A_Z] ) + Dt*Dt*( dp[A_X]*dp[A_X] + dp[A_Y]*dp[A_Y] + dp[A_Z]*dp[A_Z] ) )
						-mi*G*Dt*( dq[A_X]*nG[A_X] + dq[A_Y]*nG[A_Y] + dq[A_Z]*nG[A_Z] );
					h=-1*( (qn[A_X]+Dt*dq[A_X]-a[A_X])*n[A_X] + (qn[A_Y]+Dt*dq[A_Y]-a[A_Y])*n[A_Y] + (qn[A_Z]+Dt*dq[A_Z]-a[A_Z])*n[A_Z]);
					if(h+theta>0)	Tr+=0.5*r*(h+theta)*(h+theta);
//					cout<<"h"<<count<<"="<<h<<endl;

					dTr[0]=Dt/mi*(pn[A_X]+Dt*dp[A_X]);
					dTr[1]=Dt/mi*(pn[A_Y]+Dt*dp[A_Y]);
					dTr[2]=Dt/mi*(pn[A_Z]+Dt*dp[A_Z]);
					dTr[3]=-mi*G*Dt*nG[A_X];
					dTr[4]=-mi*G*Dt*nG[A_Y];
					dTr[5]=-mi*G*Dt*nG[A_Z];
					if(h+theta>0)
					{
						dTr[3]+=-r*Dt*n[A_X]*(h+theta);
						dTr[4]+=-r*Dt*n[A_Y]*(h+theta);
						dTr[5]+=-r*Dt*n[A_Z]*(h+theta);
					}
			
					double old_E=E;

					E=sqrt(dTr[0]*dTr[0]+dTr[1]*dTr[1]+dTr[2]*dTr[2]+dTr[3]*dTr[3]+dTr[4]*dTr[4]+dTr[5]*dTr[5]);

					cout<<"E"<<count<<"="<<E<<endl;

					stringstream se;
					se<<"./E/E_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fe(se.str(), ios::app);
					fe<<count<<","<<E<<endl;
					fe.close();

					stringstream ssp;
					ssp<<"./dp/dp_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fp(ssp.str(), ios::app);
					fp<<count<<","<<dp[A_X]<<","<<dp[A_Y]<<","<<dp[A_Z]<<endl;
					fp.close();
	
					stringstream ssq;
					ssq<<"./dq/dq_t"<<t<<"cmin"<<count_min<<".csv";
					ofstream fq(ssq.str(), ios::app);
					fq<<count<<","<<dq[A_X]<<","<<dq[A_Y]<<","<<dq[A_Z]<<endl;
					fq.close();

					if(E<ep)	break;
					/*
					if(E>old_E && h<=0)
					{
						dp[A_X]=dp_k[A_X];	dp[A_Y]=dp_k[A_Y];	dp[A_X]=dp_k[A_Z];
						dq[A_X]=dq_k[A_X];	dq[A_Y]=dq_k[A_Y];	dq[A_X]=dq_k[A_Z];
						break;
					}*/
				

					double old_E_min=E_min;

					double s[6]={dp[A_X]-dp_k[A_X],dp[A_Y]-dp_k[A_Y], dp[A_Z]-dp_k[A_Z], dq[A_X]-dq_k[A_X], dq[A_Y]-dq_k[A_Y], dq[A_Z]-dq_k[A_Z]};
					double y[6]={dTr[0]-dTr_k[0],dTr[1]-dTr_k[1], dTr[2]-dTr_k[2], dTr[3]-dTr_k[3], dTr[4]-dTr_k[4], dTr[5]-dTr_k[5]};

					double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3]+y[4]*s[4]+y[5]*s[5];
					double sigma=(s[0]*B[0*Nx*3+0]+s[1]*B[1*Nx*3+0]+s[2]*B[2*Nx*3+0]+s[3]*B[3*Nx*3+0]+s[4]*B[4*Nx*3+0]+s[5]*B[5*Nx*3+0])*s[0]
								+(s[0]*B[0*Nx*3+1]+s[1]*B[1*Nx*3+1]+s[2]*B[2*Nx*3+1]+s[3]*B[3*Nx*3+1]+s[4]*B[4*Nx*3+1]+s[5]*B[5*Nx*3+1])*s[1]
								+(s[0]*B[0*Nx*3+2]+s[1]*B[1*Nx*3+2]+s[2]*B[2*Nx*3+2]+s[3]*B[3*Nx*3+2]+s[4]*B[4*Nx*3+2]+s[5]*B[5*Nx*3+2])*s[2]
								+(s[0]*B[0*Nx*3+3]+s[1]*B[1*Nx*3+3]+s[2]*B[2*Nx*3+3]+s[3]*B[3*Nx*3+3]+s[4]*B[4*Nx*3+3]+s[5]*B[5*Nx*3+3])*s[3]
								+(s[0]*B[0*Nx*3+4]+s[1]*B[1*Nx*3+4]+s[2]*B[2*Nx*3+4]+s[3]*B[3*Nx*3+4]+s[4]*B[4*Nx*3+4]+s[5]*B[5*Nx*3+4])*s[4]
								+(s[0]*B[0*Nx*3+5]+s[1]*B[1*Nx*3+5]+s[2]*B[2*Nx*3+5]+s[3]*B[3*Nx*3+5]+s[4]*B[4*Nx*3+5]+s[5]*B[5*Nx*3+5])*s[5];

					cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

					if(beta>0)//(beta>=0.2*sigma)
					{
						double bs[6]={  B[0*3*Nx+0]*s[0]+B[0*3*Nx+1]*s[1]+B[0*3*Nx+2]*s[2]+B[0*3*Nx+3]*s[3]+B[0*3*Nx+4]*s[4]+B[0*3*Nx+5]*s[5],
										B[1*3*Nx+0]*s[0]+B[1*3*Nx+1]*s[1]+B[1*3*Nx+2]*s[2]+B[1*3*Nx+3]*s[3]+B[1*3*Nx+4]*s[4]+B[1*3*Nx+5]*s[5],
										B[2*3*Nx+0]*s[0]+B[2*3*Nx+1]*s[1]+B[2*3*Nx+2]*s[2]+B[2*3*Nx+3]*s[3]+B[2*3*Nx+4]*s[4]+B[2*3*Nx+5]*s[5],
										B[3*3*Nx+0]*s[0]+B[3*3*Nx+1]*s[1]+B[3*3*Nx+2]*s[2]+B[3*3*Nx+3]*s[3]+B[3*3*Nx+4]*s[4]+B[3*3*Nx+5]*s[5],
										B[4*3*Nx+0]*s[0]+B[4*3*Nx+1]*s[1]+B[4*3*Nx+2]*s[2]+B[4*3*Nx+3]*s[3]+B[4*3*Nx+4]*s[4]+B[4*3*Nx+5]*s[5],
										B[5*3*Nx+0]*s[0]+B[5*3*Nx+1]*s[1]+B[5*3*Nx+2]*s[2]+B[5*3*Nx+3]*s[3]+B[5*3*Nx+4]*s[4]+B[5*3*Nx+5]*s[5]};

						double bss[36]={bs[0]*s[0],bs[0]*s[1],bs[0]*s[2],bs[0]*s[3],bs[0]*s[4],bs[0]*s[5],
										bs[1]*s[0],bs[1]*s[1],bs[1]*s[2],bs[1]*s[3],bs[1]*s[4],bs[1]*s[5],
										bs[2]*s[0],bs[2]*s[1],bs[2]*s[2],bs[2]*s[3],bs[2]*s[4],bs[2]*s[5],
										bs[3]*s[0],bs[3]*s[1],bs[3]*s[2],bs[3]*s[3],bs[3]*s[4],bs[3]*s[5],
										bs[4]*s[0],bs[4]*s[1],bs[4]*s[2],bs[4]*s[3],bs[4]*s[4],bs[4]*s[5],
										bs[5]*s[0],bs[5]*s[1],bs[5]*s[2],bs[5]*s[3],bs[5]*s[4],bs[5]*s[5]};

						B[0*3*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+0]+bss[0*3*Nx+1]*B[1*3*Nx+0]+bss[0*3*Nx+2]*B[2*3*Nx+0]+bss[0*3*Nx+3]*B[3*3*Nx+0]+bss[0*3*Nx+4]*B[4*3*Nx+0]+bss[0*3*Nx+5]*B[5*3*Nx+0]);
						B[0*3*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+1]+bss[0*3*Nx+1]*B[1*3*Nx+1]+bss[0*3*Nx+2]*B[2*3*Nx+1]+bss[0*3*Nx+3]*B[3*3*Nx+1]+bss[0*3*Nx+4]*B[4*3*Nx+1]+bss[0*3*Nx+5]*B[5*3*Nx+1]);
						B[0*3*Nx+2]+=1/beta*y[0]*y[2]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+2]+bss[0*3*Nx+1]*B[1*3*Nx+2]+bss[0*3*Nx+2]*B[2*3*Nx+2]+bss[0*3*Nx+3]*B[3*3*Nx+2]+bss[0*3*Nx+4]*B[4*3*Nx+2]+bss[0*3*Nx+5]*B[5*3*Nx+2]);
						B[0*3*Nx+3]+=1/beta*y[0]*y[3]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+3]+bss[0*3*Nx+1]*B[1*3*Nx+3]+bss[0*3*Nx+2]*B[2*3*Nx+3]+bss[0*3*Nx+3]*B[3*3*Nx+3]+bss[0*3*Nx+4]*B[4*3*Nx+3]+bss[0*3*Nx+5]*B[5*3*Nx+3]);
						B[0*3*Nx+4]+=1/beta*y[0]*y[4]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+4]+bss[0*3*Nx+1]*B[1*3*Nx+4]+bss[0*3*Nx+2]*B[2*3*Nx+4]+bss[0*3*Nx+3]*B[3*3*Nx+4]+bss[0*3*Nx+4]*B[4*3*Nx+4]+bss[0*3*Nx+5]*B[5*3*Nx+4]);
						B[0*3*Nx+5]+=1/beta*y[0]*y[5]-1/sigma*(bss[0*3*Nx+0]*B[0*3*Nx+5]+bss[0*3*Nx+1]*B[1*3*Nx+5]+bss[0*3*Nx+2]*B[2*3*Nx+5]+bss[0*3*Nx+3]*B[3*3*Nx+5]+bss[0*3*Nx+4]*B[4*3*Nx+5]+bss[0*3*Nx+5]*B[5*3*Nx+5]);
				
						B[1*3*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+0]+bss[1*3*Nx+1]*B[1*3*Nx+0]+bss[1*3*Nx+2]*B[2*3*Nx+0]+bss[1*3*Nx+3]*B[3*3*Nx+0]+bss[1*3*Nx+4]*B[4*3*Nx+0]+bss[1*3*Nx+5]*B[5*3*Nx+0]);
						B[1*3*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+1]+bss[1*3*Nx+1]*B[1*3*Nx+1]+bss[1*3*Nx+2]*B[2*3*Nx+1]+bss[1*3*Nx+3]*B[3*3*Nx+1]+bss[1*3*Nx+4]*B[4*3*Nx+1]+bss[1*3*Nx+5]*B[5*3*Nx+1]);
						B[1*3*Nx+2]+=1/beta*y[1]*y[2]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+2]+bss[1*3*Nx+1]*B[1*3*Nx+2]+bss[1*3*Nx+2]*B[2*3*Nx+2]+bss[1*3*Nx+3]*B[3*3*Nx+2]+bss[1*3*Nx+4]*B[4*3*Nx+2]+bss[1*3*Nx+5]*B[5*3*Nx+2]);
						B[1*3*Nx+3]+=1/beta*y[1]*y[3]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+3]+bss[1*3*Nx+1]*B[1*3*Nx+3]+bss[1*3*Nx+2]*B[2*3*Nx+3]+bss[1*3*Nx+3]*B[3*3*Nx+3]+bss[1*3*Nx+4]*B[4*3*Nx+3]+bss[1*3*Nx+5]*B[5*3*Nx+3]);
						B[1*3*Nx+4]+=1/beta*y[1]*y[4]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+4]+bss[1*3*Nx+1]*B[1*3*Nx+4]+bss[1*3*Nx+2]*B[2*3*Nx+4]+bss[1*3*Nx+3]*B[3*3*Nx+4]+bss[1*3*Nx+4]*B[4*3*Nx+4]+bss[1*3*Nx+5]*B[5*3*Nx+4]);
						B[1*3*Nx+5]+=1/beta*y[1]*y[5]-1/sigma*(bss[1*3*Nx+0]*B[0*3*Nx+5]+bss[1*3*Nx+1]*B[1*3*Nx+5]+bss[1*3*Nx+2]*B[2*3*Nx+5]+bss[1*3*Nx+3]*B[3*3*Nx+5]+bss[1*3*Nx+4]*B[4*3*Nx+5]+bss[1*3*Nx+5]*B[5*3*Nx+5]);

						B[2*3*Nx+0]+=1/beta*y[2]*y[0]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+0]+bss[2*3*Nx+1]*B[1*3*Nx+0]+bss[2*3*Nx+2]*B[2*3*Nx+0]+bss[2*3*Nx+3]*B[3*3*Nx+0]+bss[2*3*Nx+4]*B[4*3*Nx+0]+bss[2*3*Nx+5]*B[5*3*Nx+0]);
						B[2*3*Nx+1]+=1/beta*y[2]*y[1]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+1]+bss[2*3*Nx+1]*B[1*3*Nx+1]+bss[2*3*Nx+2]*B[2*3*Nx+1]+bss[2*3*Nx+3]*B[3*3*Nx+1]+bss[2*3*Nx+4]*B[4*3*Nx+1]+bss[2*3*Nx+5]*B[5*3*Nx+1]);
						B[2*3*Nx+2]+=1/beta*y[2]*y[2]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+2]+bss[2*3*Nx+1]*B[1*3*Nx+2]+bss[2*3*Nx+2]*B[2*3*Nx+2]+bss[2*3*Nx+3]*B[3*3*Nx+2]+bss[2*3*Nx+4]*B[4*3*Nx+2]+bss[2*3*Nx+5]*B[5*3*Nx+2]);
						B[2*3*Nx+3]+=1/beta*y[2]*y[3]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+3]+bss[2*3*Nx+1]*B[1*3*Nx+3]+bss[2*3*Nx+2]*B[2*3*Nx+3]+bss[2*3*Nx+3]*B[3*3*Nx+3]+bss[2*3*Nx+4]*B[4*3*Nx+3]+bss[2*3*Nx+5]*B[5*3*Nx+3]);
						B[2*3*Nx+4]+=1/beta*y[2]*y[4]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+4]+bss[2*3*Nx+1]*B[1*3*Nx+4]+bss[2*3*Nx+2]*B[2*3*Nx+4]+bss[2*3*Nx+3]*B[3*3*Nx+4]+bss[2*3*Nx+4]*B[4*3*Nx+4]+bss[2*3*Nx+5]*B[5*3*Nx+4]);
						B[2*3*Nx+5]+=1/beta*y[2]*y[5]-1/sigma*(bss[2*3*Nx+0]*B[0*3*Nx+5]+bss[2*3*Nx+1]*B[1*3*Nx+5]+bss[2*3*Nx+2]*B[2*3*Nx+5]+bss[2*3*Nx+3]*B[3*3*Nx+5]+bss[2*3*Nx+4]*B[4*3*Nx+5]+bss[2*3*Nx+5]*B[5*3*Nx+5]);

						B[3*3*Nx+0]+=1/beta*y[3]*y[0]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+0]+bss[3*3*Nx+1]*B[1*3*Nx+0]+bss[3*3*Nx+2]*B[2*3*Nx+0]+bss[3*3*Nx+3]*B[3*3*Nx+0]+bss[3*3*Nx+4]*B[4*3*Nx+0]+bss[3*3*Nx+5]*B[5*3*Nx+0]);
						B[3*3*Nx+1]+=1/beta*y[3]*y[1]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+1]+bss[3*3*Nx+1]*B[1*3*Nx+1]+bss[3*3*Nx+2]*B[2*3*Nx+1]+bss[3*3*Nx+3]*B[3*3*Nx+1]+bss[3*3*Nx+4]*B[4*3*Nx+1]+bss[3*3*Nx+5]*B[5*3*Nx+1]);
						B[3*3*Nx+2]+=1/beta*y[3]*y[2]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+2]+bss[3*3*Nx+1]*B[1*3*Nx+2]+bss[3*3*Nx+2]*B[2*3*Nx+2]+bss[3*3*Nx+3]*B[3*3*Nx+2]+bss[3*3*Nx+4]*B[4*3*Nx+2]+bss[3*3*Nx+5]*B[5*3*Nx+2]);
						B[3*3*Nx+3]+=1/beta*y[3]*y[3]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+3]+bss[3*3*Nx+1]*B[1*3*Nx+3]+bss[3*3*Nx+2]*B[2*3*Nx+3]+bss[3*3*Nx+3]*B[3*3*Nx+3]+bss[3*3*Nx+4]*B[4*3*Nx+3]+bss[3*3*Nx+5]*B[5*3*Nx+3]);
						B[3*3*Nx+4]+=1/beta*y[3]*y[4]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+4]+bss[3*3*Nx+1]*B[1*3*Nx+4]+bss[3*3*Nx+2]*B[2*3*Nx+4]+bss[3*3*Nx+3]*B[3*3*Nx+4]+bss[3*3*Nx+4]*B[4*3*Nx+4]+bss[3*3*Nx+5]*B[5*3*Nx+4]);
						B[3*3*Nx+5]+=1/beta*y[3]*y[5]-1/sigma*(bss[3*3*Nx+0]*B[0*3*Nx+5]+bss[3*3*Nx+1]*B[1*3*Nx+5]+bss[3*3*Nx+2]*B[2*3*Nx+5]+bss[3*3*Nx+3]*B[3*3*Nx+5]+bss[3*3*Nx+4]*B[4*3*Nx+5]+bss[3*3*Nx+5]*B[5*3*Nx+5]);

						B[4*3*Nx+0]+=1/beta*y[4]*y[0]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+0]+bss[4*3*Nx+1]*B[1*3*Nx+0]+bss[4*3*Nx+2]*B[2*3*Nx+0]+bss[4*3*Nx+3]*B[3*3*Nx+0]+bss[4*3*Nx+4]*B[4*3*Nx+0]+bss[4*3*Nx+5]*B[5*3*Nx+0]);
						B[4*3*Nx+1]+=1/beta*y[4]*y[1]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+1]+bss[4*3*Nx+1]*B[1*3*Nx+1]+bss[4*3*Nx+2]*B[2*3*Nx+1]+bss[4*3*Nx+3]*B[3*3*Nx+1]+bss[4*3*Nx+4]*B[4*3*Nx+1]+bss[4*3*Nx+5]*B[5*3*Nx+1]);
						B[4*3*Nx+2]+=1/beta*y[4]*y[2]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+2]+bss[4*3*Nx+1]*B[1*3*Nx+2]+bss[4*3*Nx+2]*B[2*3*Nx+2]+bss[4*3*Nx+3]*B[3*3*Nx+2]+bss[4*3*Nx+4]*B[4*3*Nx+2]+bss[4*3*Nx+5]*B[5*3*Nx+2]);
						B[4*3*Nx+3]+=1/beta*y[4]*y[3]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+3]+bss[4*3*Nx+1]*B[1*3*Nx+3]+bss[4*3*Nx+2]*B[2*3*Nx+3]+bss[4*3*Nx+3]*B[3*3*Nx+3]+bss[4*3*Nx+4]*B[4*3*Nx+3]+bss[4*3*Nx+5]*B[5*3*Nx+3]);
						B[4*3*Nx+4]+=1/beta*y[4]*y[4]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+4]+bss[4*3*Nx+1]*B[1*3*Nx+4]+bss[4*3*Nx+2]*B[2*3*Nx+4]+bss[4*3*Nx+3]*B[3*3*Nx+4]+bss[4*3*Nx+4]*B[4*3*Nx+4]+bss[4*3*Nx+5]*B[5*3*Nx+4]);
						B[4*3*Nx+5]+=1/beta*y[4]*y[5]-1/sigma*(bss[4*3*Nx+0]*B[0*3*Nx+5]+bss[4*3*Nx+1]*B[1*3*Nx+5]+bss[4*3*Nx+2]*B[2*3*Nx+5]+bss[4*3*Nx+3]*B[3*3*Nx+5]+bss[4*3*Nx+4]*B[4*3*Nx+5]+bss[4*3*Nx+5]*B[5*3*Nx+5]);

						B[5*3*Nx+0]+=1/beta*y[5]*y[0]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+0]+bss[5*3*Nx+1]*B[1*3*Nx+0]+bss[5*3*Nx+2]*B[2*3*Nx+0]+bss[5*3*Nx+3]*B[3*3*Nx+0]+bss[5*3*Nx+4]*B[4*3*Nx+0]+bss[5*3*Nx+5]*B[5*3*Nx+0]);
						B[5*3*Nx+1]+=1/beta*y[5]*y[1]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+1]+bss[5*3*Nx+1]*B[1*3*Nx+1]+bss[5*3*Nx+2]*B[2*3*Nx+1]+bss[5*3*Nx+3]*B[3*3*Nx+1]+bss[5*3*Nx+4]*B[4*3*Nx+1]+bss[5*3*Nx+5]*B[5*3*Nx+1]);
						B[5*3*Nx+2]+=1/beta*y[5]*y[2]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+2]+bss[5*3*Nx+1]*B[1*3*Nx+2]+bss[5*3*Nx+2]*B[2*3*Nx+2]+bss[5*3*Nx+3]*B[3*3*Nx+2]+bss[5*3*Nx+4]*B[4*3*Nx+2]+bss[5*3*Nx+5]*B[5*3*Nx+2]);
						B[5*3*Nx+3]+=1/beta*y[5]*y[3]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+3]+bss[5*3*Nx+1]*B[1*3*Nx+3]+bss[5*3*Nx+2]*B[2*3*Nx+3]+bss[5*3*Nx+3]*B[3*3*Nx+3]+bss[5*3*Nx+4]*B[4*3*Nx+3]+bss[5*3*Nx+5]*B[5*3*Nx+3]);
						B[5*3*Nx+4]+=1/beta*y[5]*y[4]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+4]+bss[5*3*Nx+1]*B[1*3*Nx+4]+bss[5*3*Nx+2]*B[2*3*Nx+4]+bss[5*3*Nx+3]*B[3*3*Nx+4]+bss[5*3*Nx+4]*B[4*3*Nx+4]+bss[5*3*Nx+5]*B[5*3*Nx+4]);
						B[5*3*Nx+5]+=1/beta*y[5]*y[5]-1/sigma*(bss[5*3*Nx+0]*B[0*3*Nx+5]+bss[5*3*Nx+1]*B[1*3*Nx+5]+bss[5*3*Nx+2]*B[2*3*Nx+5]+bss[5*3*Nx+3]*B[3*3*Nx+5]+bss[5*3*Nx+4]*B[4*3*Nx+5]+bss[5*3*Nx+5]*B[5*3*Nx+5]);
					}
				}
				double old_E_min=E_min;
				E_min=sqrt( (old_dp[A_X]-dp[A_X])*(old_dp[A_X]-dp[A_X]) + (old_dp[A_Y]-dp[A_Y])*(old_dp[A_Y]-dp[A_Y]) + (old_dp[A_Z]-dp[A_Z])*(old_dp[A_Z]-dp[A_Z]) +
							(old_dq[A_X]-dq[A_X])*(old_dq[A_X]-dq[A_X]) + (old_dq[A_Y]-dq[A_Y])*(old_dq[A_Y]-dq[A_Y]) + (old_dq[A_Z]-dq[A_Z])*(old_dq[A_Z]-dq[A_Z]));	
				if(E_min<ep_min*1000)	r*=4;
				if(h+theta>0)	theta+=h;
	
//				cout<<"E_min"<<count_min<<"="<<E_min<<endl;
				
				stringstream se_min;
				se_min<<"./E_min"<<t<<".csv";
				ofstream f_emin(se_min.str(), ios::app);
				f_emin<<E_min<<endl;
				f_emin.close();

				stringstream sh;
				sh<<"./h"<<t<<".csv";
				ofstream fh(sh.str(), ios::app);
				fh<<h<<endl;
				fh.close();

				stringstream se_c;
				se_c<<"./count"<<t<<".csv";
				ofstream fc(se_c.str(), ios::app);
				fc<<count_min<<","<<count<<endl;
				fc.close();

		
				stringstream se_p;
				se_p<<"./dp"<<t<<".csv";
				ofstream fp(se_p.str(), ios::app);
				fp<<count_min<<","<<dp[A_X]<<","<<dp[A_Y]<<","<<dp[A_Z]<<endl;
				fp.close();

				stringstream se_q;
				se_q<<"./dq"<<t<<".csv";
				ofstream fq(se_q.str(), ios::app);
				fq<<count_min<<","<<dq[A_X]<<","<<dq[A_Y]<<","<<dq[A_Z]<<endl;
				fq.close();

/*				if(E_min>old_E_min)
				{
					dp[A_X]=old_dp[A_X];	dp[A_Y]=old_dp[A_Y];	dp[A_X]=old_dp[A_Z];
					dq[A_X]=old_dq[A_X];	dq[A_Y]=old_dq[A_Y];	dq[A_X]=old_dq[A_Z];
					break;
				}*/
			}
//			cout<<"dp"<<count_min<<"="<<dp[0]<<", "<<dp[1]<<", "<<dp[2]<<endl;
//			cout<<"dq"<<count_min<<"="<<dq[0]<<", "<<dq[1]<<", "<<dq[2]<<endl;
		
			qn[A_X]+=Dt*dq[A_X];	qn[A_Y]+=Dt*dq[A_Y];	qn[A_Z]+=Dt*dq[A_Z];
			pn[A_X]+=Dt*dp[A_X];	pn[A_Y]+=Dt*dp[A_Y];	pn[A_Z]+=Dt*dp[A_Z];

			cout<<"pn"<<t<<"="<<pn[0]<<", "<<pn[1]<<", "<<pn[2]<<endl;
			cout<<"qn"<<t<<"="<<qn[0]<<", "<<qn[1]<<", "<<qn[2]<<endl;
			cout<<"------------------------------OK"<<endl;

			ofstream fcmin("cmin.csv", ios::app);
			fcmin<<t<<","<<count_min<<endl;
			fcmin.close();
	
		}
		ofstream fp("p.csv", ios::app);
		fp<<pn[A_X]<<","<<pn[A_Y]<<","<<pn[A_Z]<<endl;
		fp.close();
		ofstream fq("q.csv", ios::app);
		fq<<qn[A_X]<<","<<qn[A_Y]<<","<<qn[A_Z]<<endl;
		fq.close();
		ofstream fh("h.csv", ios::app);
		fh<<h<<endl;
		fh.close();
		ofstream fEn("Energy.csv", ios::app);
		double En=0.5/mi*(pn[A_X]*pn[A_X]+pn[A_Y]*pn[A_Y]+pn[A_Z]*pn[A_Z])-mi*G*Dt*(qn[A_X]*nG[A_X]+qn[A_Y]*nG[A_Y]+qn[A_Z]*nG[A_Z]);
		double error=(En-E0)/E0;
		fEn<<En<<","<<error<<endl;
		fEn.close();


		double *q=new double [DIMENSION];
		double *p=new double [DIMENSION];

		p[A_X]=pn[A_X];	p[A_Y]=pn[A_Y];	p[A_Z]=pn[A_Z];
		q[A_X]=qn[A_X];	q[A_Y]=qn[A_Y];	q[A_Z]=qn[A_Z];

		PM_AVS(t_max,t,q);
		MM_AVS(t_max,t,Dt,mi,q, p);
		
		delete[]	q;
		delete[]	p;

		clock_t t1=clock();
		ofstream t_log("time_log.dat", ios::app);
		t_log<<"step="<<t<<", time="<<(t1-t0)/CLOCKS_PER_SEC<<"[sec]"<<endl;
		t_log.close();


	}
	
	//v0　計算
//	int nv=0;
//	cout<<"fc="<<c0<<", "<<c1<<endl;
	

	//double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);
	//double v1=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	//cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	//cout<<"c="<<c0<<endl;
	//cout<<"v="<<v0<<endl;
	//cout<<"fx="<<fx<<endl;
	//cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/

	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}


int QP_with_2_newton_way()	
{
	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double *x=new double [Nx];
	double *old_x=new double [Nx];
	double x0=-20,x1=-19;
	
	x[0]=x0;
	x[1]=x1;


	double c0=0;
	double fx=0;
	double Txr=0;

	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dTxr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rTxr[i*Nx+j]=0;
		}
	}


	int calc_way=1;	//最適解探索	Newton 0, 準Newton 1, 最急降下法2

	double r=1;
	double ep=1e-10;
	double theta0=0, theta1=0, theta2=0, theta3=0;
	double E_min=1;
	int count_min=0;

	double *B=new double[Nx*Nx];
	double *d=new double [Nx];	
	double *Nr=new double [Nx];


		
	ofstream init_emin("E.csv", ios::trunc);
	init_emin.close();					
		
	ofstream init_x("x.csv", ios::trunc);
	init_x.close();
	ofstream init_fx("fx.csv", ios::trunc);
	init_fx.close();
	ofstream init_t("time_log.csv", ios::trunc);
	init_t.close();

	while(E_min>ep*100000)
	{
		count_min++;
		clock_t t0=clock();


		if(count_min>5000)	break;
		for(int i=0;i<Nx;i++)
		{
			old_x[i]=x[i];
			d[i]=0;
			Nr[i]=0;
			for(int j=0;j<Nx;j++)
			{
				if(j==i)	B[i*Nx+j]=1;
				else
				{
					B[i*Nx+j]=0;
				}
			}
		}


		double E=1;
		int count=0;

		
		stringstream se_min;
		se_min<<"./E"<<count_min<<".csv";
		ofstream init_emin(se_min.str(), ios::trunc);
		init_emin.close();					
		
		stringstream se_x;
		se_x<<"./x"<<count_min<<".csv";
		ofstream init_x(se_x.str(), ios::trunc);
		init_x.close();
	
		stringstream se_fx;
		se_fx<<"./fx"<<count_min<<".csv";
		ofstream init_fx(se_fx.str(), ios::trunc);
		init_fx.close();

		if(calc_way==0)
		{
			while(E>ep)
			{
				count++;

				stringstream se_min;
				se_min<<"./E"<<count_min<<".csv";
				ofstream init_emin(se_min.str(), ios::app);
				init_emin<<E<<endl;
				init_emin.close();					
		
				stringstream se_x;
				se_x<<"./x"<<count_min<<".csv";
				ofstream init_x(se_x.str(), ios::app);
				init_x<<x[0]<<","<<x[1]<<endl;
				init_x.close();
	
				stringstream se_fx;
				se_fx<<"./fx"<<count_min<<".csv";
				ofstream init_fx(se_fx.str(), ios::app);
				init_fx<<fx<<endl;
				init_fx.close();

				fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
				dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);
				rfx[0*Nx+0]=2;	rfx[0*Nx+1]=0;	rfx[1*Nx+0]=0;	rfx[1*Nx+1]=2;

				c0=x[0]*x[0]+x[1]*x[1]-1;
				dc0[0]=2*x[0];	dc0[1]=2*x[1];
				rc0[0*Nx+0]=2;	rc0[0*Nx+1]=0;	rc0[1*Nx+0]=0;	rc0[1*Nx+1]=2;

				Txr=fx;
				if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);
	
				//cout<<"Txr="<<Txr<<endl;

				for(int i=0;i<Nx;i++)
				{
					dTxr[i]=dfx[i];
					if(c0+theta0>0)	dTxr[i]+=r*(c0+theta0)*dc0[i];
				}

				//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;

				for(int i=0;i<Nx;i++)
				{
					for(int j=0;j<Nx;j++)
					{
						rTxr[i*Nx+j]=rfx[i*Nx+j];
						if(c0+theta0>0)	rTxr[i*Nx+j]+=r*dc0[j]*dc0[i]+r*(c0+theta0)*rc0[i*Nx+j];
					}
				}
				//cout<<"rTxr="<<rTxr[0]<<", "<<rTxr[1]<<", "<<rTxr[2]<<", "<<rTxr[3]<<endl;
				//cout<<rTxr[4]<<", "<<rTxr[5]<<", "<<rTxr[6]<<", "<<rTxr[7]<<endl;
				//cout<<rTxr[8]<<", "<<rTxr[9]<<", "<<rTxr[10]<<", "<<rTxr[11]<<endl;
				//cout<<rTxr[12]<<", "<<rTxr[13]<<", "<<rTxr[14]<<", "<<rTxr[15]<<endl;

				gauss(rTxr,dTxr,Nx);
				//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;
				for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];


				E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
//				cout<<"E"<<count<<"="<<E<<endl;
			}
		}
		else if(calc_way==1)
		{
	
			fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
			dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

			c0=x[0]*x[0]+x[1]*x[1]-1;
			dc0[0]=2*x[0];	dc0[1]=2*x[1];
			
			Txr=fx;
			if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);			
			dTxr[0]=dfx[0];
			dTxr[1]=dfx[1];
			if(c0+theta0>0)
			{
				dTxr[0]+=r*(c0+theta0)*dc0[0];
				dTxr[1]+=r*(c0+theta0)*dc0[1];
			}
			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
//			cout<<"E0="<<E<<endl;
			if(E<ep)	return 0;
			while(E>ep)
			{
				count++;

				stringstream se_min;
				se_min<<"./E"<<count_min<<".csv";
				ofstream init_emin(se_min.str(), ios::app);
				init_emin<<E<<endl;
				init_emin.close();					
		
				stringstream se_x;
				se_x<<"./x"<<count_min<<".csv";
				ofstream init_x(se_x.str(), ios::app);
				init_x<<x[0]<<","<<x[1]<<endl;
				init_x.close();
	
				stringstream se_fx;
				se_fx<<"./fx"<<count_min<<".csv";
				ofstream init_fx(se_fx.str(), ios::app);
				init_fx<<fx<<endl;
				init_fx.close();

				if(count>500)	break;
				double x_k[2]={x[0], x[1]};
				double dTxr_k[2]={dTxr[0],dTxr[1]};

				Nr[0]=dTxr[0];
				Nr[1]=dTxr[1];
				gauss(B,Nr,Nx);
				d[0]=-1*Nr[0];
				d[1]=-1*Nr[1];

//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;

				double Txr_min=Txr;
				double a_min=1e-3;

				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;
					double x0_a=x[0]+d[0]*alpha;
					double x1_a=x[1]+d[1]*alpha;

					double Txr_a=(x0_a-1)*(x0_a-1)+(x1_a-2)*(x1_a-2);
					if(x0_a*x0_a+x1_a*x1_a-1+theta0>0)	Txr_a+=0.5*r*(x0_a*x0_a+x1_a*x1_a-1+theta0)*(x0_a*x0_a+x1_a*x1_a-1+theta0);

					if(Txr_a<Txr_min)
					{
						Txr_min=Txr_a;
						a_min=alpha;
					}
				}
//				cout<<"Txr"<<count<<"="<<Txr_min<<", alpha="<<a_min<<endl;

				double d0=d[0]*a_min;
				double d1=d[1]*a_min;
				x[0]+=d0;
				x[1]+=d1;

//				cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<endl;


				fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
				dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

				c0=x[0]*x[0]+x[1]*x[1]-1;
				dc0[0]=2*x[0];	dc0[1]=2*x[1];

				Txr=fx;
				if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);
	
				//cout<<"Txr="<<Txr<<endl;

				for(int i=0;i<Nx;i++)
				{
					dTxr[i]=dfx[i];
					if(c0+theta0>0)	dTxr[i]+=r*(c0+theta0)*dc0[i];
				}
				
				E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
//				if(count%100==0)cout<<"E"<<count<<"="<<E<<endl;
				//cout<<endl;
				if(E<ep)	break;

				double s[2]={x[0]-x_k[0],x[1]-x_k[1]};
				double y[2]={dTxr[0]-dTxr_k[0],dTxr[1]-dTxr_k[1]};

				double beta=y[0]*s[0]+y[1]*s[1];
				double sigma=(s[0]*B[0]+s[1]*B[2])*s[0]+(s[0]*B[2]+s[1]*B[3])*s[1];

//				cout<<"beta"<<count<<"="<<beta<<", sigma"<<count<<"="<<sigma<<endl;

				if(beta>0)//(beta>=0.2*sigma)
				{
					double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
					double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

					B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
					B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
					B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
					B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
				}
	/*			else if(beta>0)
				{
					double theta=0.8*sigma/(sigma-beta);
					cout<<"theta="<<theta<<endl;
					double y0_s=y[0]*theta+(1-theta)*(B[0]*s[0]+B[1]*s[1]);
					double y1_s=y[1]*theta+(1-theta)*(B[2]*s[0]+B[3]*s[1]);

					beta=y0_s*s[0]+y1_s*s[1];
					double bs[2]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1], B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]};
					double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

					B[0*Nx+0]+=1/beta*y0_s*y0_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]);
					B[0*Nx+1]+=1/beta*y0_s*y1_s-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]);
					B[1*Nx+0]+=1/beta*y1_s*y0_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]);
					B[1*Nx+1]+=1/beta*y1_s*y1_s-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]);
				}*/
			}
		}
		if(calc_way==2)
		{
			fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
			dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

			c0=x[0]*x[0]+x[1]*x[1]-1;
			dc0[0]=2*x[0];	dc0[1]=2*x[1];

			Txr=fx;			
			if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);			
			dTxr[0]=dfx[0];
			dTxr[1]=dfx[1];
		
			if(c0+theta0>0)
			{
				dTxr[0]+=r*(c0+theta0)*dc0[0];
				dTxr[1]+=r*(c0+theta0)*dc0[1];
			}

			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
//			cout<<"E0="<<E<<endl;
			if(E<ep)	return 0;
			while(E>ep)
			{
				count++;
				stringstream se_min;
				se_min<<"./E"<<count_min<<".csv";
				ofstream init_emin(se_min.str(), ios::app);
				init_emin<<E<<endl;
				init_emin.close();					
		
				stringstream se_x;
				se_x<<"./x"<<count_min<<".csv";
				ofstream init_x(se_x.str(), ios::app);
				init_x<<x[0]<<","<<x[1]<<endl;
				init_x.close();
	
				stringstream se_fx;
				se_fx<<"./fx"<<count_min<<".csv";
				ofstream init_fx(se_fx.str(), ios::app);
				init_fx<<fx<<endl;
				init_fx.close();

				if(count>500)	break;

				d[0]=-1*dTxr[0];
				d[1]=-1*dTxr[1];

//				cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;

				double Txr_min=Txr;
				double a_min=1e-3;

				for(int i=0;i<1000;i++)
				{
					double alpha=(i+1)*1e-3;
					double x0_a=x[0]+d[0]*alpha;
					double x1_a=x[1]+d[1]*alpha;

					double Txr_a=(x0_a-1)*(x0_a-1)+(x1_a-2)*(x1_a-2);
					if(x0_a*x0_a+x1_a*x1_a-1+theta0>0)	Txr+=0.5*r*(x0_a*x0_a+x1_a*x1_a-1+theta0)*(x0_a*x0_a+x1_a*x1_a-1+theta0);

					if(Txr_a<Txr_min)
					{
						Txr_min=Txr_a;
						a_min=alpha;
					}
				}
//				cout<<"Txr"<<count<<"="<<Txr_min<<", alpha="<<a_min<<endl;

				double d0=d[0]*a_min;
				double d1=d[1]*a_min;
				x[0]+=d0;
				x[1]+=d1;

	//			cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<endl;


				fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
				dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);

				c0=x[0]*x[0]+x[1]*x[1]-1;
				dc0[0]=2*x[0];	dc0[1]=2*x[1];

				Txr=fx;
				if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);
	
				//cout<<"Txr="<<Txr<<endl;

				for(int i=0;i<Nx;i++)
				{
					dTxr[i]=dfx[i];
					if(c0+theta0>0)	dTxr[i]+=r*(c0+theta0)*dc0[i];
				}
				
				E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
//				if(count%100==0)	cout<<"E"<<count<<"="<<E<<endl;
				//cout<<endl;
			}
		}

		E_min=sqrt( (old_x[0]-x[0])*(old_x[0]-x[0]) + (old_x[1]-x[1])*(old_x[1]-x[1]));
		if(E_min<ep*1000)	r*=4;
		if(c0+theta0>0)	theta0+=c0;

		clock_t t1=clock();
		ofstream t_log("time_log.csv", ios::app);
		t_log<<count_min<<","<<(long double)(t1-t0)/CLOCKS_PER_SEC<<endl;
		t_log.close();

		ofstream femin("E.csv", ios::app);
		femin<<E<<","<<count<<endl;
		femin.close();					
		
		ofstream ffx("x.csv", ios::app);
		ffx<<x[0]<<","<<x[1]<<endl;
		ffx.close();
		ofstream fffx("fx.csv", ios::app);
		fffx<<fx<<endl;
		fffx.close();
	
		
		//		cout<<"E_min"<<count_min<<"="<<E_min<<endl<<endl;
	}

	//v0　計算
	int nv=0;
	cout<<"fc="<<c0<<endl;
	double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	cout<<"c="<<c0<<endl;
	cout<<"v="<<v0<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;

	ofstream fs("final_data.csv");
	fs<<"x0="<<x0<<","<<x1<<endl;
	fs<<"x="<<x[0]<<","<<x[1]<<endl;
	fs<<"c="<<c0<<endl;
	fs<<"v="<<v0<<endl;
	fs<<"fx="<<fx<<endl;
	fs<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/


	delete[]	x;
	delete[]	old_x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rTxr;
	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}

int MM_method()
{
	////例題		数理計画入門　システム制御情報ライブラリー15　問題（4.78）	p.138
	int Nx=4;
	double *x=new double [Nx];
	double *old_x=new double [Nx];
	x[0]=10;
	x[1]=-10;
	x[2]=10;
	x[3]=-10;


	double c0=10, c1=10,c2=10;
	double fx=0;
	double Txr=0;


	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rc1=new double [Nx*Nx];
	double *rc2=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rc1[i*Nx+j]=0;
			rc2[i*Nx+j]=0;
		}
	}


	double r=0.1;
	double ep=1e-10;
	double theta0=0, theta1=0, theta2=0;
	double E_min=1;
	int count_min=0;

	while(E_min>ep)
	{
		count_min++;
		for(int i=0;i<Nx;i++)	old_x[i]=x[i];

		double E=1;
		int count=0;
		while(E>ep)
		{
			count++;

			fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
			dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;
			rfx[0*Nx+0]=2;	rfx[1*Nx+1]=2;	rfx[2*Nx+2]=4;	rfx[3*Nx+3]=2;

			c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
			c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
			c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	
			dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
			dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
			dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;
			rc0[0*Nx+0]=2;	rc0[1*Nx+1]=2;	rc0[2*Nx+2]=2;	rc0[3*Nx+3]=2;
			rc1[0*Nx+0]=2;	rc1[1*Nx+1]=4;	rc1[2*Nx+2]=2;	rc1[3*Nx+3]=4;
			rc2[0*Nx+0]=4;	rc2[1*Nx+1]=2;	rc2[2*Nx+2]=2;	rc2[3*Nx+3]=0;
		
			Txr=fx;
			if(c0+theta0>0)	Txr+=0.5*r*(c0+theta0)*(c0+theta0);
			if(c1+theta1>0)	Txr+=0.5*r*(c1+theta1)*(c1+theta1);
			if(c2+theta2>0)	Txr+=0.5*r*(c2+theta2)*(c2+theta2);
	
			//cout<<"Txr="<<Txr<<endl;

			for(int i=0;i<Nx;i++)
			{
				dTxr[i]=dfx[i];
				if(c0+theta0>0)	dTxr[i]+=r*(c0+theta0)*dc0[i];
				if(c1+theta1>0)	dTxr[i]+=r*(c1+theta1)*dc1[i];
				if(c2+theta2>0)	dTxr[i]+=r*(c2+theta2)*dc2[i];
			}

			//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;

			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Nx;j++)
				{
					rTxr[i*Nx+j]=rfx[i*Nx+j];
					if(c0+theta0>0)	rTxr[i*Nx+j]+=r*dc0[j]*dc0[i]+r*(c0+theta0)*rc0[i*Nx+j];
					if(c1+theta1>0)	rTxr[i*Nx+j]+=r*dc1[j]*dc1[i]+r*(c1+theta1)*rc1[i*Nx+j];
					if(c2+theta2>0)	rTxr[i*Nx+j]+=r*dc2[j]*dc2[i]+r*(c2+theta2)*rc2[i*Nx+j];
				}
			}
			//cout<<"rTxr="<<rTxr[0]<<", "<<rTxr[1]<<", "<<rTxr[2]<<", "<<rTxr[3]<<endl;
			//cout<<rTxr[4]<<", "<<rTxr[5]<<", "<<rTxr[6]<<", "<<rTxr[7]<<endl;
			//cout<<rTxr[8]<<", "<<rTxr[9]<<", "<<rTxr[10]<<", "<<rTxr[11]<<endl;
			//cout<<rTxr[12]<<", "<<rTxr[13]<<", "<<rTxr[14]<<", "<<rTxr[15]<<endl;

			gauss(rTxr,dTxr,Nx);
			//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;
			for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];


			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]+dTxr[2]*dTxr[2]+dTxr[3]*dTxr[3]);
			cout<<"E"<<count<<"="<<E<<endl;
		}

		E_min=sqrt( (old_x[0]-x[0])*(old_x[0]-x[0]) + (old_x[1]-x[1])*(old_x[1]-x[1]) + (old_x[2]-x[2])*(old_x[2]-x[2])+ (old_x[3]-x[3])*(old_x[3]-x[3]));
		if(E_min<ep*100)	r*=4;
		if(c0+theta0>0)	theta0+=c0;
		if(c1+theta1>0)	theta1+=c1;
		if(c2+theta2>0)	theta2+=c2;

		cout<<"E_min"<<count_min<<"="<<E_min<<endl<<endl;
	}

	double *v=new double [3];
	v[0]=0;	v[1]=0;	v[2]=0;

	int *Nv_n=new int [3];
	Nv_n[0]=0;	Nv_n[1]=0;	Nv_n[2]=0;

	//v0　計算
	int nv=0;
	cout<<"fc="<<c0<<", "<<c1<<", "<<c2<<endl;
	if(c0>ep)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(c1>ep)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(c2>ep)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	cout<<"nv="<<nv<<endl;
	if(nv>1)
	{
		double  *Nrv=new double [nv];
		double  *Nlv=new double [nv*nv];

		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[i];
				else if(jv==1)	Nlv[i*nv+j]=dc1[i];
				else if(jv==2)	Nlv[i*nv+j]=dc2[i];
				cout<<"Nlv["<<i<<","<<j<<"]="<<Nlv[i*nv+j]<<endl;
			}
			//Nrv[i]=-1*(dTxr[i]+rTxr[i*Nx+0]*d[0]+rTxr[i*Nx+1]*d[1]+rTxr[i*Nx+2]*d[2]+rTxr[i*Nx+3]*d[3]);
			Nrv[i]=-dfx[i];
			cout<<"Nrv["<<i<<"]="<<Nrv[i]<<endl;
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
		delete[]	Nrv;
		delete[]	Nlv;
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	}
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;

	delete[]	v;
	delete[]	Nv_n;


	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;


	delete[]	x;
	delete[]	old_x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rc1;
	delete[]	rc2;
	delete[]	rTxr;


	return 0;
}

int MM_method2()
{
	////例題		システム工学第2版　森北出版（株）　演習問題5の4	p.197
	int Nx=2;
	double *x=new double [Nx];
	double *old_x=new double [Nx];
	x[0]=10;
	x[1]=-10;


	double c0=0;
	double fx=0;
	double Txr=0;


	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dTxr[i]=0;
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rTxr[i*Nx+j]=0;
		}
	}


	double r=1;
	double ep=1e-10;
	double theta0=0, theta1=0, theta2=0, theta3=0;
	double E_min=1;
	int count_min=0;

	while(E_min>ep)
	{
		count_min++;
		for(int i=0;i<Nx;i++)	old_x[i]=x[i];

		double E=1;
		int count=0;
		while(E>ep)
		{
			count++;

			fx=(x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
			dfx[0]=2*(x[0]-1);	dfx[1]=2*(x[1]-2);
			rfx[0*Nx+0]=2;	rfx[0*Nx+1]=0;	rfx[1*Nx+0]=0;	rfx[1*Nx+1]=2;

			c0=x[0]*x[0]+x[1]*x[1]-1;
			dc0[0]=2*x[0];	dc0[1]=2*x[1];
			rc0[0*Nx+0]=2;	rc0[0*Nx+1]=0;	rc0[1*Nx+0]=0;	rc0[1*Nx+1]=2;

			Txr=fx;
			Txr+=0.5*r*(c0+theta0)*(c0+theta0);
	
			//cout<<"Txr="<<Txr<<endl;

			for(int i=0;i<Nx;i++)
			{
				dTxr[i]=dfx[i];
				dTxr[i]+=r*(c0+theta0)*dc0[i];
			}

			//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;

			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Nx;j++)
				{
					rTxr[i*Nx+j]=rfx[i*Nx+j];
					rTxr[i*Nx+j]+=r*dc0[j]*dc0[i]+r*(c0+theta0)*rc0[i*Nx+j];
				}
			}
			//cout<<"rTxr="<<rTxr[0]<<", "<<rTxr[1]<<", "<<rTxr[2]<<", "<<rTxr[3]<<endl;
			//cout<<rTxr[4]<<", "<<rTxr[5]<<", "<<rTxr[6]<<", "<<rTxr[7]<<endl;
			//cout<<rTxr[8]<<", "<<rTxr[9]<<", "<<rTxr[10]<<", "<<rTxr[11]<<endl;
			//cout<<rTxr[12]<<", "<<rTxr[13]<<", "<<rTxr[14]<<", "<<rTxr[15]<<endl;

			gauss(rTxr,dTxr,Nx);
			//cout<<"dTxr="<<dTxr[0]<<", "<<dTxr[1]<<", "<<dTxr[2]<<", "<<dTxr[3]<<endl;
			for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];


			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
			cout<<"E"<<count<<"="<<E<<endl;
		}

		E_min=sqrt( (old_x[0]-x[0])*(old_x[0]-x[0]) + (old_x[1]-x[1])*(old_x[1]-x[1]));
		if(E_min<ep*1000)	r*=4;
		theta0+=c0;

		cout<<"E_min"<<count_min<<"="<<E_min<<endl<<endl;
	}

	//v0　計算
	int nv=0;
	cout<<"fc="<<c0<<endl;
	double v0=-1/( dc0[0]+dc0[1])*(dfx[0]+dfx[1]);

		//if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
		//else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dTxr[0]+dTxr[1]+dTxr[2]+dTxr[3] + rTxr[0*Nx+0]*d[0]+rTxr[0*Nx+1]*d[1]+rTxr[0*Nx+2]*d[2]+rTxr[0*Nx+3]*d[3]);
	cout<<"x="<<x[0]<<", "<<x[1]<<endl;
	cout<<"c="<<c0<<endl;
	cout<<"v="<<v0<<endl;
	cout<<"fx="<<fx<<endl;
	cout<<"r="<<r<<endl<<endl;


/*	double v0=1;
	double v1=0;
	double v2=2;
	x[0]=0;	x[1]=1;	x[2]=2;	x[3]=-1;
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

	cout<<"理論解"<<endl;
	cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"v="<<v0<<", "<<v1<<", "<<v2<<endl;
	cout<<"vc="<<c0*v0<<", "<<c1*v1<<", "<<c2*v2<<endl;*/


	delete[]	x;
	delete[]	old_x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rTxr;


	return 0;
}

int SUMT_IPM()
{
	int Nx=2;
	double *x=new double [Nx];
	x[0]=1;
	x[1]=5;

	double c0=10, c1=10,c2=10;
	double fx=0;
	double Txr=0;


	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rc1=new double [Nx*Nx];
	double *rc2=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rc1[i*Nx+j]=0;
			rc2[i*Nx+j]=0;
		}
	}


	double r=10.0;
	double ep=1e-10;
	
	while(r>ep/(-c0-c1-c2) && r>0)
	{
		
		double E=1;
		int count=0;

		while(E>ep)
		{
			count++;
			fx=3*x[0]+2*x[1];

			c0=x[0]+x[1]-8;
			c1=-7*x[0]-3*x[1]+21;
			c2=-2*x[0]-7*x[1]+14;	

			Txr=fx-r*(1/c0+1/c1+1/c2);
	
			dfx[0]=3;	dfx[1]=2;
	
			dc0[0]=1;	dc0[1]=1;
			dc1[0]=-7;	dc1[1]=-3;
			dc2[0]=-2;	dc2[1]=-7;
			
			for(int i=0;i<Nx;i++)
			{
				dTxr[i]=dfx[i]+r*(1/(c0*c0)*dc0[i]+1/(c1*c1)*dc1[i]+1/(c2*c2)*dc2[i]);
			}


			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Nx;j++)
				{
					rTxr[i*Nx+j]=rfx[i*Nx+j]+r*(-2/(c0*c0*c0)*dc0[i]*dc0[j]+1/(c0*c0)*rc0[i*Nx+j]-2/(c1*c1*c1)*dc1[i]*dc1[j]+1/(c1*c1)*rc1[i*Nx+j]-2/(c2*c2*c2)*dc2[i]*dc2[j]+1/(c2*c2)*rc2[i*Nx+j]);
				}
			}
			gauss(rTxr,dTxr,Nx);
			for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];

			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]);
			cout<<"E"<<count<<"="<<E<<endl;
		}
		cout<<"x="<<x[0]<<", "<<x[1]<<endl;
		cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
		cout<<"fx="<<fx<<endl;
		cout<<"r="<<r<<endl;

		r-=0.001;
	}

	delete[]	x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rc1;
	delete[]	rc2;
	delete[]	rTxr;


	/*
	int Nx=4;
	double *x=new double [Nx];
	x[0]=10;
	x[1]=-10;
	x[2]=10;
	x[3]=-10;

	double c0=10, c1=10, c2=10;
	double fx=0;
	double Txr=0;


	double *dfx=new double [Nx];
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];
	double *dTxr=new double [Nx];

	double *rfx=new double [Nx*Nx];
	double *rc0=new double [Nx*Nx];
	double *rc1=new double [Nx*Nx];
	double *rc2=new double [Nx*Nx];
	double *rTxr=new double [Nx*Nx];

	double r=10.0;
	double ep=1e-10;
	
	while(r>ep/(-c0-c1-c2) && r>0)
	{
		r-=0.1;
		
		double E=1;
		int count=0;

		while(E>ep)
		{
			count++;
			fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];

			c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
			c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
			c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	

			if(c0>0 || c1>0 || c2>0)
			{
				for(int i=0;i<Nx;i++)	x[i]+=dTxr[i];
				break;
			}

			Txr=fx-r*(1/c0+1/c1+1/c2);
	
			dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;
	
			dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
			dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
			dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;

			for(int i=0;i<Nx;i++)
			{
				dTxr[i]=dfx[i]+r*(1/(c0*c0)*dc0[i]+1/(c1*c1)*dc1[i]+1/(c2*c2)*dc2[i]);
			}

			rfx[0*Nx+0]=2;	rfx[1*Nx+1]=2;	rfx[2*Nx+2]=4;	rfx[3*Nx+3]=2;
			rc0[0*Nx+0]=2;	rc0[1*Nx+1]=2;	rc0[2*Nx+2]=2;	dc0[3*Nx+3]=2;
			rc1[0*Nx+0]=2;	rc1[1*Nx+1]=4;	rc1[2*Nx+2]=2;		dc1[3*Nx+3]=4;
			rc2[0*Nx+0]=4;	rc2[1*Nx+1]=2;	rc2[2*Nx+2]=2;		dc2[3*Nx+3]=0;

			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Nx;j++)
				{
					if(j!=i)
					{
						rfx[i*Nx+j]=0;
						rc0[i*Nx+j]=0;
						rc1[i*Nx+j]=0;
						rc2[i*Nx+j]=0;
					}
				}
			}

			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Nx;j++)
				{
					rTxr[i*Nx+j]=rfx[i*Nx+j]+r*(-2/(c0*c0*c0)*dc0[i]*dc0[j]+1/(c0*c0)*rc0[i*Nx+j]-2/(c1*c1*c1)*dc1[i]*dc1[j]+1/(c1*c1)*rc1[i*Nx+j]-2/(c2*c2*c2)*dc2[i]*dc2[j]+1/(c2*c2)*rc2[i*Nx+j]);
				}
			}
			gauss(rTxr,dTxr,Nx);
			for(int i=0;i<Nx;i++)	x[i]-=dTxr[i];

			E=sqrt(dTxr[0]*dTxr[0]+dTxr[1]*dTxr[1]+dTxr[2]*dTxr[2]+dTxr[3]*dTxr[3]);
			cout<<"E"<<count<<"="<<E<<endl;
		}
		cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
		cout<<"c="<<c0<<", "<<c1<<", "<<c2<<endl;
		cout<<"fx="<<fx<<endl;
		cout<<"r="<<r<<endl;

	}

	delete[]	x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	dTxr;
	delete[]	rfx;
	delete[]	rc0;
	delete[]	rc1;
	delete[]	rc2;
	delete[]	rTxr;
	*/
	return 0;
}

int SQP_method()
{
	
	

	int Nx=4;
	int Nv=3;

	double *x=new double [Nx];

	double fx=0;
	double *dfx=new double [Nx];

	double c0=0, c1=0, c2=0;
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];

	double *d=new double [Nx];

	double *dL=new double [Nx];
	double *B=new double [Nx*Nx];
	double *Nr=new double [Nx];
		
	double p=10.0;
	double Fp=0;
	double *dFp=new double[Nx];		

	double fc0=0, fc1=0, fc2=0;
	double *v=new double [Nv];
	int nv=0;
	int *Nv_n=new int [Nv];
	double *Nrv=new double [Nv];
	double *Nlv=new double [Nv*Nv];

	//初期化
	x[0]=10;
	x[1]=-10;
	x[2]=10;
	x[3]=-10;

	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dc1[i]=0;
		dc2[i]=0;
		for(int j=0;j<Nx;j++)
		{
			if(j==i)	B[i*Nx+j]=1;
			else
			{
				B[i*Nx+j]=0;
			}
		}
		d[i]=0;
		dL[i]=0;
		Nr[i]=0;
		dFp[i]=0;
	}
	for(int i=0;i<Nv;i++)
	{
		v[i]=0;
		Nv_n[i]=0;
		Nrv[i]=0;
		for(int j=0;j<Nv;j++)	Nlv[i*Nv+j]=0;
	}

	//fx関係　計算
	fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
	dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

	//c関係　計算
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	
	dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
	dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
	dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;
	
	for(int i=0;i<Nx;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];

	//Fp　計算
	Fp=fx;
	for(int i=0;i<Nx;i++)	dFp[i]+=dfx[i];
	if(c0>0)
	{
		Fp+=c0*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc0[i]*p;
	}
	else if(c0==0)
	{
		for(int i=0;i<Nx;i++)	if(dc0[i]>0)	dFp[i]+=dc0[i]*p;
	}
	if(c1>0)
	{
		Fp+=c1*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc1[i]*p;
	}
	else if(c1==0)
	{
		for(int i=0;i<Nx;i++)	if(dc1[i]>0)	dFp[i]+=dc1[i]*p;
	}
	if(c2>0)
	{
		Fp+=c2*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc2[i]*p;
	}
	else if(c2==0)
	{
		for(int i=0;i<Nx;i++)	if(dc2[i]>0)	dFp[i]+=dc2[i]*p;
	}

	double E=1;
	double ep=1e-10;

	//Nr[0]=dL[0];	Nr[1]=dL[1];	Nr[2]=dL[2];	Nr[3]=dL[3];	
	Nr[0]=dfx[0];	Nr[1]=dfx[1];	Nr[2]=dfx[2];	Nr[3]=dfx[3];	
	//Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
	//Nr[0]=dfx[0]+v[0]*dc0[0]+v[1]*dc1[0]+v[2]*dc2[0];	Nr[1]=dfx[1]+v[0]*dc0[1]+v[1]*dc1[1]+v[2]*dc2[1];	Nr[2]=dfx[2]+v[0]*dc0[2]+v[1]*dc1[2]+v[2]*dc2[2];	Nr[3]=dfx[3]+v[0]*dc0[3]+v[1]*dc1[3]+v[2]*dc2[3];	
	//qa																																																																																																																																																																																																																																																																																																																																	Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
	gauss(B,Nr,Nx);
	d[0]=-Nr[0];	d[1]=-Nr[1];	d[2]=-Nr[2];	d[3]=-Nr[3];	

	//v0　計算
	nv=0;
	fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
	fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
	fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
	if(fc0!=0)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(fc1!=0)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(fc2!=0)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	if(nv>1)
	{
		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[j];
				else if(jv==1)	Nlv[i*nv+j]=dc1[j];
				else if(jv==2)	Nlv[i*nv+j]=dc2[j];
			}
			Nrv[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
	}

	E=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
	//E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);

	//出力
	cout<<"x0="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c0="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"fc0="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
	cout<<"v0="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"Fp0="<<Fp<<endl;
	//cout<<"dFp0="<<dFp[0]<<", "<<dFp[1]<<", "<<dFp[2]<<", "<<dFp[3]<<endl;
	cout<<"dfx0="<<dfx[0]<<", "<<dfx[1]<<", "<<dfx[2]<<", "<<dfx[3]<<endl;
	cout<<"d0="<<d[0] <<", "<<d[1] <<", "<<d[2] <<", "<<d[3] <<endl;
	cout<<"d0/answer="<<d[0]/(-15) <<", "<<d[1]/25 <<", "<<d[2]/(-19) <<", "<<d[3]/13 <<endl;
	//cout<<"fx0="<<fx<<endl;
	cout<<"E0="<<E<<endl;


	if(E<ep)	return 0;
	else
	{
		int count=0;
		while(E>ep)
		{
			double x_k[4]={x[0], x[1], x[2], x[3]};
			double dFp_k[4]={dFp[0], dFp[1], dFp[2], dFp[3]};
			double dfx_k[4]={dfx[0], dfx[1], dfx[2], dfx[3]};
			double dc0_k[4]={dc0[0], dc0[1], dc0[2], dc0[3]};
			double dc1_k[4]={dc1[0], dc1[1], dc1[2], dc1[3]};
			double dc2_k[4]={dc2[0], dc2[1], dc2[2], dc2[3]};
			

			double Fp_min=Fp;
			double a_min=1;
			double x_a[4]={0,0,0,0};
			double fx_a=0;
			double c0_a=0, c1_a=0, c2_a=0;
			double Fp_a=0;
			if(count==0)	a_min=0.5;
			else
			{
				int N=0;
				double er=0;
				double r=10;
				while (r>ep)
				{
					r-=0.1;
					Fp_a=fx-r*(1/c0+1/c1+1/c2);
					if(Fp_a<Fp_min)
					{
						Nr=dfx;
						Fp_min=Fp_a;
				}
				for(int i=0;i<10000+N;i++)	
				{
					double alpha=(i+1)*1e-5;
					//double alpha=(i+1);
					x_a[0]=x[0]+d[0]*alpha; x_a[1]=x[1]+d[1]*alpha; x_a[2]=x[2]+d[2]*alpha; x_a[3]=x[3]+d[3]*alpha;
			
					c0_a=x_a[0]*x_a[0] +	x_a[1]*x_a[1] + x_a[2]*x_a[2] + x_a[3]*x_a[3] + x_a[0] -x_a[1] + x_a[2] -x_a[3] -8;
					c1_a=x_a[0]*x_a[0] +	2*x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[3]*x_a[3] -x_a[0] -x_a[3] -10;
					c2_a=2*x_a[0]*x_a[0] + x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[0] -x_a[1] -x_a[3] -5;	
					double dc0_a[4]={2*x_a[0]+1, 2*x_a[1]-1, 2*x_a[2]+1, 2*x_a[3]-1};
					double dc1_a[4]={2*x_a[0]-1, 4*x_a[1], 2*x_a[2], 4*x_a[3]-1};
					double dc2_a[4]={4*x_a[0]+2, 2*x_a[1]-1, 2*x_a[2], -1};//*/
					fx_a=x_a[0]*x_a[0] + x_a[1]*x_a[1] + 2*x_a[2]*x_a[2] + x_a[3]*x_a[3] -5*x_a[0] -5*x_a[1] -21*x_a[2] + 7*x_a[3];

					Fp_a=fx_a;
					if(c0_a>0)	Fp_a+=c0_a*p;
					if(c1_a>0)	Fp_a+=c1_a*p;
					if(c2_a>0)	Fp_a+=c2_a*p;

					double fc0_a=c0_a+dc0_a[0]*d[0]+dc0_a[1]*d[1]+dc0_a[2]*d[2]+dc0_a[3]*d[3];
					double fc1_a=c1_a+dc1_a[0]*d[0]+dc1_a[1]*d[1]+dc1_a[2]*d[2]+dc1_a[3]*d[3];
					double fc2_a=c2_a+dc2_a[0]*d[0]+dc2_a[1]*d[1]+dc2_a[2]*d[2]+dc2_a[3]*d[3];//*/
	//				cout<<"fc="<<fc0_a<<", "<<fc1_a<<", "<<fc2_a;
					if(fc0_a<=0 && fc1_a<=0 && fc2_a<=0)
					{
						if(Fp_a<Fp_min)
						{
							Fp_min=Fp_a;
							a_min=alpha;
						}
//						cout<<", Fp="<<Fp_a<<", a="<<alpha;
					}
			//		cout<<endl;
					//if(i==99+N && (fc0_a>0 || fc1_a>0 || fc2_a>0))	N+=100;
				}
			}
			cout<<"da"<<count<<"="<<d[0]*a_min<<", "<<d[1]*a_min<<", "<<d[2]*a_min<<", "<<d[3]*a_min<<endl<<endl;
						//Bに関する出力
			count++;
 			for(int i=0;i<Nx;i++)	x[i]+=d[i]*a_min;

			fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
			dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

			c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
			c1=x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
			c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;
			dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
			dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
			dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;

			Fp=fx;
			for(int i=0;i<Nx;i++)	dFp[i]+=dfx[i];
			if(c0>0)
			{
				Fp+=c0*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc0[i]*p;
			}
			else if(c0==0)
			{
				for(int i=0;i<Nx;i++)	if(dc0[i]>0)	dFp[i]+=dc0[i]*p;
			}
			if(c1>0)
			{
				Fp+=c1*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc1[i]*p;
			}
			else if(c1==0)
			{
				for(int i=0;i<Nx;i++)	if(dc1[i]>0)	dFp[i]+=dc1[i]*p;
			}
			if(c2>0)
			{
				Fp+=c2*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc2[i]*p;
			}
			else if(c2==0)
			{
				for(int i=0;i<Nx;i++)	if(dc2[i]>0)	dFp[i]+=dc2[i]*p;
			}

			//v0　計算
			nv=0;
			fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
			fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
			fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
			if(fc0!=0)	v[0]=0;
			else
			{
				Nv_n[nv]=0;
				nv++;
			}
			if(fc1!=0)	v[1]=0;
			else
			{
				Nv_n[nv]=1;
				nv++;
			}
			if(fc2!=0)	v[2]=0;
			else
			{
				Nv_n[nv]=2;
				nv++;
			}
			if(nv>1)
			{
				for(int i=0;i<nv;i++)
				{
					for(int j=0;j<nv;j++)
					{
						int jv=Nv_n[j];
						if(jv==0)	Nlv[i*nv+j]=dc0[j];
						else if(jv==1)	Nlv[i*nv+j]=dc1[j];
						else if(jv==2)	Nlv[i*nv+j]=dc2[j];
					}
					Nrv[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
				}
				gauss(Nlv,Nrv,nv);
				for(int i=0;i<nv;i++)
				{
					int iv=Nv_n[i];
					if(iv==0)	v[0]=Nrv[i];
					else if(iv==1)	v[1]=Nrv[i];
					else if(iv==2)	v[2]=Nrv[i];
				}
			}
			else if(nv==1)
			{
				int iv=Nv_n[nv-1];
				if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
				else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
				else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
			}

			//E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);
			E=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
			
			//出力
			cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
			cout<<"c"<<count<<"="<<c0<<", "<<c1<<", "<<c2<<endl;
			cout<<"v"<<count<<"="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
			cout<<"fc"<<count<<"="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
			cout<<"Fp"<<count<<"="<<Fp<<", alpha="<<a_min<<endl;
			//cout<<"fx"<<count<<"="<<f<<endl;
			cout<<"E"<<count<<"="<<E<<endl;
			
			if(E<ep)	break;

			for(int i=0;i<Nx;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];
			double dL_k[4]={dfx_k[0]+v[0]*dc0_k[0]+v[1]*dc1_k[0]+v[2]*dc2_k[0],
				dfx_k[1]+v[0]*dc0_k[1]+v[1]*dc1_k[1]+v[2]*dc2_k[1],
				dfx_k[2]+v[0]*dc0_k[2]+v[1]*dc1_k[2]+v[2]*dc2_k[2],
				dfx_k[3]+v[0]*dc0_k[3]+v[1]*dc1_k[3]+v[2]*dc2_k[3]};

			double s[4]={x[0]-x_k[0],x[1]-x_k[1], x[2]-x_k[2], x[3]-x_k[3]};
			//double y[4]={dFp[0]-dFp_k[0], dFp[1]-dFp_k[1], dFp[2]-dFp_k[2], dFp[3]-dFp_k[3]};
			double y[4]={dL[0]-dL_k[0], dL[1]-dL_k[1], dL[2]-dL_k[2], dL[3]-dL_k[3]};

			double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];

			double sigma=(s[0]*B[0*Nx+0]+s[1]*B[1*Nx+0]+s[2]*B[2*Nx+0]+s[3]*B[3*Nx+0])*s[0]
			+(s[0]*B[0*Nx+1]+s[1]*B[1*Nx+1]+s[2]*B[2*Nx+1]+s[3]*B[3*Nx+1])*s[1]
			+(s[0]*B[0*Nx+2]+s[1]*B[1*Nx+2]+s[2]*B[2*Nx+2]+s[3]*B[3*Nx+2])*s[2]
			+(s[0]*B[0*Nx+3]+s[1]*B[1*Nx+3]+s[2]*B[2*Nx+3]+s[3]*B[3*Nx+3])*s[3];
				
			double theta=1.0;
			if(beta>0)
			{
				if(beta<0.2*sigma)
				{
					theta=0.8*sigma/(sigma-beta);
					
					y[0]*=theta;
					y[1]*=theta;
					y[2]*=theta;
					y[3]*=theta;
					y[0]+=(1-theta)*(B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3]);
					y[1]+=(1-theta)*(B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3]);
					y[2]+=(1-theta)*(B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3]);
					y[3]+=(1-theta)*(B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3]);
				
					beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];
				}
				double bs[4]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3],
					B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3],
					B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3],
					B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3]};
				double bss[16]={bs[0]*s[0], bs[0]*s[1], bs[0]*s[2], bs[0]*s[3],
				bs[1]*s[0], bs[1]*s[1], bs[1]*s[2], bs[1]*s[3],
				bs[2]*s[0], bs[2]*s[1], bs[2]*s[2], bs[2]*s[3],
				bs[3]*s[0], bs[3]*s[1], bs[3]*s[2], bs[3]*s[3]};

				B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]+bss[0*Nx+2]*B[2*Nx+0]+bss[0*Nx+3]*B[3*Nx+0]);
				B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]+bss[0*Nx+2]*B[2*Nx+1]+bss[0*Nx+3]*B[3*Nx+1]);
				B[0*Nx+2]+=1/beta*y[0]*y[2]-1/sigma*(bss[0*Nx+0]*B[0*Nx+2]+bss[0*Nx+1]*B[1*Nx+2]+bss[0*Nx+2]*B[2*Nx+2]+bss[0*Nx+3]*B[3*Nx+2]);
				B[0*Nx+3]+=1/beta*y[0]*y[3]-1/sigma*(bss[0*Nx+0]*B[0*Nx+3]+bss[0*Nx+1]*B[1*Nx+3]+bss[0*Nx+2]*B[2*Nx+3]+bss[0*Nx+3]*B[3*Nx+3]);

				B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]+bss[1*Nx+2]*B[2*Nx+0]+bss[1*Nx+3]*B[3*Nx+0]);
				B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]+bss[1*Nx+2]*B[2*Nx+1]+bss[1*Nx+3]*B[3*Nx+1]);
				B[1*Nx+2]+=1/beta*y[1]*y[2]-1/sigma*(bss[1*Nx+0]*B[0*Nx+2]+bss[1*Nx+1]*B[1*Nx+2]+bss[1*Nx+2]*B[2*Nx+2]+bss[1*Nx+3]*B[3*Nx+2]);
				B[1*Nx+3]+=1/beta*y[1]*y[3]-1/sigma*(bss[1*Nx+0]*B[0*Nx+3]+bss[1*Nx+1]*B[1*Nx+3]+bss[1*Nx+2]*B[2*Nx+3]+bss[1*Nx+3]*B[3*Nx+3]);

				B[2*Nx+0]+=1/beta*y[2]*y[0]-1/sigma*(bss[2*Nx+0]*B[0*Nx+0]+bss[2*Nx+1]*B[1*Nx+0]+bss[2*Nx+2]*B[2*Nx+0]+bss[2*Nx+3]*B[3*Nx+0]);
				B[2*Nx+1]+=1/beta*y[2]*y[1]-1/sigma*(bss[2*Nx+0]*B[0*Nx+1]+bss[2*Nx+1]*B[1*Nx+1]+bss[2*Nx+2]*B[2*Nx+1]+bss[2*Nx+3]*B[3*Nx+1]);
				B[2*Nx+2]+=1/beta*y[2]*y[2]-1/sigma*(bss[2*Nx+0]*B[0*Nx+2]+bss[2*Nx+1]*B[1*Nx+2]+bss[2*Nx+2]*B[2*Nx+2]+bss[2*Nx+3]*B[3*Nx+2]);
				B[2*Nx+3]+=1/beta*y[2]*y[3]-1/sigma*(bss[2*Nx+0]*B[0*Nx+3]+bss[2*Nx+1]*B[1*Nx+3]+bss[2*Nx+2]*B[2*Nx+3]+bss[2*Nx+3]*B[3*Nx+3]);

				B[3*Nx+0]+=1/beta*y[3]*y[0]-1/sigma*(bss[3*Nx+0]*B[0*Nx+0]+bss[3*Nx+1]*B[1*Nx+0]+bss[3*Nx+2]*B[2*Nx+0]+bss[3*Nx+3]*B[3*Nx+0]);
				B[3*Nx+1]+=1/beta*y[3]*y[1]-1/sigma*(bss[3*Nx+0]*B[0*Nx+1]+bss[3*Nx+1]*B[1*Nx+1]+bss[3*Nx+2]*B[2*Nx+1]+bss[3*Nx+3]*B[3*Nx+1]);
				B[3*Nx+2]+=1/beta*y[3]*y[2]-1/sigma*(bss[3*Nx+0]*B[0*Nx+2]+bss[3*Nx+1]*B[1*Nx+2]+bss[3*Nx+2]*B[2*Nx+2]+bss[3*Nx+3]*B[3*Nx+2]);
				B[3*Nx+3]+=1/beta*y[3]*y[3]-1/sigma*(bss[3*Nx+0]*B[0*Nx+3]+bss[3*Nx+1]*B[1*Nx+3]+bss[3*Nx+2]*B[2*Nx+3]+bss[3*Nx+3]*B[3*Nx+3]);
			}			
			cout<<"beta"<<count<<"="<<beta<<endl;

			//Nr[0]=dL[0];	Nr[1]=dL[1];	Nr[2]=dL[2];	Nr[3]=dL[3];	
			Nr[0]=dfx[0];	Nr[1]=dfx[1];	Nr[2]=dfx[2];	Nr[3]=dfx[3];	
			//Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
			//Nr[0]=dfx[0]+v[0]*dc0[0]+v[1]*dc1[0]+v[2]*dc2[0];	Nr[1]=dfx[1]+v[0]*dc0[1]+v[1]*dc1[1]+v[2]*dc2[1];	Nr[2]=dfx[2]+v[0]*dc0[2]+v[1]*dc1[2]+v[2]*dc2[2];	Nr[3]=dfx[3]+v[0]*dc0[3]+v[1]*dc1[3]+v[2]*dc2[3];	
			gauss(B,Nr,Nx);
			d[0]=-Nr[0];	d[1]=-Nr[1];	d[2]=-Nr[2];	d[3]=-Nr[3];	
			
			cout<<"d"<<count<<"="<<d[0] <<", "<<d[1] <<", "<<d[2] <<", "<<d[3] <<endl;
			cout<<"d1/answer="<<d[0]/(4.12442) <<", "<<d[1]/(-6.87404) <<", "<<d[2]/(20.28458) <<", "<<d[3]/(-3.57450) <<endl;
		}
	}
	}

	/*
	int Nx=4;
	int Nv=3;

	double *x=new double [Nx];

	double fx=0;
	double *dfx=new double [Nx];

	double c0=0, c1=0, c2=0;
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];

	double *d=new double [Nx];

	double *dL=new double [Nx];
	double *B=new double [Nx*Nx];
	double *Nr=new double [Nx];
		
	double p=10.0;
	double Fp=0;
	double *dFp=new double[Nx];		

	double fc0=0, fc1=0, fc2=0;
	double *v=new double [Nv];
	int nv=0;
	int *Nv_n=new int [Nv];
	double *Nrv=new double [Nv];
	double *Nlv=new double [Nv*Nv];

	//初期化
	x[0]=10;
	x[1]=-10;
	x[2]=10;
	x[3]=-10;

	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dc1[i]=0;
		dc2[i]=0;
		for(int j=0;j<Nx;j++)
		{
			if(j==i)	B[i*Nx+j]=1;
			else
			{
				B[i*Nx+j]=0;
			}
		}
		d[i]=0;
		dL[i]=0;
		Nr[i]=0;
		dFp[i]=0;
	}
	for(int i=0;i<Nv;i++)
	{
		v[i]=0;
		Nv_n[i]=0;
		Nrv[i]=0;
		for(int j=0;j<Nv;j++)	Nlv[i*Nv+j]=0;
	}

	//fx関係　計算
	fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
	dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

	//c関係　計算
	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
	c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;	
	dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
	dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
	dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;
	
	for(int i=0;i<Nx;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];

	//Fp　計算
	Fp=fx;
	for(int i=0;i<Nx;i++)	dFp[i]+=dfx[i];
	if(c0>0)
	{
		Fp+=c0*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc0[i]*p;
	}
	else if(c0==0)
	{
		for(int i=0;i<Nx;i++)	if(dc0[i]>0)	dFp[i]+=dc0[i]*p;
	}
	if(c1>0)
	{
		Fp+=c1*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc1[i]*p;
	}
	else if(c1==0)
	{
		for(int i=0;i<Nx;i++)	if(dc1[i]>0)	dFp[i]+=dc1[i]*p;
	}
	if(c2>0)
	{
		Fp+=c2*p;
		for(int i=0;i<Nx;i++)	dFp[i]+=dc2[i]*p;
	}
	else if(c2==0)
	{
		for(int i=0;i<Nx;i++)	if(dc2[i]>0)	dFp[i]+=dc2[i]*p;
	}

	double E=1;
	double ep=1e-10;

	//Nr[0]=dL[0];	Nr[1]=dL[1];	Nr[2]=dL[2];	Nr[3]=dL[3];	
	Nr[0]=dfx[0];	Nr[1]=dfx[1];	Nr[2]=dfx[2];	Nr[3]=dfx[3];	
	//Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
	//Nr[0]=dfx[0]+v[0]*dc0[0]+v[1]*dc1[0]+v[2]*dc2[0];	Nr[1]=dfx[1]+v[0]*dc0[1]+v[1]*dc1[1]+v[2]*dc2[1];	Nr[2]=dfx[2]+v[0]*dc0[2]+v[1]*dc1[2]+v[2]*dc2[2];	Nr[3]=dfx[3]+v[0]*dc0[3]+v[1]*dc1[3]+v[2]*dc2[3];	
	//qa																																																																																																																																																																																																																																																																																																																																	Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
	gauss(B,Nr,Nx);
	d[0]=-Nr[0];	d[1]=-Nr[1];	d[2]=-Nr[2];	d[3]=-Nr[3];	

	//v0　計算
	nv=0;
	fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
	fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
	fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
	if(fc0!=0)	v[0]=0;
	else
	{
		Nv_n[nv]=0;
		nv++;
	}
	if(fc1!=0)	v[1]=0;
	else
	{
		Nv_n[nv]=1;
		nv++;
	}
	if(fc2!=0)	v[2]=0;
	else
	{
		Nv_n[nv]=2;
		nv++;
	}
	if(nv>1)
	{
		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nlv[i*nv+j]=dc0[j];
				else if(jv==1)	Nlv[i*nv+j]=dc1[j];
				else if(jv==2)	Nlv[i*nv+j]=dc2[j];
			}
			Nrv[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
		}
		gauss(Nlv,Nrv,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nrv[i];
			else if(iv==1)	v[1]=Nrv[i];
			else if(iv==2)	v[2]=Nrv[i];
		}
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
	}

	E=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
	//E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);

	//出力
	cout<<"x0="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c0="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"fc0="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
	cout<<"v0="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"Fp0="<<Fp<<endl;
	//cout<<"dFp0="<<dFp[0]<<", "<<dFp[1]<<", "<<dFp[2]<<", "<<dFp[3]<<endl;
	cout<<"dfx0="<<dfx[0]<<", "<<dfx[1]<<", "<<dfx[2]<<", "<<dfx[3]<<endl;
	cout<<"d0="<<d[0] <<", "<<d[1] <<", "<<d[2] <<", "<<d[3] <<endl;
	cout<<"d0/answer="<<d[0]/(-15) <<", "<<d[1]/25 <<", "<<d[2]/(-19) <<", "<<d[3]/13 <<endl;
	//cout<<"fx0="<<fx<<endl;
	cout<<"E0="<<E<<endl;


	if(E<ep)	return 0;
	else
	{
		int count=0;
		while(E>ep)
		{
			double x_k[4]={x[0], x[1], x[2], x[3]};
			double dFp_k[4]={dFp[0], dFp[1], dFp[2], dFp[3]};
			double dfx_k[4]={dfx[0], dfx[1], dfx[2], dfx[3]};
			double dc0_k[4]={dc0[0], dc0[1], dc0[2], dc0[3]};
			double dc1_k[4]={dc1[0], dc1[1], dc1[2], dc1[3]};
			double dc2_k[4]={dc2[0], dc2[1], dc2[2], dc2[3]};
			

			double Fp_min=Fp;
			double a_min=1;
			double x_a[4]={0,0,0,0};
			double fx_a=0;
			double c0_a=0, c1_a=0, c2_a=0;
			double Fp_a=0;
			if(count==0)	a_min=0.5;
			else
			{
				int N=0;
				for(int i=0;i<10000+N;i++)	
				{
					double alpha=(i+1)*1e-5;
					//double alpha=(i+1);
					x_a[0]=x[0]+d[0]*alpha; x_a[1]=x[1]+d[1]*alpha; x_a[2]=x[2]+d[2]*alpha; x_a[3]=x[3]+d[3]*alpha;
			
					c0_a=x_a[0]*x_a[0] +	x_a[1]*x_a[1] + x_a[2]*x_a[2] + x_a[3]*x_a[3] + x_a[0] -x_a[1] + x_a[2] -x_a[3] -8;
					c1_a=x_a[0]*x_a[0] +	2*x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[3]*x_a[3] -x_a[0] -x_a[3] -10;
					c2_a=2*x_a[0]*x_a[0] + x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[0] -x_a[1] -x_a[3] -5;	
					double dc0_a[4]={2*x_a[0]+1, 2*x_a[1]-1, 2*x_a[2]+1, 2*x_a[3]-1};
					double dc1_a[4]={2*x_a[0]-1, 4*x_a[1], 2*x_a[2], 4*x_a[3]-1};
					double dc2_a[4]={4*x_a[0]+2, 2*x_a[1]-1, 2*x_a[2], -1};//
					fx_a=x_a[0]*x_a[0] + x_a[1]*x_a[1] + 2*x_a[2]*x_a[2] + x_a[3]*x_a[3] -5*x_a[0] -5*x_a[1] -21*x_a[2] + 7*x_a[3];

					Fp_a=fx_a;
					if(c0_a>0)	Fp_a+=c0_a*p;
					if(c1_a>0)	Fp_a+=c1_a*p;
					if(c2_a>0)	Fp_a+=c2_a*p;

					double fc0_a=c0_a+dc0_a[0]*d[0]+dc0_a[1]*d[1]+dc0_a[2]*d[2]+dc0_a[3]*d[3];
					double fc1_a=c1_a+dc1_a[0]*d[0]+dc1_a[1]*d[1]+dc1_a[2]*d[2]+dc1_a[3]*d[3];
					double fc2_a=c2_a+dc2_a[0]*d[0]+dc2_a[1]*d[1]+dc2_a[2]*d[2]+dc2_a[3]*d[3];//
	//				cout<<"fc="<<fc0_a<<", "<<fc1_a<<", "<<fc2_a;
					if(fc0_a<=0 && fc1_a<=0 && fc2_a<=0)
					{
						if(Fp_a<Fp_min)
						{
							Fp_min=Fp_a;
							a_min=alpha;
						}
//						cout<<", Fp="<<Fp_a<<", a="<<alpha;
					}
			//		cout<<endl;
					//if(i==99+N && (fc0_a>0 || fc1_a>0 || fc2_a>0))	N+=100;
				}
			}
			cout<<"da"<<count<<"="<<d[0]*a_min<<", "<<d[1]*a_min<<", "<<d[2]*a_min<<", "<<d[3]*a_min<<endl<<endl;
						//Bに関する出力
			count++;
 			for(int i=0;i<Nx;i++)	x[i]+=d[i]*a_min;

			fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
			dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

			c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
			c1=x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
			c2=2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] -x[1] -x[3] -5;
			dc0[0]=2*x[0]+1;	dc0[1]=2*x[1]-1;	dc0[2]=2*x[2]+1;	dc0[3]=2*x[3]-1;
			dc1[0]=2*x[0]-1;	dc1[1]=4*x[1];		dc1[2]=2*x[2];		dc1[3]=4*x[3]-1;
			dc2[0]=4*x[0]+2;	dc2[1]=2*x[1]-1;	dc2[2]=2*x[2];		dc2[3]=-1;

			Fp=fx;
			for(int i=0;i<Nx;i++)	dFp[i]+=dfx[i];
			if(c0>0)
			{
				Fp+=c0*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc0[i]*p;
			}
			else if(c0==0)
			{
				for(int i=0;i<Nx;i++)	if(dc0[i]>0)	dFp[i]+=dc0[i]*p;
			}
			if(c1>0)
			{
				Fp+=c1*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc1[i]*p;
			}
			else if(c1==0)
			{
				for(int i=0;i<Nx;i++)	if(dc1[i]>0)	dFp[i]+=dc1[i]*p;
			}
			if(c2>0)
			{
				Fp+=c2*p;
				for(int i=0;i<Nx;i++)	dFp[i]+=dc2[i]*p;
			}
			else if(c2==0)
			{
				for(int i=0;i<Nx;i++)	if(dc2[i]>0)	dFp[i]+=dc2[i]*p;
			}

			//v0　計算
			nv=0;
			fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
			fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
			fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
			if(fc0!=0)	v[0]=0;
			else
			{
				Nv_n[nv]=0;
				nv++;
			}
			if(fc1!=0)	v[1]=0;
			else
			{
				Nv_n[nv]=1;
				nv++;
			}
			if(fc2!=0)	v[2]=0;
			else
			{
				Nv_n[nv]=2;
				nv++;
			}
			if(nv>1)
			{
				for(int i=0;i<nv;i++)
				{
					for(int j=0;j<nv;j++)
					{
						int jv=Nv_n[j];
						if(jv==0)	Nlv[i*nv+j]=dc0[j];
						else if(jv==1)	Nlv[i*nv+j]=dc1[j];
						else if(jv==2)	Nlv[i*nv+j]=dc2[j];
					}
					Nrv[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
				}
				gauss(Nlv,Nrv,nv);
				for(int i=0;i<nv;i++)
				{
					int iv=Nv_n[i];
					if(iv==0)	v[0]=Nrv[i];
					else if(iv==1)	v[1]=Nrv[i];
					else if(iv==2)	v[2]=Nrv[i];
				}
			}
			else if(nv==1)
			{
				int iv=Nv_n[nv-1];
				if(iv==0)	v[0]=-1/( dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
				else if(iv==1)	v[1]=-1/( dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
				else if(iv==2)	v[2]=-1/( dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
			}

			//E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);
			E=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3]);
			
			//出力
			cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
			cout<<"c"<<count<<"="<<c0<<", "<<c1<<", "<<c2<<endl;
			cout<<"v"<<count<<"="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
			cout<<"fc"<<count<<"="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
			cout<<"Fp"<<count<<"="<<Fp<<", alpha="<<a_min<<endl;
			//cout<<"fx"<<count<<"="<<f<<endl;
			cout<<"E"<<count<<"="<<E<<endl;
			
			if(E<ep)	break;

			for(int i=0;i<Nx;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];
			double dL_k[4]={dfx_k[0]+v[0]*dc0_k[0]+v[1]*dc1_k[0]+v[2]*dc2_k[0],
				dfx_k[1]+v[0]*dc0_k[1]+v[1]*dc1_k[1]+v[2]*dc2_k[1],
				dfx_k[2]+v[0]*dc0_k[2]+v[1]*dc1_k[2]+v[2]*dc2_k[2],
				dfx_k[3]+v[0]*dc0_k[3]+v[1]*dc1_k[3]+v[2]*dc2_k[3]};

			double s[4]={x[0]-x_k[0],x[1]-x_k[1], x[2]-x_k[2], x[3]-x_k[3]};
			//double y[4]={dFp[0]-dFp_k[0], dFp[1]-dFp_k[1], dFp[2]-dFp_k[2], dFp[3]-dFp_k[3]};
			double y[4]={dL[0]-dL_k[0], dL[1]-dL_k[1], dL[2]-dL_k[2], dL[3]-dL_k[3]};

			double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];

			double sigma=(s[0]*B[0*Nx+0]+s[1]*B[1*Nx+0]+s[2]*B[2*Nx+0]+s[3]*B[3*Nx+0])*s[0]
			+(s[0]*B[0*Nx+1]+s[1]*B[1*Nx+1]+s[2]*B[2*Nx+1]+s[3]*B[3*Nx+1])*s[1]
			+(s[0]*B[0*Nx+2]+s[1]*B[1*Nx+2]+s[2]*B[2*Nx+2]+s[3]*B[3*Nx+2])*s[2]
			+(s[0]*B[0*Nx+3]+s[1]*B[1*Nx+3]+s[2]*B[2*Nx+3]+s[3]*B[3*Nx+3])*s[3];
				
			double theta=1.0;
			if(beta>0)
			{
				if(beta<0.2*sigma)
				{
					theta=0.8*sigma/(sigma-beta);
					
					y[0]*=theta;
					y[1]*=theta;
					y[2]*=theta;
					y[3]*=theta;
					y[0]+=(1-theta)*(B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3]);
					y[1]+=(1-theta)*(B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3]);
					y[2]+=(1-theta)*(B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3]);
					y[3]+=(1-theta)*(B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3]);
				
					beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];
				}
				double bs[4]={B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3],
					B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3],
					B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3],
					B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3]};
				double bss[16]={bs[0]*s[0], bs[0]*s[1], bs[0]*s[2], bs[0]*s[3],
				bs[1]*s[0], bs[1]*s[1], bs[1]*s[2], bs[1]*s[3],
				bs[2]*s[0], bs[2]*s[1], bs[2]*s[2], bs[2]*s[3],
				bs[3]*s[0], bs[3]*s[1], bs[3]*s[2], bs[3]*s[3]};

				B[0*Nx+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*Nx+0]*B[0*Nx+0]+bss[0*Nx+1]*B[1*Nx+0]+bss[0*Nx+2]*B[2*Nx+0]+bss[0*Nx+3]*B[3*Nx+0]);
				B[0*Nx+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*Nx+0]*B[0*Nx+1]+bss[0*Nx+1]*B[1*Nx+1]+bss[0*Nx+2]*B[2*Nx+1]+bss[0*Nx+3]*B[3*Nx+1]);
				B[0*Nx+2]+=1/beta*y[0]*y[2]-1/sigma*(bss[0*Nx+0]*B[0*Nx+2]+bss[0*Nx+1]*B[1*Nx+2]+bss[0*Nx+2]*B[2*Nx+2]+bss[0*Nx+3]*B[3*Nx+2]);
				B[0*Nx+3]+=1/beta*y[0]*y[3]-1/sigma*(bss[0*Nx+0]*B[0*Nx+3]+bss[0*Nx+1]*B[1*Nx+3]+bss[0*Nx+2]*B[2*Nx+3]+bss[0*Nx+3]*B[3*Nx+3]);

				B[1*Nx+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*Nx+0]*B[0*Nx+0]+bss[1*Nx+1]*B[1*Nx+0]+bss[1*Nx+2]*B[2*Nx+0]+bss[1*Nx+3]*B[3*Nx+0]);
				B[1*Nx+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*Nx+0]*B[0*Nx+1]+bss[1*Nx+1]*B[1*Nx+1]+bss[1*Nx+2]*B[2*Nx+1]+bss[1*Nx+3]*B[3*Nx+1]);
				B[1*Nx+2]+=1/beta*y[1]*y[2]-1/sigma*(bss[1*Nx+0]*B[0*Nx+2]+bss[1*Nx+1]*B[1*Nx+2]+bss[1*Nx+2]*B[2*Nx+2]+bss[1*Nx+3]*B[3*Nx+2]);
				B[1*Nx+3]+=1/beta*y[1]*y[3]-1/sigma*(bss[1*Nx+0]*B[0*Nx+3]+bss[1*Nx+1]*B[1*Nx+3]+bss[1*Nx+2]*B[2*Nx+3]+bss[1*Nx+3]*B[3*Nx+3]);

				B[2*Nx+0]+=1/beta*y[2]*y[0]-1/sigma*(bss[2*Nx+0]*B[0*Nx+0]+bss[2*Nx+1]*B[1*Nx+0]+bss[2*Nx+2]*B[2*Nx+0]+bss[2*Nx+3]*B[3*Nx+0]);
				B[2*Nx+1]+=1/beta*y[2]*y[1]-1/sigma*(bss[2*Nx+0]*B[0*Nx+1]+bss[2*Nx+1]*B[1*Nx+1]+bss[2*Nx+2]*B[2*Nx+1]+bss[2*Nx+3]*B[3*Nx+1]);
				B[2*Nx+2]+=1/beta*y[2]*y[2]-1/sigma*(bss[2*Nx+0]*B[0*Nx+2]+bss[2*Nx+1]*B[1*Nx+2]+bss[2*Nx+2]*B[2*Nx+2]+bss[2*Nx+3]*B[3*Nx+2]);
				B[2*Nx+3]+=1/beta*y[2]*y[3]-1/sigma*(bss[2*Nx+0]*B[0*Nx+3]+bss[2*Nx+1]*B[1*Nx+3]+bss[2*Nx+2]*B[2*Nx+3]+bss[2*Nx+3]*B[3*Nx+3]);

				B[3*Nx+0]+=1/beta*y[3]*y[0]-1/sigma*(bss[3*Nx+0]*B[0*Nx+0]+bss[3*Nx+1]*B[1*Nx+0]+bss[3*Nx+2]*B[2*Nx+0]+bss[3*Nx+3]*B[3*Nx+0]);
				B[3*Nx+1]+=1/beta*y[3]*y[1]-1/sigma*(bss[3*Nx+0]*B[0*Nx+1]+bss[3*Nx+1]*B[1*Nx+1]+bss[3*Nx+2]*B[2*Nx+1]+bss[3*Nx+3]*B[3*Nx+1]);
				B[3*Nx+2]+=1/beta*y[3]*y[2]-1/sigma*(bss[3*Nx+0]*B[0*Nx+2]+bss[3*Nx+1]*B[1*Nx+2]+bss[3*Nx+2]*B[2*Nx+2]+bss[3*Nx+3]*B[3*Nx+2]);
				B[3*Nx+3]+=1/beta*y[3]*y[3]-1/sigma*(bss[3*Nx+0]*B[0*Nx+3]+bss[3*Nx+1]*B[1*Nx+3]+bss[3*Nx+2]*B[2*Nx+3]+bss[3*Nx+3]*B[3*Nx+3]);
			}			
			cout<<"beta"<<count<<"="<<beta<<endl;

			//Nr[0]=dL[0];	Nr[1]=dL[1];	Nr[2]=dL[2];	Nr[3]=dL[3];	
			Nr[0]=dfx[0];	Nr[1]=dfx[1];	Nr[2]=dfx[2];	Nr[3]=dfx[3];	
			//Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
			//Nr[0]=dfx[0]+v[0]*dc0[0]+v[1]*dc1[0]+v[2]*dc2[0];	Nr[1]=dfx[1]+v[0]*dc0[1]+v[1]*dc1[1]+v[2]*dc2[1];	Nr[2]=dfx[2]+v[0]*dc0[2]+v[1]*dc1[2]+v[2]*dc2[2];	Nr[3]=dfx[3]+v[0]*dc0[3]+v[1]*dc1[3]+v[2]*dc2[3];	
			gauss(B,Nr,Nx);
			d[0]=-Nr[0];	d[1]=-Nr[1];	d[2]=-Nr[2];	d[3]=-Nr[3];	
			
			cout<<"d"<<count<<"="<<d[0] <<", "<<d[1] <<", "<<d[2] <<", "<<d[3] <<endl;
			cout<<"d1/answer="<<d[0]/(4.12442) <<", "<<d[1]/(-6.87404) <<", "<<d[2]/(20.28458) <<", "<<d[3]/(-3.57450) <<endl;
		}
	}

	*/

	delete[]	x;
	delete[]	dfx;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	d;
	delete[]	B;
	delete[]	Nr;
	delete[]	dFp;
	delete[]	v;
	delete[]	Nv_n;
	delete[]	Nrv;
	delete[]	Nlv;
	return 0;
}

int q_newton()
{
	int N=2;
	double *x=new double [N];
	double *df=new double [N];
	double *B=new double[N*N];
	double *d=new double [N];	
	double *Nr=new double [N];

	for(int i=0;i<N;i++)
	{
		x[i]=0;
		df[i]=0;
		d[i]=0;
		Nr[i]=0;
		for(int j=0;j<N;j++)
		{
			
			if(j==i)	B[i*N+j]=1;
			else
			{
				B[i*N+j]=0;
			}
		}

	}

	//x_kとB_kの初期設定

	double f=(x[0]-1)*(x[0]-1)+10*(x[0]*x[0]-x[1])*(x[0]*x[0]-x[1]);

	//df_kの初期設定
	df[0]=2*(x[0]-1)+20*(x[0]*x[0]-x[1])*2*x[0];
	df[1]=-20*(x[0]*x[0]-x[1]);

	double E=1;
	double ep=1e-10;
	
	////df_k=0なら計算終了
	E=sqrt(df[0]*df[0]+df[1]*df[1]);
	cout<<"E0="<<E<<endl;
	if(E<ep)	return 0;

	else
	{
		int count=0;
		while(E>ep)
		{
			count++;

			double x_k[2]={x[0], x[1]};
			double df_k[2]={df[0],df[1]};

			Nr[0]=df[0];
			Nr[1]=df[1];
			gauss(B, Nr,N);
			d[0]=-Nr[0];
			d[1]=-Nr[1];

			cout<<"d"<<count<<"="<<d[0]<<", "<<d[1]<<endl;

			double f_min=f;
			double a_min=0;

			for(int i=0;i<1000;i++)
			{
				double alpha=(i+1)*1e-3;
				double x0_a=x[0]+d[0]*alpha;
				double x1_a=x[1]+d[1]*alpha;
				double f_a=(x0_a-1)*(x0_a-1)+10*(x0_a*x0_a-x1_a)*(x0_a*x0_a-x1_a);
				if(f_a<f_min)
				{
					f_min=f_a;
					a_min=alpha;
				}
			}
			cout<<"f"<<count<<"="<<f_min<<", alpha="<<a_min<<endl;

			double d0=d[0]*a_min;
			double d1=d[1]*a_min;
			x[0]+=d0;
			x[1]+=d1;

			cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<endl;

			f=(x[0]-1)*(x[0]-1)+10*(x[0]*x[0]-x[1])*(x[0]*x[0]-x[1]);
			df[0]=2*(x[0]-1)+20*(x[0]*x[0]-x[1])*2*x[0];
			df[1]=-20*(x[0]*x[0]-x[1]);
			cout<<"f"<<count<<"="<<f<<endl;

			E=sqrt(df[0]*df[0]+df[1]*df[1]);
			cout<<"E"<<count<<"="<<E<<endl;
			cout<<endl;

			if(E<ep)	break;

			double s[2]={x[0]-x_k[0],x[1]-x_k[1]};
			double y[2]={df[0]-df_k[0],df[1]-df_k[1]};

			double beta=y[0]*s[0]+y[1]*s[1];

			cout<<"beta"<<count<<"="<<beta<<endl;

			if(beta>0)
			{

				double sigma=(s[0]*B[0]+s[1]*B[2])*s[0]+(s[0]*B[2]+s[1]*B[3])*s[1];
				double bs[2]={B[0*N+0]*s[0]+B[0*N+1]*s[1], B[1*N+0]*s[0]+B[1*N+1]*s[1]};
				double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

				B[0*N+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*N+0]*B[0*N+0]+bss[0*N+1]*B[1*N+0]);
				B[0*N+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*N+0]*B[0*N+1]+bss[0*N+1]*B[1*N+1]);
				B[1*N+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*N+0]*B[0*N+0]+bss[1*N+1]*B[1*N+0]);
				B[1*N+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*N+0]*B[0*N+1]+bss[1*N+1]*B[1*N+1]);
			}
		}
	}

	delete[]	x;
	delete[]	df;
	delete[]	B;
	delete[]	d;
	delete[]	Nr;

	return 0;
}


int newton()
{
	int N=2;
	double fx=0;
	double *dfx=new double [N];
	double *rfx=new double [N*N];

	double x0=0;
	double x1=0;

	int count=0;
	double E=1;
	double ep=0;
	while(E>ep)
	{
		count++;
		fx=(x0-1)*(x0-1)+10*(x0*x0-x1)*(x0*x0-x1);
		dfx[0]=2*(x0-1)+20*(x0*x0-x1)*2*x0;
		dfx[1]=-20*(x0*x0-x1);

		rfx[0*N+0]=120*x0*x0-40*x1+2;
		rfx[0*N+1]=-40*x0;
		rfx[1*N+0]=-40*x0;
		rfx[1*N+1]=20;

		gauss(rfx,dfx,N);

		x0+=-1*dfx[0];
		x1+=-1*dfx[1];
		E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]);

		cout<<"k="<<count-1<<", E="<<E<<endl;
		cout<<"fx="<<fx<<endl;
		cout<<"x0="<<x0<<", x1="<<x1<<endl;
		cout<<endl;
	}
	cout<<"x0="<<x0<<", x1="<<x1<<endl;

	delete[]	dfx;
	delete[]	rfx;

	return 0;
}


void gauss(double *matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k*N+k];
		
		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i*N+k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i*N+j]-=A*matrix[k*N+j];					
				}
				B[i]-=A*B[k];				
			}
		}
	}
	for(int k=0;k<N;k++) B[k]/=matrix[k*N+k];

}


void PM_AVS(int t_max,int t, double *qn)
{
	//t:タイムステップ　T:総合時間
	double TIME=t/t_max;///時間更新

	if(t==1)
	{
		ofstream fout("particle_movie.mgf");
		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<t_max+1<<endl;//microAVSに出力する総ステップ数。ファイル出力はCON.get_interval()回に1回と最初に行う。		
		fout.close();
	}

	ofstream avs("particle_movie.mgf",ios :: app);
	avs<<"step"<<t<<endl;
	avs<<"sphere"<<endl;
	avs<<"time="<<TIME<<endl;
	avs<<"color"<<endl;

	double red,green,blue;	//粒子の色を表現する3原色

	int num=1;//表示する粒子数

	avs<<num<<endl;

	red=0;
	green=1;
	blue=1;

	avs<<qn[A_X]<<" "<<qn[A_Y]<<" "<<qn[A_Z]<<" ";//座標出力
				
	avs<<1<<" ";//粒子の大きさ出力
	
	avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
		
	avs.close();
		////////////////////
}

void MM_AVS(int t_max, int t,double Dt, double mi, double *qn,double *pn)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	double TIME=t_max*Dt;
	int STEP=t_max+1;		//出力する総ステップ数
	int step=t;
	int n=3;	//出力次元


	if(t==1) 
	{
		ofstream fp("momentum.inp", ios::trunc);			
		fp<<STEP<<endl;//総ステップ数
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//mainファイル書き込み
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<step<<" TIME="<<TIME<<endl;

	fp<<1<<" "<<1<<endl;	//節点数と要素数出力
	
	//節点番号とその座標の出力 
	fp<<0<<" "<<qn[A_X]<<" "<<qn[A_Y]<<" "<<qn[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	fp<<0<<"  0 pt "<<0<<endl;

	fp<<"3 0"<<endl;
	fp<<"3 1 1 1"<<endl;
	fp<<"v_x,"<<endl;
	fp<<"v_y,"<<endl;
	fp<<"v_z,"<<endl;

	//各節点の情報値入力
	fp<<0<<" "<<1/mi*pn[A_X]<<" "<<1/mi*pn[A_Y]<<" "<<1/mi*pn[A_Z]<<endl;
	fp.close();
}
