#pragma once
using namespace std;

#include <iostream>
#include <stdio.h>
#include <math.h>

void gauss(double *matrix,double *B,int N);
int newton();
int q_newton();

int main()
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
		
	double Fp=0;
	double p=10.0;
	double *dFp=new double[Nx];		

	double *v=new double [Nv];
	int *Nv_n=new int [Nv];


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
	}

	fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
	dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

	c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
	c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
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
	int nv=0;
	double fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
	double fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
	double fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
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
		double *Nr=new double [nv];
		double *Nl=new double [nv*nv];
		for(int i=0;i<nv;i++)
		{
			for(int j=0;j<nv;j++)
			{
				int jv=Nv_n[j];
				if(jv==0)	Nl[i*nv+j]=dc0[j];
				else if(jv==1)	Nl[i*nv+j]=dc1[j];
				else if(jv==2)	Nl[i*nv+j]=dc2[j];
			}
			Nr[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
		}
		gauss(Nl,Nr,nv);
		for(int i=0;i<nv;i++)
		{
			int iv=Nv_n[i];
			if(iv==0)	v[0]=Nr[i];
			else if(iv==1)	v[1]=Nr[i];
			else if(iv==2)	v[2]=Nr[i];
		}
	}
	else if(nv==1)
	{
		int iv=Nv_n[nv-1];
		if(iv==0)	v[0]=-1/(dc0[0]+dc0[1]+dc0[2]+dc0[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==1)	v[1]=-1/(dc1[0]+dc1[1]+dc1[2]+dc1[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
		else if(iv==2)	v[2]=-1/(dc2[0]+dc2[1]+dc2[2]+dc2[3])*(dfx[0]+dfx[1]+dfx[2]+dfx[3] + B[0*Nx+0]*d[0]+B[0*Nx+1]*d[1]+B[0*Nx+2]*d[2]+B[0*Nx+3]*d[3]);
	}

	double E=1;
	double ep=1e-10;

	//E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);
	E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]+dfx[2]*dfx[2]+dfx[3]*dfx[3]);

	//出力
	cout<<"x0="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
	cout<<"c0="<<c0<<", "<<c1<<", "<<c2<<endl;
	cout<<"fc0="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
	cout<<"v0="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
	cout<<"Fp0="<<Fp<<endl;
	//cout<<"dFp0="<<dFp[0]<<", "<<dFp[1]<<", "<<dFp[2]<<", "<<dFp[3]<<endl;
	cout<<"dfx0="<<dfx[0]<<", "<<dfx[1]<<", "<<dfx[2]<<", "<<dfx[3]<<endl;
	//cout<<"fx0="<<fx<<endl;
	cout<<"E0="<<E<<endl<<endl;

	if(E<ep)	return 0;
	else
	{
		int count=0;
		while(E>ep)
		{
			count++;
			double x_k[4]={x[0], x[1],x[2],x[3]};
			double dFp_k[4]={dFp[0],dFp[1],dFp[2],dFp[3]};
			double dL_k[4]={dfx[0]+v[0]*dc0[0]+v[1]*dc1[0]+v[2]*dc2[0],
				dfx[1]+v[0]*dc0[1]+v[1]*dc1[1]+v[2]*dc2[1],
				dfx[2]+v[0]*dc0[2]+v[1]*dc1[2]+v[2]*dc2[2],
				dfx[3]+v[0]*dc0[3]+v[1]*dc1[3]+v[2]*dc2[3]};
			
			Nr[0]=dfx[0];	Nr[1]=dfx[1];	Nr[2]=dfx[2];	Nr[3]=dfx[3];	
			//Nr[0]=dFp[0];	Nr[1]=dFp[1];	Nr[2]=dFp[2];	Nr[3]=dFp[3];	
			gauss(B,Nr,Nx);
			d[0]=-Nr[0];	d[1]=-Nr[1];	d[2]=-Nr[2];	d[3]=-Nr[3];	

			double Fp_min=Fp;
			double a_min=0;
			for(int i=0;i<100000;i++)
			{
				double alpha=(i+1)*1e-5;
				
				double x_a[4]={x[0]+d[0]*alpha, x[1]+d[1]*alpha, x[2]+d[2]*alpha, x[3]+d[3]*alpha};
			
				double c0_a=x_a[0]*x_a[0] +	x_a[1]*x_a[1] + x_a[2]*x_a[2] + x_a[3]*x_a[3] + x_a[0] -x_a[1] + x_a[2] -x_a[3] -8;
				double c1_a=x_a[0]*x_a[0] +	2*x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[3]*x_a[3] -x_a[0] -x_a[3] -10;
				double c2_a=2*x_a[0]*x_a[0] + x_a[1]*x_a[1] + x_a[2]*x_a[2] + 2*x_a[0] -x_a[1] -x_a[3] -5;	
				
				double fx_a=x_a[0]*x_a[0] + x_a[1]*x_a[1] + 2*x_a[2]*x_a[2] + x_a[3]*x_a[3] -5*x_a[0] -5*x_a[1] -21*x_a[2] + 7*x_a[3];

				double Fp_a=fx_a;
				if(c0_a>0)	Fp_a+=c0_a*p;
				if(c1_a>0)	Fp_a+=c1_a*p;
				if(c2_a>0)	Fp_a+=c2_a*p;

				if(Fp_a<Fp_min)
				{
					Fp_min=Fp_a;
					a_min=alpha;
				}
			}

			for(int i=0;i<Nx;i++)	x[i]+=d[i]*a_min;

			fx=x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] -5*x[0] -5*x[1] -21*x[2] + 7*x[3];
			dfx[0]=2*x[0]-5;	dfx[1]=2*x[1]-5;	dfx[2]=4*x[2]-21;	dfx[3]=2*x[3]+7;

			c0=x[0]*x[0] +	x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] -x[1] + x[2] -x[3] -8;
			c1=x[0]*x[0] +	2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] -x[0] -x[3] -10;
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
			int nv=0;
			double fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
			double fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
			double fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
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
				double *Nr=new double [nv];
				double *Nl=new double [nv*nv];
				for(int i=0;i<nv;i++)
				{
					for(int j=0;j<nv;j++)
					{
						int jv=Nv_n[j];
						if(jv==0)	Nl[i*nv+j]=dc0[j];
						else if(jv==1)	Nl[i*nv+j]=dc1[j];
						else if(jv==2)	Nl[i*nv+j]=dc2[j];
					}
					Nr[i]=-1*(dfx[i]+B[i*Nx+0]*d[0]+B[i*Nx+1]*d[1]+B[i*Nx+2]*d[2]+B[i*Nx+3]*d[3]);
				}
				gauss(Nl,Nr,nv);
				for(int i=0;i<nv;i++)
				{
					int iv=Nv_n[i];
					if(iv==0)	v[0]=Nr[i];
					else if(iv==1)	v[1]=Nr[i];
					else if(iv==2)	v[2]=Nr[i];
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
			E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]+dfx[2]*dfx[2]+dfx[3]*dfx[3]);
			
			//出力
			cout<<"x"<<count<<"="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
			cout<<"fc"<<count<<"="<<fc0<<", "<<fc1<<", "<<fc2<<endl;
			cout<<"v"<<count<<"="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
			cout<<"Fp"<<count<<"="<<Fp_min<<", alpha="<<a_min<<endl;
			//cout<<"fx"<<count<<"="<<f<<endl;
			cout<<"d"<<count<<"="<<d[0] <<", "<<d[1] <<", "<<d[2] <<", "<<d[3] <<endl;
			cout<<"da"<<count<<"="<<d[0]*a_min<<", "<<d[1]*a_min<<", "<<d[2]*a_min<<", "<<d[3]*a_min<<endl;
			cout<<"E"<<count<<"="<<E<<endl;
			cout<<endl;
			if(E<ep)	break;

			for(int i=0;i<Nx;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];

			double s[4]={x[0]-x_k[0],x[1]-x_k[1], x[2]-x_k[2], x[3]-x_k[3]};
			//double y[4]={dFp[0]-dFp_k[0], dFp[1]-dFp_k[1], dFp[2]-dFp_k[2], dFp[3]-dFp_k[3]};
			double y[4]={dL[0]-dL_k[0], dL[1]-dL_k[1], dL[2]-dL_k[2], dL[3]-dL_k[3]};

			double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];

			if(beta>0)
			{

				double sigma=(s[0]*B[0*Nx+0]+s[1]*B[1*Nx+0]+s[2]*B[2*Nx+0]+s[3]*B[3*Nx+0])*s[0]
				+(s[0]*B[0*Nx+1]+s[1]*B[1*Nx+1]+s[2]*B[2*Nx+1]+s[3]*B[3*Nx+1])*s[1]
				+(s[0]*B[0*Nx+2]+s[1]*B[1*Nx+2]+s[2]*B[2*Nx+2]+s[3]*B[3*Nx+2])*s[2]
				+(s[0]*B[0*Nx+3]+s[1]*B[1*Nx+3]+s[2]*B[2*Nx+3]+s[3]*B[3*Nx+3])*s[3];
				
				double seta=1.0;
				if(beta<0.2*sigma)
				{
					seta=0.8*sigma/(sigma-beta);
					
					y[0]*=seta;
					y[1]*=seta;
					y[2]*=seta;
					y[3]*=seta;
					y[0]+=(1-seta)*(B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3]);
					y[1]+=(1-seta)*(B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3]);
					y[2]+=(1-seta)*(B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3]);
					y[3]+=(1-seta)*(B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3]);
				
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
			//Bに関する出力
			cout<<"beta"<<count<<"="<<beta<<endl;
			//*/
		}
	}

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

			for(int i=0;i<100000;i++)
			{
				double alpha=(i+1)*1e-5;
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



