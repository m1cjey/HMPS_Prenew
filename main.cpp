#pragma once
using namespace std;

#include <iostream>
#include <stdio.h>
#include <math.h>

void gauss(double *matrix,double *B,int N);
int newton();
int newton2();





int main()
{
	int Nx=4;
	int Nv=3;

	double fx=0;
	double *dfx=new double [Nx];
	double *rfx=new double [Nx*Nx];

	double x0=10, x1=-10, x2=10, x3=-10;

	double c0=0, c1=0, c2=0;
	double *dc0=new double [Nx];
	double *dc1=new double [Nx];
	double *dc2=new double [Nx];
	double *rc0=new double [Nx*Nx];
	double *rc1=new double [Nx*Nx];
	double *rc2=new double [Nx*Nx];

	double *d=new double [Nx];
	
	double v0=0, v1=0, v2=0;
	int *Nv_n=new int [Nv];

	double *B=new double [Nx*Nx];
	double *Nr=new double [Nx];
	
	//初期化
	for(int i=0;i<Nx;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dc1[i]=0;
		dc2[i]=0;
		for(int j=0;j<Nx;j++)
		{
			rfx[i*Nx+j]=0;
			rc0[i*Nx+j]=0;
			rc1[i*Nx+j]=0;
			rc2[i*Nx+j]=0;
			if(j==i)	B[i*Nx+j]=1;
			else
			{
				B[i*Nx+j]=0;
			}
		}
		d[i]=0;
		Nr[i]=0;
	}
	Nv_n[0]=0;
	Nv_n[1]=0;
	Nv_n[2]=0;


	int count=0;
	double E=1;
	double	ep=1e-10;

	fx=x0*x0 + x1*x1 + 2*x2*x2 + x3*x3 -5*x0 -5*x1 -21*x2 + 7*x3;

	dfx[0]=2*x0-5;
	dfx[1]=2*x1-5;
	dfx[2]=4*x2-21;
	dfx[3]=2*x3+7;

	rfx[0*Nx+0]=2;
	rfx[1*Nx+1]=2;
	rfx[2*Nx+2]=4;
	rfx[3*Nx+3]=2;

	c0=x0*x0 + x1*x1 + x2*x2 + x3*x3 + x0 -x1 + x2 -x3 -8;
	c1=x0*x0 + 2*x1*x1 + x2*x2 + 2*x3*x3 -x0 -x3 -10;
	c2=2*x0*x0 + x1*x1 + x2*x2 + 2*x0 -x1 -x3 -5;

	dc0[0]=2*x0+1;
	dc0[1]=2*x1-1;
	dc0[2]=2*x2+1;
	dc0[3]=2*x3-1;

	dc1[0]=2*x0-1;
	dc1[1]=4*x1;
	dc1[2]=2*x2;
	dc1[3]=4*x3-1;

	dc2[0]=4*x0+2;
	dc2[1]=2*x1-1;
	dc2[2]=2*x2;
	dc2[3]=-1;

	rc0[0*Nx+0]=2;
	rc0[1*Nx+1]=2;
	rc0[2*Nx+2]=2;
	rc0[3*Nx+3]=2;

	rc1[0*Nx+0]=2;
	rc1[1*Nx+1]=4;
	rc1[2*Nx+2]=2;
	rc1[3*Nx+3]=4;

	rc2[0*Nx+0]=4;
	rc2[1*Nx+1]=2;
	rc2[2*Nx+2]=2;
	rc2[3*Nx+3]=0;

	cout<<"x0="<<x0<<", x1="<<x1<<", x2="<<x2<<", x3="<<x3<<endl;
	cout<<endl;

	while(E>ep)
	{
		count++;

		gauss(rfx,dfx,Nx);

		//d更新
		for(int i=0;i<Nx;i++)	d[i]=-1*dfx[i];
		
		E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]+dfx[2]*dfx[2]+dfx[3]*dfx[3]);
		cout<<"E"<<count<<"="<<E<<endl;
		if(E<=ep)	break;

		//x更新
	/*	x0+=d[0];
		x1+=d[1];
		x2+=d[2];
		x3+=d[3];*/
		x0=-5;//10;//1.25787;
		x1=15;//-10;//2.24825;
		x2=-9;//10;//4.38011;
		x3=3;//-10;//-1.95909;
		
		d[0]=-0.87558-x0;//-15;//0.47116-x0;
		d[1]=8.12596-x1;//25;//0.07174-x1;
		d[2]=11.28458-x2;//-19;//3.18360-x2;
		d[3]=-0.57450-x3;//13;//-0.95163-x3;

		cout<<"x0="<<x0<<", x1="<<x1<<", x2="<<x2<<", x3="<<x3<<endl;
		cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<endl;

		//f更新
		fx=x0*x0 + x1*x1 + 2*x2*x2 + x3*x3 -5*x0 -5*x1 -21*x2 + 7*x3;


		//df更新
		dfx[0]=2*x0-5;
		dfx[1]=2*x1-5;
		dfx[2]=4*x2-21;
		dfx[3]=2*x3+7;
		
		Nr[0]=dfx[0];
		Nr[1]=dfx[1];
		Nr[2]=dfx[2];
		Nr[3]=dfx[3];

		cout<<endl;

		cout<<"alpha(dfx)="<<-1*d[0]/dfx[0]<<", "<<-1*d[1]/dfx[1]<<", "<<-1*d[2]/dfx[2]<<", "<<-1*d[3]/dfx[3]<<endl;
	
		gauss(rfx,Nr,Nx);
	//	cout<<"df="<<-1*Nr[0]<<", "<<-1*Nr[1]<<", "<<-1*Nr[2]<<", "<<-1*Nr[3]<<endl;

		cout<<"alpha(dfx/rfx)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		//c更新
		c0=x0*x0 + x1*x1 + x2*x2 + x3*x3 + x0 -x1 + x2 -x3 -8;
		c1=x0*x0 + 2*x1*x1 + x2*x2 + 2*x3*x3 -x0 -x3 -10;
		c2=2*x0*x0 + x1*x1 + x2*x2 + 2*x0 -x1 -x3 -5;
		
		//dc更新
		dc0[0]=2*x0+1;
		dc0[1]=2*x1-1;
		dc0[2]=2*x2+1;
		dc0[3]=2*x3-1;

		dc1[0]=2*x0-1;
		dc1[1]=4*x1;
		dc1[2]=2*x2;
		dc1[3]=4*x3-1;

		dc2[0]=4*x0+2;
		dc2[1]=2*x1-1;
		dc2[2]=2*x2;
		dc2[3]=-1;

		//Fp計算
		double Fp=fx;
		double *dFp=new double[Nx];		
		double *rFp=new double [Nx*Nx];
		for(int i=0;i<Nx;i++)
		{
			dFp[i]=0;
			for(int j=0;j<Nx;j++)	rFp[i*Nx+j]=0;
		}


		for(int i=0;i<Nx;i++)
		{
			dFp[i]=dfx[i];
			rFp[i*Nx+i]=rfx[i*Nx+i];
		}

		if(c0>0)
		{
			Fp+=c0*10;
			for(int i=0;i<Nx;i++)	dFp[i]+=10*dc0[i];
			for(int i=0;i<Nx;i++)	rFp[i*Nx+i]+=10*rc0[i*Nx+i];
		}
		else
		{
			for(int i=0;i<Nx;i++)	if(dc0[i]>0)	dFp[i]+=10*dc0[i];
			for(int i=0;i<Nx;i++)	if(rc0[i*Nx+i]>0)	rFp[i*Nx+i]+=10*rc0[i*Nx+i];
		}
		if(c1>0)
		{
			Fp+=c1*10;
			for(int i=0;i<Nx;i++)	dFp[i]+=10*dc1[i];
			for(int i=0;i<Nx;i++)	rFp[i*Nx+i]+=10*rc1[i*Nx+i];
		}
		else
		{
			for(int i=0;i<Nx;i++)	if(dc1[i]>0)	dFp[i]+=10*dc1[i];
			for(int i=0;i<Nx;i++)	if(rc1[i*Nx+i]>0)	rFp[i*Nx+i]+=10*rc1[i*Nx+i];
		}
		if(c2>0)
		{
			Fp+=c2*10;
			for(int i=0;i<Nx;i++)	dFp[i]+=10*dc2[i];
			for(int i=0;i<Nx;i++)	rFp[i*Nx+i]+=10*rc2[i*Nx+i];
		}
		else
		{
			for(int i=0;i<Nx;i++)	if(dc2[i]>0)	dFp[i]+=10*dc2[i];
			for(int i=0;i<Nx;i++)	if(rc2[i*Nx+i]>0)	rFp[i*Nx+i]+=10*rc2[i*Nx+i];
		}

		Nr[0]=dFp[0];
		Nr[1]=dFp[1];
		Nr[2]=dFp[2];
		Nr[3]=dFp[3];

		cout<<"alpha(dFp)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		gauss(rFp,Nr,Nx);

//		cout<<"dFp="<<-1*Nr[0]<<", "<<-1*Nr[1]<<", "<<-1*Nr[2]<<", "<<-1*Nr[3]<<endl;

		cout<<"alpha(dFp/rFp)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		Nr[0]=dfx[0];
		Nr[1]=dfx[1];
		Nr[2]=dfx[2];
		Nr[3]=dfx[3];

		gauss(rFp,Nr,Nx);

		cout<<"alpha(dfx/rFp)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		double *rL=new double [Nx*Nx];
		double *dL=new double [Nx];
		for(int i=0;i<Nx;i++)
		{
			dL[i]=0;
			for(int j=0;j<Nx;j++)	rL[i*Nx+j]=0;
		}

		v0=0;
		v1=0;
		v2=0;

		for(int i=0;i<Nx;i++)
		{
			rL[i*Nx+i]=rfx[i*Nx+i]+rc0[i*Nx+i]*v0+rc1[i*Nx+i]*v1+rc2[i*Nx+i]*v2;
			dL[i]=dfx[i]+dc0[i]*v0+dc1[i]*v1+dc2[i]*v2;
		}


		Nr[0]=dL[0];
		Nr[1]=dL[1];
		Nr[2]=dL[2];
		Nr[3]=dL[3];

		cout<<"alpha(dL)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		gauss(rL,Nr,Nx);

		cout<<"alpha(dL/rL)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		Nr[0]=dfx[0];
		Nr[1]=dfx[1];
		Nr[2]=dfx[2];
		Nr[3]=dfx[3];

		gauss(rL,Nr,Nx);

		cout<<"alpha(dfx/rL)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		Nr[0]=dFp[0];
		Nr[1]=dFp[1];
		Nr[2]=dFp[2];
		Nr[3]=dFp[3];

		gauss(rL,Nr,Nx);

		cout<<"alpha(dFp/rL)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		cout<<"Fp="<<Fp<<endl;

		//Bの更新
		double x0_k=10;//-0.87558;
		double x1_k=-10;//8.12596;
		double x2_k=10;//11.28458;
		double x3_k=-10;//-0.57450;

		double x0_k1=-5;//-0.87558;
		double x1_k1=15;//8.12596;
		double x2_k1=-9;//11.28458;
		double x3_k1=3;//-0.57450;
		
		double s[4]={x0_k1-x0_k, x1_k1-x1_k,x2_k1-x2_k,x3_k1-x3_k};

		double dfx_k[4]={2*x0_k-5, 2*x1_k-5, 4*x2_k-21, 2*x3_k+7};
		double dc0_k[4]={2*x0_k+1, 2*x1_k-1, 2*x2_k+1, 2*x3_k-1};
		double dc1_k[4]={2*x0_k-1, 4*x1_k, 2*x2_k, 4*x3_k-1};
		double dc2_k[4]={4*x0_k+2, 2*x1_k-1, 2*x2_k, -1};

		double dfx_k1[4]={2*x0_k1-5, 2*x1_k1-5, 4*x2_k1-21, 2*x3_k1+7};
		double dc0_k1[4]={2*x0_k1+1, 2*x1_k1-1, 2*x2_k1+1, 2*x3_k1-1};
		double dc1_k1[4]={2*x0_k1-1, 4*x1_k1, 2*x2_k1, 4*x3_k1-1};
		double dc2_k1[4]={4*x0_k1+2, 2*x1_k1-1, 2*x2_k1, -1};

		double dL_k1[4]={dfx_k1[0]+dc0_k1[0]*v0+dc1_k1[0]*v1+dc2_k1[0]*v2, dfx_k1[1]+dc0_k1[1]*v0+dc1_k1[1]*v1+dc2_k1[1]*v2, dfx_k1[2]+dc0_k1[2]*v0+dc1_k1[2]*v1+dc2_k1[2]*v2, dfx_k1[3]+dc0_k1[3]*v0+dc1_k1[3]*v1+dc2_k1[3]*v2};
		double dL_k[4]={dfx_k[0]+dc0_k[0]*v0+dc1_k[0]*v1+dc2_k[0]*v2, dfx_k[1]+dc0_k[1]*v0+dc1_k[1]*v1+dc2_k[1]*v2, dfx_k[2]+dc0_k[2]*v0+dc1_k[2]*v1+dc2_k[2]*v2, dfx_k[3]+dc0_k[3]*v0+dc1_k[3]*v1+dc2_k[3]*v2};

		double y[4]={dL_k1[0]-dL_k[0], dL_k1[1]-dL_k[1], dL_k1[2]-dL_k[2], dL_k1[3]-dL_k[3]};

	
		double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];
		double sigma=(s[0]*B[0*Nx+0]+s[1]*B[1*Nx+0]+s[2]*B[2*Nx+0]+s[3]*B[3*Nx+0])*s[0]+(s[0]*B[0*Nx+1]+s[1]*B[1*Nx+1]+s[2]*B[2*Nx+1]+s[3]*B[3*Nx+1])*s[1]+(s[0]*B[0*Nx+2]+s[1]*B[1*Nx+2]+s[2]*B[2*Nx+2]+s[3]*B[3*Nx+2])*s[2]+(s[0]*B[0*Nx+3]+s[1]*B[1*Nx+3]+s[2]*B[2*Nx+3]+s[3]*B[3*Nx+3])*s[3];	

		double B3_f[4]={1/sigma*(B[0*Nx+0]*s[0]+B[0*Nx+1]*s[1]+B[0*Nx+2]*s[2]+B[0*Nx+3]*s[3]),1/sigma*(B[1*Nx+0]*s[0]+B[1*Nx+1]*s[1]+B[1*Nx+2]*s[2]+B[1*Nx+3]*s[3]),1/sigma*(B[2*Nx+0]*s[0]+B[2*Nx+1]*s[1]+B[2*Nx+2]*s[2]+B[2*Nx+3]*s[3]),1/sigma*(B[3*Nx+0]*s[0]+B[3*Nx+1]*s[1]+B[3*Nx+2]*s[2]+B[3*Nx+3]*s[3])};
		double B3_f2[16]={B3_f[0]*s[0],B3_f[0]*s[1],B3_f[0]*s[2],B3_f[0]*s[3],B3_f[1]*s[0],B3_f[1]*s[1],B3_f[1]*s[2],B3_f[1]*s[3],B3_f[2]*s[0],B3_f[2]*s[1],B3_f[2]*s[2],B3_f[2]*s[3],B3_f[3]*s[0],B3_f[3]*s[1],B3_f[3]*s[2],B3_f[3]*s[3]};
		double B3[16]={B3_f2[0*Nx+0]*B[0*Nx+0]+B3_f2[0*Nx+1]*B[1*Nx+0]+B3_f2[0*Nx+2]*B[2*Nx+0]+B3_f2[0*Nx+3]*B[3*Nx+0], B3_f2[0*Nx+0]*B[0*Nx+1]+B3_f2[0*Nx+1]*B[1*Nx+1]+B3_f2[0*Nx+2]*B[2*Nx+1]+B3_f2[0*Nx+3]*B[3*Nx+1], B3_f2[0*Nx+0]*B[0*Nx+2]+B3_f2[0*Nx+1]*B[1*Nx+2]+B3_f2[0*Nx+2]*B[2*Nx+2]+B3_f2[0*Nx+3]*B[3*Nx+2], B3_f2[0*Nx+0]*B[0*Nx+3]+B3_f2[0*Nx+1]*B[1*Nx+3]+B3_f2[0*Nx+2]*B[2*Nx+3]+B3_f2[0*Nx+3]*B[3*Nx+3],
			B3_f2[1*Nx+0]*B[0*Nx+0]+B3_f2[1*Nx+1]*B[1*Nx+0]+B3_f2[1*Nx+2]*B[2*Nx+0]+B3_f2[1*Nx+3]*B[3*Nx+0], B3_f2[1*Nx+0]*B[0*Nx+1]+B3_f2[1*Nx+1]*B[1*Nx+1]+B3_f2[1*Nx+2]*B[2*Nx+1]+B3_f2[1*Nx+3]*B[3*Nx+1], B3_f2[1*Nx+0]*B[0*Nx+2]+B3_f2[1*Nx+1]*B[1*Nx+2]+B3_f2[1*Nx+2]*B[2*Nx+2]+B3_f2[1*Nx+3]*B[3*Nx+2], B3_f2[1*Nx+0]*B[0*Nx+3]+B3_f2[1*Nx+1]*B[1*Nx+3]+B3_f2[1*Nx+2]*B[2*Nx+3]+B3_f2[1*Nx+3]*B[3*Nx+3],
			B3_f2[2*Nx+0]*B[0*Nx+0]+B3_f2[2*Nx+1]*B[1*Nx+0]+B3_f2[2*Nx+2]*B[2*Nx+0]+B3_f2[2*Nx+3]*B[3*Nx+0], B3_f2[2*Nx+0]*B[0*Nx+1]+B3_f2[2*Nx+1]*B[1*Nx+1]+B3_f2[2*Nx+2]*B[2*Nx+1]+B3_f2[2*Nx+3]*B[3*Nx+1], B3_f2[2*Nx+0]*B[0*Nx+2]+B3_f2[2*Nx+1]*B[1*Nx+2]+B3_f2[2*Nx+2]*B[2*Nx+2]+B3_f2[2*Nx+3]*B[3*Nx+2], B3_f2[2*Nx+0]*B[0*Nx+3]+B3_f2[2*Nx+1]*B[1*Nx+3]+B3_f2[2*Nx+2]*B[2*Nx+3]+B3_f2[2*Nx+3]*B[3*Nx+3],
			B3_f2[3*Nx+0]*B[0*Nx+0]+B3_f2[3*Nx+1]*B[1*Nx+0]+B3_f2[3*Nx+2]*B[2*Nx+0]+B3_f2[3*Nx+3]*B[3*Nx+0], B3_f2[3*Nx+0]*B[0*Nx+1]+B3_f2[3*Nx+1]*B[1*Nx+1]+B3_f2[3*Nx+2]*B[2*Nx+1]+B3_f2[3*Nx+3]*B[3*Nx+1], B3_f2[3*Nx+0]*B[0*Nx+2]+B3_f2[3*Nx+1]*B[1*Nx+2]+B3_f2[3*Nx+2]*B[2*Nx+2]+B3_f2[3*Nx+3]*B[3*Nx+2], B3_f2[3*Nx+0]*B[0*Nx+3]+B3_f2[3*Nx+1]*B[1*Nx+3]+B3_f2[3*Nx+2]*B[2*Nx+3]+B3_f2[3*Nx+3]*B[3*Nx+3]};


		B[0*Nx+0]+=1/beta*y[0]*y[0]-B3[0*Nx+0];
		B[0*Nx+1]+=1/beta*y[0]*y[1]-B3[0*Nx+1];
		B[0*Nx+2]+=1/beta*y[0]*y[2]-B3[0*Nx+2];
		B[0*Nx+3]+=1/beta*y[0]*y[3]-B3[0*Nx+3];
		B[1*Nx+0]+=1/beta*y[1]*y[0]-B3[1*Nx+0];
		B[1*Nx+1]+=1/beta*y[1]*y[1]-B3[1*Nx+1];
		B[1*Nx+2]+=1/beta*y[1]*y[2]-B3[1*Nx+2];
		B[1*Nx+3]+=1/beta*y[1]*y[3]-B3[1*Nx+3];
		B[2*Nx+0]+=1/beta*y[2]*y[0]-B3[2*Nx+0];
		B[2*Nx+1]+=1/beta*y[2]*y[1]-B3[2*Nx+1];
		B[2*Nx+2]+=1/beta*y[2]*y[2]-B3[2*Nx+2];
		B[2*Nx+3]+=1/beta*y[2]*y[3]-B3[2*Nx+3];
		B[3*Nx+0]+=1/beta*y[3]*y[0]-B3[3*Nx+0];
		B[3*Nx+1]+=1/beta*y[3]*y[1]-B3[3*Nx+1];
		B[3*Nx+2]+=1/beta*y[3]*y[2]-B3[3*Nx+2];
		B[3*Nx+3]+=1/beta*y[3]*y[3]-B3[3*Nx+3];

		Nr[0]=dfx[0];
		Nr[1]=dfx[1];
		Nr[2]=dfx[2];
		Nr[3]=dfx[3];

		gauss(B,Nr,Nx);

		cout<<"alpha(dfx/B)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;

		Nr[0]=dL[0];
		Nr[1]=dL[1];
		Nr[2]=dL[2];
		Nr[3]=dL[3];
		gauss(B,Nr,Nx);
		cout<<"alpha(dL/B)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;
		
		Nr[0]=dFp[0];
		Nr[1]=dFp[1];
		Nr[2]=dFp[2];
		Nr[3]=dFp[3];
		gauss(B,Nr,Nx);
		cout<<"alpha(dFp/B)="<<-1*d[0]/Nr[0]<<", "<<-1*d[1]/Nr[1]<<", "<<-1*d[2]/Nr[2]<<", "<<-1*d[3]/Nr[3]<<endl;







		delete[]	dFp;
		delete[]	rFp;
		delete[]	rL;
		delete[]	dL;

		int nv=0;
		double fc0=c0+(dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3]);
		double fc1=c1+(dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3]);
		double fc2=c2+(dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3]);
		cout<<"fc0="<<fc0<<", fc1="<<fc1<<", fc2="<<fc2<<endl;
		if(fc0!=0)	v0=0;
		else
		{
			Nv_n[nv]=0;
			nv++;
		}
		if(fc1!=0)	v1=0;
		else
		{
			Nv_n[nv]=1;
			nv++;
		}
		if(fc2!=0)	v2=0;
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
					if(jv==0)	Nl[i*nv+0]=rc0[i]*d[i]+dc0[i];
					else if(jv==1)	Nl[i*nv+1]=rc1[i]*d[i]+dc1[i];
					else if(jv==2)	Nl[i*nv+2]=rc2[i]*d[i]+dc2[i];
				}
				int iv=Nv_n[i];
				Nr[iv]=-1*(dfx[iv]+rfx[iv*Nx+iv]*d[iv]);
			}
			gauss(Nl,Nr,nv);
			for(int i=0;i<nv;i++)
			{
				int iv=Nv_n[i];
				if(iv==0)	v0=Nr[i];
				else if(iv==1)	v1=Nr[i];
				else if(iv==2)	v2=Nr[i];
			}
		}
		else if(nv==1)
		{
			int iv=Nv_n[nv-1];
			if(iv==0)	v0=-1/(dc0[0]+rc0[0]*d[0])*(dfx[0]+rfx[0]*d[0]);
			else if(iv==1)	v1=-1/(dc1[0]+rc1[0]*d[0])*(dfx[0]+rfx[0]*d[0]);
			else if(iv==2)	v2=-1/(dc2[0]+rc2[0]*d[0])*(dfx[0]+rfx[0]*d[0]);
		}
		cout<<"v0="<<v0<<", v1="<<v1<<", v2="<<v2<<endl;

	}




	delete[]	dfx;
	delete[]	rfx;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	rc0;
	delete[]	rc1;
	delete[]	rc2;
	delete[]	d;
	delete[]	Nv_n;
	delete[]	Nr;
	delete[]	B;
	return 0;
}


/*
//準ニュートン法
int newton2()
{
	int N=2;

	//x_kとB_kの初期設定
	double x0_k=0, x1_k=1;

	double *B_k=new double[N*N];
	B_k[0*N+0]=1;
	B_k[0*N+1]=0;
	B_k[1*N+0]=0;
	B_k[1*N+1]=1;


	//df_kの初期設定
	double *df_k=new double [N];
	df_k[0]=2*(x0_k-1)+20*(x0_k*x0_k-x1_k)*2*x0_k;
	df_k[1]=-20*(x0_k*x0_k-x1_k);

	double E=1;
	double ep=0;

	
	////df_k=0なら計算終了
	E=sqrt(df_k[0]*df_k[0]+df_k[1]*df_k[1]);
	if(E<ep)	return 0;



	////df_k>0なら反復計算
	int count=0;

	//dの初期設定
	double *d_k=new double [N];	
	d_k[0]=0;
	d_k[1]=0;

	while(E>ep)
	{
		count++;
		
		//dの更新
		gauss(B_k,df_k,N);
		d_k[0]=-1*df_k[0];
		d_k[1]=-1*df_k[1];

		//a_kの更新
		double Es=1;
		while(Es>ep)
		{
			double s=0.5;
			double f_k=(x0_k-1)*(x0_k-1)+10*(x0_k*x0_k-x1_k)*(x0_k*x0_k-x1_k);
			double f_kd=(x0_k+d_k[0]-1)*(x0_k+d_k[0]-1)+10*(x0_k+d_k[0]*x0_k+d_k[0]-x1_k+d_k[1])*(x0_k+d_k[0]*x0_k+d_k[0]-x1_k+d_k[1]);
			double f_k_s=(1-s)*f_k+f_kd;

			double f_sk=((x0_k+d_k[0]*s)-1)*((x0_k+d_k[0]*s)-1)+10*((x0_k+d_k[0]*s)*(x0_k+d_k[0]*s)-x1_k+d_k[1]*s)*((x0_k+d_k[0]*s)*(x0_k+d_k[0]*s)-(x1_k+d_k[1]*s));
			Es=f_k_s-f_sk;
			if(Es<ep)	break;
			d_k[0]*=s;
			d_k[1]*=s;
		}

		//x_kの更新
		x0_k+=d_k[0];
		x1_k+=d_k[1];

		//B_kの更新

		//d_f_kの更新

		double dfxl[2]={df_k[0], df_k[1]};




		double s[2]={0, 0};
		double y[2]={0,0};

		double beta=y[0]*s[0]+y[1]*s[1];
		cout<<"beta"<<beta<<endl;
		if(beta>0)
		{
			double sigma=(s[0]*B_k[0]+s[1]*B_k[2])*s[0]+(s[0]*B_k[2]+s[1]*B_k[3])*s[1];

			double bs[2]={B_k[0*N+0]*s[0]+B_k[0*N+1]*s[1], B_k[1*N+0]*s[0]+B_k[1*N+1]*s[1]};
			double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

			B_k[0*N+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*N+0]*B_k[0*N+0]+bss[0*N+1]*B_k[1*N+0]);
			B_k[0*N+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*N+0]*B_k[0*N+1]+bss[0*N+1]*B_k[1*N+1]);
			B_k[1*N+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*N+0]*B_k[0*N+0]+bss[1*N+1]*B_k[1*N+0]);
			B_k[1*N+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*N+0]*B_k[0*N+1]+bss[1*N+1]*B_k[1*N+1]);
		}

		/*

		int Na=100;
		double fxa_min=fx;
		double alpha_min=alpha;
		for(int i=0;i<Na;i++)
		{
			double a=(i+1)*alpha/Na;
			double x0a=x0+a*d[0];
			double x1a=x0+a*d[1];
			double fxa=(x0a-1)*(x0a-1)+10*(x0a*x0a-x1a)*(x0a*x0a-x1a);
			if(fxa<fxa_min)
			{
				fxa_min=fxa;
				alpha_min=a;
			}
		}//

		x0l=x0;
		x1l=x1;

		x0+=alpha_min*d[0];
		x1+=alpha_min*d[1];
		cout<<"alpha="<<alpha_min<<endl;
		cout<<"x0="<<x0<<", x1="<<x1<<endl;

		cout<<endl;
		*/
	/*}

	delete[]	d_k;
	delete[]	df_k;
	delete[]	B_k;

	return 0;
}*/
/*
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
*/

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



