#pragma once
using namespace std;

#include <iostream>



void gauss(double *matrix,double *B,int N);
int newton();
//int newton2();
int SQP();

//準ニュートン法
int newton2()
{
	int N=2;
	double fx=0;
	double *d=new double [N];
	double *dfx=new double [N];
	double *B=new double[N*N];
	double x0=0, x1=0;


	double x0l=0, x1l=0;

	int count=0;
	double E=1;
	double ep=0;

	while(E>ep)
	{
		count++;
		fx=(x0-1)*(x0-1)+10*(x0*x0-x1)*(x0*x0-x1);

		double dfxl[2]={dfx[0], dfx[1]};

		dfx[0]=2*(x0-1)+20*(x0*x0-x1)*2*x0;
		dfx[1]=-20*(x0*x0-x1);

		E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]);
		if(E<ep)	break;

		d[0]=-1*dfx[0];
		d[1]=-1*dfx[1];
		


		double s[2]={x0-x0l, x1-x1l};
		double y[2]={dfx[0]-dfxl[0], dfx[1]-dfxl[1]};

		double beta=y[0]*s[0]+y[1]*s[1];
		double sigma=s[0]*(B[0*N+0]+B[1*N+0])*s[0]+s[1]*(B[0*N+1]+B[1*N+1])*s[1];

		double bs[2]={B[0*N+0]*s[0]+B[0*N+1]*s[1], B[1*N+0]*s[0]+B[1*N+1]*s[1]};
		double bss[4]={bs[0]*s[0],bs[0]*s[1],bs[1]*s[0],bs[1]*s[1]};

		B[0*N+0]+=1/beta*y[0]*y[0]-1/sigma*(bss[0*N+0]*B[0*N+0]+bss[0*N+1]*B[1*N+0]);
		B[0*N+1]+=1/beta*y[0]*y[1]-1/sigma*(bss[0*N+0]*B[0*N+1]+bss[0*N+1]*B[1*N+1]);
		B[1*N+0]+=1/beta*y[1]*y[0]-1/sigma*(bss[1*N+0]*B[0*N+0]+bss[1*N+1]*B[1*N+0]);
		B[1*N+1]+=1/beta*y[1]*y[1]-1/sigma*(bss[1*N+0]*B[0*N+1]+bss[1*N+1]*B[1*N+1]);

		gauss(B,dfx,N);

		x0+=-1*dfx[0];
		x1+=-1*dfx[1];
		E=sqrt(dfx[0]*dfx[0]+dfx[1]*dfx[1]);

		cout<<"k="<<count-1<<", E="<<E<<endl;
		cout<<"fx="<<fx<<endl;
		cout<<"x0="<<x0<<", x1="<<x1<<endl;
		cout<<endl;
	}
	cout<<"x0="<<x0<<", x1="<<x1<<endl;

	delete[]	d;
	delete[]	dfx;
	delete[]	B;

	return 0;
}

int SQP()
{
	int N=4;
	double fx=0;
	double *dfx=new double [N];
	double *rfx=new double [N*N];

	double x0=10, x1=-10, x2=10, x3=-10;

	double c0=0, c1=0, c2=0;
	double *dc0=new double [N];
	double *dc1=new double [N];
	double *dc2=new double [N];
	double *rc0=new double [N*N];
	double *rc1=new double [N*N];
	double *rc2=new double [N*N];

	double Fp=0;
	double *dFp=new double [N];
	double *rFp=new double [N*N];
	
	double *d=new double [N];

	double p=10;

	double v0=0, v1=0, v2=0;
	
	int Nc=3;
	int Nv=0;
	int *Nv_n=new int [Nc];

	double fc0=0, fc1=0, fc2=0;

	//初期化
	for(int i=0;i<N;i++)
	{
		dfx[i]=0;
		dc0[i]=0;
		dc1[i]=0;
		dc2[i]=0;
		dFp[i]=0;
		d[i]=0;
		for(int j=0;j<N;j++)
		{
			rfx[i*N+j]=0;
			rc0[i*N+j]=0;
			rc1[i*N+j]=0;
			rc2[i*N+j]=0;
			rFp[i*N+j]=0;
		}
	}
	for(int i=0;i<Nc;i++)	Nv_n[i]=0;

	int count=0;
	double E=1;
	double	ep=1e-10;

	while(E>ep)
	{
		count++;
		fx=x0*x0+x1*x1+2*x2*x2+x3*x3-5*x0-5*x1-21*x2+7*x3;

		dfx[0]=2*x0-5;
		dfx[1]=2*x1-5;
		dfx[2]=4*x2-21;
		dfx[3]=2*x3+7;

		rfx[0*N+0]=2;
		rfx[1*N+1]=2;
		rfx[2*N+2]=4;
		rfx[3*N+3]=2;

		c0=x0*x0+x1*x1+x2*x2+x3*x3+x0-x1+x2-x3-8;
		c1=x0*x0+2*x1*x1+x2*x2+2*x3*x3-x0-x3-10;
		c2=2*x0*x0+x1*x1+x2*x2+2*x0-x1-x3-5;

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

		rc0[0*N+0]=2;
		rc0[1*N+1]=2;
		rc0[2*N+2]=2;
		rc0[3*N+3]=2;

		rc1[0*N+0]=2;
		rc1[1*N+1]=4;
		rc1[2*N+2]=2;
		rc1[3*N+3]=4;

		rc2[0*N+0]=4;
		rc2[1*N+1]=2;
		rc2[2*N+2]=2;
		rc2[3*N+3]=0;

		Fp=fx;
		for(int i=0;i<N;i++)
		{
			dFp[i]=dfx[i];
			rFp[i*N+i]=rfx[i*N+i];
		}
		if(c0>0)
		{
			Fp+=p*c0;
			for(int i=0;i<N;i++)
			{
				dFp[i]+=p*dc0[i];
				rFp[i*N+i]+=p*rc0[i*N+i];
			}
		}
		else if(c0==0)
		{
			for(int i=0;i<N;i++)
			{
				if(dc0[i]>0)	dFp[i]+=p*dc0[i];
				if(rc0[i*N+i]>0)	rFp[i*N+i]+=p*rc0[i*N+i];
			}
		}
		if(c1>0)
		{
			Fp+=p*c1;
			for(int i=0;i<N;i++)
			{
				dFp[i]+=p*dc1[i];
				rFp[i*N+i]+=p*rc1[i*N+i];
			}
		}
		else if(c1==0)
		{
			for(int i=0;i<N;i++)
			{
				if(dc1[i]>0)	dFp[i]+=p*dc1[i];
				if(rc1[i*N+i]>0)	rFp[i*N+i]+=p*rc1[i*N+i];
			}
		}
		if(c2>0)
		{
			Fp+=p*c2;
			for(int i=0;i<N;i++)
			{
				dFp[i]+=p*dc2[i];
				rFp[i*N+i]+=p*rc2[i*N+i];
			}
		}
		else if(c2==0)
		{
			for(int i=0;i<N;i++)
			{
				if(dc2[i]>0)	dFp[i]+=p*dc2[i];
				if(rc2[i*N+i]>0)	rFp[i*N+i]+=p*rc2[i*N+i];
			}
		}
		cout<<"Fp="<<Fp<<endl;

		gauss(rFp,dFp,N);

		E=sqrt(dFp[0]*dFp[0]+dFp[1]*dFp[1]+dFp[2]*dFp[2]+dFp[3]*dFp[3]);
		cout<<"E"<<count-1<<"="<<E<<endl;

		if(E<ep)	break;
		cout<<endl;


		x0+=-1*dFp[0];
		x1+=-1*dFp[1];
		x2+=-1*dFp[2];
		x3+=-1*dFp[3];
		cout<<"x0="<<x0<<", x1="<<x1<<", x2="<<x2<<", x3="<<x3<<endl;

		for(int i=0;i<N;i++)	d[i]=-1*dFp[i];
		cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<endl;
		

		fc0=c0+dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3];
		fc1=c1+dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3];
		fc2=c2+dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3];
		
		Nv=0;

		if(fc0<0)	v0=0;
		else
		{
			Nv_n[Nv]=0;
			Nv++;
		}
		if(fc1<0)	v1=0;
		else
		{
			Nv_n[Nv]=1;
			Nv++;
		}
		if(fc2<0)	v2=0;
		else
		{
			Nv_n[Nv]=2;
			Nv++;
		}
		
		if(Nv==1)
		{
			if(Nv_n[Nv-1]==0)	v0=(-dfx[0]-1*(rfx[0*N+0]*d[0]+rfx[0*N+1]*d[1]+rfx[0*N+2]*d[2]+rfx[0*N+3]*d[3]))/2/dc0[0];
			if(Nv_n[Nv-1]==1)	v1=(-dfx[0]-1*(rfx[0*N+0]*d[0]+rfx[0*N+1]*d[1]+rfx[0*N+2]*d[2]+rfx[0*N+3]*d[3]))/2/dc1[0];
			if(Nv_n[Nv-1]==2)	v2=(-dfx[0]-1*(rfx[0*N+0]*d[0]+rfx[0*N+1]*d[1]+rfx[0*N+2]*d[2]+rfx[0*N+3]*d[3]))/2/dc2[0];
		}
		else if(Nv>1)
		{
			double *Nvr=new double[Nv];
			double *Nvl=new double[Nv*Nv];
			for(int i=0;i<Nv;i++)
			{
				Nvr[i]=-dfx[i]-1*(rfx[i*N+0]*d[0]+rfx[i*N+1]*d[1]+rfx[i*N+2]*d[2]+rfx[i*N+3]*d[3]);
				for(int j=0;j<Nv;j++)
				{
					int jv=Nv_n[j];
					if(jv==0)	Nvl[i*Nv+j]=2*dc0[i];
					else if(jv==1)	Nvl[i*Nv+j]=2*dc1[i];
					else if(jv==2)	Nvl[i*Nv+j]=2*dc2[i];
				}
			}
			gauss(Nvl,Nvr,Nv);
			for(int i=0;i<Nv;i++)
			{
				int iv=Nv_n[i];
				if(iv==0)	v0=Nvr[i];
				else if(iv==1)	v1=Nvr[i];
				else if(iv==2)	v2=Nvr[i];
			}
			delete[]	Nvr;
			delete[]	Nvl;
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
	delete[]	dFp;
	delete[]	rFp;
	delete[]	d;
	delete[]	Nv_n;

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

