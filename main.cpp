#pragma once
using namespace std;

#include <iostream>



void gauss(double *matrix,double *B,int N);

int main()
{
	int N=4;
	int Nc=3;

	double p=10;

	double fx=0;
	double *dfx=new double [N];
	double *rfx=new double [N*N];

	double L=0;
	double *dL=new double [N];
	double *rL=new double [N*N];

	
	double *v=new double [Nc];

	double c0=0;
	double c1=0;
	double c2=0;

	double fc0=0;
	double fc1=0;
	double fc2=0;


	double *dc0=new double [N];
	double *dc1=new double [N];
	double *dc2=new double [N];

	double *rc0=new double [N*N];
	double *rc1=new double [N*N];
	double *rc2=new double [N*N];

	double *x=new double [N];
	double *d=new double [N];

	double *Nr=new double [N];
	double *Nl=new double [N*N];

	double *Nvr=new double [Nc];
	double *Nvl=new double [Nc*Nc];

	double Fp=0;
	double last_Fp=0;
	double *dFp=new double [N];
	double *rFp=new double [N*N];

	double *B=new double [N*N];
	double *s=new double [N];
	double *y=new double [N];

	for(int i=0;i<N;i++)
	{
		dfx[i]=0;
		dL[i]=0;
		Nr[i]=0;
		for(int j=0;j<N;j++)
		{
			rfx[i*N+j]=0;
			rL[i*N+j]=0;
			rc0[i*N+j]=0;
			rc1[i*N+j]=0;
			rc2[i*N+j]=0;
			Nl[i*N+j]=0;
			rFp[i*N+j]=0;
			B[i*N+j]=0;
		}
		dc0[i]=0;
		dc1[i]=0;
		dc2[i]=0;

		x[i]=0;
		d[i]=0;
		dFp[i]=0;
		s[i]=0;
		y[i]=0;
	}

	for(int i=0;i<Nc;i++)
	{
		v[i]=0;
		Nvr[i]=0;
		for(int j=0;j<Nc;j++)	Nvl[i*Nc+j]=0;
	}


	x[0]=-5;
	x[1]=15;
	x[2]=-9;
	x[3]=3;
	cout<<"x0="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;

//	double alpha=0.1;

	//for(int i=0;i<N;i++)	d[i]=x0[i]-x[i];

	double ep=1.0e-12;	
	double E=1;

	fx=x[0]*x[0]+x[1]*x[1]+2*x[2]*x[2]+x[3]*x[3]-5*x[0]-5*x[1]-21*x[2]+7*x[3];

	c0=x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[0]-x[1]+x[2]-x[3]-8;
	c1=x[0]*x[0]+2*x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]-x[0]-x[3]-10;
	c2=2*x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+2*x[0]-x[1]-x[3]-5;
	cout<<"c0="<<c0<<", c1="<<c1<<", c2="<<c2<<endl;


	dfx[0]=2*x[0]-5;
	dfx[1]=2*x[1]-5;
	dfx[2]=2*2*x[2]-21;
	dfx[3]=2*x[3]+7;

	rfx[0*N+0]=2;
	rfx[1*N+1]=2;
	rfx[2*N+2]=4;
	rfx[3*N+3]=2;

	dc0[0]=2*x[0]+1;
	dc0[1]=2*x[1]-1;
	dc0[2]=2*x[2]+1;
	dc0[3]=2*x[3]-1;

	dc1[0]=2*x[0]-1;
	dc1[1]=2*2*x[1];
	dc1[2]=2*x[2];
	dc1[3]=2*2*x[3]-1;

	dc2[0]=2*2*x[0]+2;
	dc2[1]=2*x[1]-1;
	dc2[2]=2*x[2];
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
		
	//Fp‚ð‹‚ß‚é
	Fp=fx;
	if(c0>0)	Fp+=p*c0;
	if(c1>0)	Fp+=p*c1;
	if(c2>0)	Fp+=p*c2;
	cout<<"Fp="<<Fp<<endl;


	for(int i=0;i<N;i++)
	{
		B[i*N+i]=1;
		s[i]=x[i]-last_x[i];
		y[i]=dL[i]-last_dL[i];
	}

	for(int i=0;i<N;i++)
	{
		s[i]=x[i]-last_x[i];
		y[i]=dL[i]-last_dL[i];
	}

	double beta=y[0]*s[0]+y[1]*s[1]+y[2]*s[2]+y[3]*s[3];
	double sigma=s[0]*(B[0*N+0]+B[1*N+0]+B[2*N+0]+B[3*N+0])*s[0]+
		s[1]*(B[0*N+1]+B[1*N+1]+B[2*N+1]+B[3*N+1])*s[1]+
		s[2]*(B[0*N+2]+B[1*N+2]+B[2*N+2]+B[3*N+2])*s[2]+
		s[3]*(B[0*N+3]+B[1*N+3]+B[2*N+3]+B[3*N+3])*s[3];

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)	B[i*N+j]+=1/beta*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3]);
	}


	//dk‚ð‹‚ß‚é
	for(int i=0;i<N;i++)
	{
		dFp[i]=dfx[i];
		if(c0>0 || (c0==0 && dc0[i]>0) )	dFp[i]+=p*dc0[i];
		if(c1>0 || (c1==0 && dc1[i]>0) )	dFp[i]+=p*dc1[i];
		if(c2>0 || (c2==0 && dc2[i]>0) )	dFp[i]+=p*dc2[i];
	}
	cout<<"dc0={"<<dc0[0]<<","<<dc0[1]<<","<<dc0[2]<<","<<dc0[3]<<"}, dc1={"<<dc1[0]<<","<<dc1[1]<<","<<dc1[2]<<","<<dc1[3]<<"}, dc2={"<<dc2[0]<<","<<dc2[1]<<","<<dc2[2]<<","<<dc2[3]<<"}"<<endl;

	for(int i=0;i<N;i++)
	{
		rFp[i*N+i]=rfx[i*N+i];
		if(c0>0 || (c0==0 && rc0[i*N+i]>0) )	rFp[i*N+i]+=p*rc0[i*N+i];
		if(c1>0 || (c1==0 && rc1[i*N+i]>0) )	rFp[i*N+i]+=p*rc1[i*N+i];
		if(c2>0 || (c2==0 && rc2[i*N+i]>0) )	rFp[i*N+i]+=p*rc2[i*N+i];
	}

	cout<<"rFp="<<"{"<<rFp[0*N+0]<<","<<rFp[0*N+1]<<","<<rFp[0*N+2]<<","<<rFp[0*N+3]<<"}, {"<<rFp[1*N+0]<<","<<rFp[1*N+1]<<","<<rFp[1*N+2]<<","<<rFp[1*N+3]<<"}, {"<<rFp[2*N+0]<<","<<rFp[2*N+1]<<","<<rFp[2*N+2]<<","<<rFp[2*N+3]<<"},  {"<<rFp[3*N+0]<<","<<rFp[3*N+1]<<","<<rFp[3*N+2]<<","<<rFp[3*N+3]<<"}"<<endl;
	cout<<"dFp="<<dFp[0]<<","<<dFp[1]<<","<<dFp[2]<<","<<dFp[3]<<endl;
	gauss(rFp,dFp,N);
	for(int i=0;i<N;i++)	d[i]=-1*dfx[i];		
	cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<endl;

	while(E>0)
	{			
		cout<<endl;
	
		//xk+1‚ð‹‚ß‚é
		for(int i=0;i<N;i++)	x[i]+=d[i];
		cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;

		fx=x[0]*x[0]+x[1]*x[1]+2*x[2]*x[2]+x[3]*x[3]-5*x[0]-5*x[1]-21*x[2]+7*x[3];
	
		c0=x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[0]-x[1]+x[2]-x[3]-8;
		c1=x[0]*x[0]+2*x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]-x[0]-x[3]-10;
		c2=2*x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+2*x[0]-x[1]-x[3]-5;

		dfx[0]=2*x[0]-5;
		dfx[1]=2*x[1]-5;
		dfx[2]=2*2*x[2]-21;
		dfx[3]=2*x[3]+7;

		rfx[0*N+0]=2;
		rfx[1*N+1]=2;
		rfx[2*N+2]=4;
		rfx[3*N+3]=2;

		dc0[0]=2*x[0]+1;
		dc0[1]=2*x[1]-1;
		dc0[2]=2*x[2]+1;
		dc0[3]=2*x[3]-1;

		dc1[0]=2*x[0]-1;
		dc1[1]=2*2*x[1];
		dc1[2]=2*x[2];
		dc1[3]=2*2*x[3]-1;

		dc2[0]=2*2*x[0]+2;
		dc2[1]=2*x[1]-1;
		dc2[2]=2*x[2];
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

		for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	rL[i*N+j]=rfx[i*N+j]+v[0]*rc0[i*N+j]+v[1]*rc1[i*N+j]+v[2]*rc2[i*N+j];

		//Fp‚ð‹‚ß‚é
		last_Fp=Fp;
		Fp=fx;
		if(c0>0)	Fp+=p*c0;
		if(c1>0)	Fp+=p*c1;
		if(c2>0)	Fp+=p*c2;
		cout<<"Fp="<<Fp<<endl;

		//vk+1‚ð‹‚ß‚é
		for(int i=0;i<Nc;i++)
		{
			Nvr[i]=-1*dfx[i]-1*(rfx[i*N+0]*d[0]+rfx[i*N+1]*d[1]+rfx[i*N+2]*d[2]+rfx[i*N+3]*d[3]);
			Nvl[i*Nc+0]=2*dc0[i];
			Nvl[i*Nc+1]=2*dc1[i];
			Nvl[i*Nc+2]=2*dc2[i];
		}
		gauss(Nvl,Nvr,Nc);
		for(int i=0;i<Nc;i++)	v[i]=Nvr[i];
		fc0=c0+dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3];
		fc1=c1+dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3];
		fc2=c2+dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3];
		if(fc0!=0)	v[0]=0;
		if(fc1!=0)	v[1]=0;
		if(fc2!=0)	v[2]=0;
		cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;

		//dk‚ð‹‚ß‚é
		for(int i=0;i<N;i++)
		{
			dFp[i]=dfx[i];
			if(c0>0 || (c0==0&&dc0[i]>0) )	dFp[i]+=p*dc0[i];
			if(c1>0 || (c1==0&&dc1[i]>0) )	dFp[i]+=p*dc1[i];
			if(c2>0 || (c2==0&&dc2[i]>0) )	dFp[i]+=p*dc2[i];
		}
		for(int i=0;i<N;i++)
		{
			rFp[i*N+i]=rfx[i*N+i];
			if(c0>0 || (c0==0 && rc0[i*N+i]>0) )	rFp[i*N+i]+=p*rc0[i*N+i];
			if(c1>0 || (c1==0 && rc1[i*N+i]>0) )	rFp[i*N+i]+=p*rc1[i*N+i];
			if(c2>0 || (c2==0 && rc2[i*N+i]>0) )	rFp[i*N+i]+=p*rc2[i*N+i];
		}
		cout<<"rFp="<<"{"<<rFp[0*N+0]<<","<<rFp[0*N+1]<<","<<rFp[0*N+2]<<","<<rFp[0*N+3]<<"}, {"<<rFp[1*N+0]<<","<<rFp[1*N+1]<<","<<rFp[1*N+2]<<","<<rFp[1*N+3]<<"}, {"<<rFp[2*N+0]<<","<<rFp[2*N+1]<<","<<rFp[2*N+2]<<","<<rFp[2*N+3]<<"},  {"<<rFp[3*N+0]<<","<<rFp[3*N+1]<<","<<rFp[3*N+2]<<","<<rFp[3*N+3]<<"}"<<endl;
		cout<<"dFp="<<dFp[0]<<","<<dFp[1]<<","<<dFp[2]<<","<<dFp[3]<<endl;
		gauss(rFp,dFp,N);
		for(int i=0;i<N;i++)	d[i]=-1*dfx[i];		
		cout<<"d="<<d[0]<<", "<<d[1]<<", "<<d[2]<<", "<<d[3]<<endl;

		//E‚ð‹‚ß‚é
		E=fabs(Fp-last_Fp);
		cout<<"E="<<E<<endl;
		/*
		fx=x[0]*x[0]+x[1]*x[1]+2*x[2]*x[2]+x[3]*x[3]-5*x[0]-5*x[1]-21*x[2]+7*x[3];
	
		dfx[0]=2*x[0]-5;
		dfx[1]=2*x[1]-5;
		dfx[2]=2*2*x[2]-21;
		dfx[3]=2*x[3]+7;

		rfx[0*N+0]=2;
		rfx[1*N+1]=2;
		rfx[2*N+2]=4;
		rfx[3*N+3]=2;

		c0=x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[0]-x[1]+x[2]-x[3]-8;
		c1=x[0]*x[0]+2*x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]-x[0]-x[3]-10;
		c2=2*x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+2*x[0]-x[1]-x[3]-5;

		dc0[0]=2*x[0]+1;
		dc0[1]=2*x[1]-1;
		dc0[2]=2*x[2]+1;
		dc0[3]=2*x[3]-1;

		dc1[0]=2*x[0]-1;
		dc1[1]=2*2*x[1];
		dc1[2]=2*x[2];
		dc1[3]=2*2*x[3]-1;

		dc2[0]=2*2*x[0]+2;
		dc2[1]=2*x[1]-1;
		dc2[2]=2*x[2];
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

		for(int i=0;i<N;i++)	dL[i]=dfx[i]+v[0]*dc0[i]+v[1]*dc1[i]+v[2]*dc2[i];
		for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	rL[i*N+j]=rfx[i*N+j]+v[0]*rc0[i*N+j]+v[1]*rc1[i*N+j]+v[2]*rc2[i*N+j];

		double last_Fp=Fp;

		Fp=fx;
		if(c0>0)	Fp+=p*c0;
		if(c1>0)	Fp+=p*c1;
		if(c2>0)	Fp+=p*c2;
		cout<<"Fp="<<Fp<<endl;

		dFp=-1*(d[0]*rL[0*N+0]*d[0]+d[1]*rL[1*N+1]*d[1]+d[2]*rL[2*N+2]*d[2]+d[3]*rL[3*N+3]*d[3])+v[0]*c0+v[1]*c1+v[2]*c2;
		if(c0>0)	dFp-=p*c0;
		if(c1>0)	dFp-=p*c1;
		if(c2>0)	dFp-=p*c2;


		for(int i=0;i<Nc;i++)
		{
			Nr[i]=-1*dfx[i]-1*(rL[i*N+0]*d[0]+rL[i*N+1]*d[1]+rL[i*N+2]*d[2]+rL[i*N+3]*d[3]);
			Nl[i*Nc+0]=dc0[i];
			Nl[i*Nc+1]=dc1[i];
			Nl[i*Nc+2]=dc2[i];
		}

		gauss(Nl,Nr,Nc);

		for(int i=0;i<Nc;i++)
		{
			v[i]=Nr[i];
		}
		
		double fc0=c0+dc0[0]*d[0]+dc0[1]*d[1]+dc0[2]*d[2]+dc0[3]*d[3];
		double fc1=c1+dc1[0]*d[0]+dc1[1]*d[1]+dc1[2]*d[2]+dc1[3]*d[3];
		double fc2=c2+dc2[0]*d[0]+dc2[1]*d[1]+dc2[2]*d[2]+dc2[3]*d[3];

		if(fc0!=0)	v[0]=0;
		if(fc1!=0)	v[1]=0;
		if(fc2!=0)	v[2]=0;

		cout<<"x="<<x[0]<<", "<<x[1]<<", "<<x[2]<<", "<<x[3]<<endl;
		cout<<"v="<<v[0]<<", "<<v[1]<<", "<<v[2]<<endl;
		E=fabs(last_Fp-Fp);
		cout<<"E="<<E<<endl;

		for(int i=0;i<N;i++)	x[i]+=alpha*d[i];
		for(int i=0;i<N;i++)	d[i]=x0[i]-x[i];*/
	}

	delete[]	dfx;
	delete[]	rfx;
	delete[]	dL;
	delete[]	rL;
	delete[]	v;
	delete[]	dc0;
	delete[]	dc1;
	delete[]	dc2;
	delete[]	rc0;
	delete[]	rc1;
	delete[]	rc2;
	delete[]	x;
	delete[]	d;
	delete[]	Nr;
	delete[]	Nl;
	delete[]	Nvr;
	delete[]	Nvl;
	delete[]	dFp;
	delete[]	rFp;
	delete[]	B;
	delete[]	s;
	delete[]	y;
}



int main2()
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

