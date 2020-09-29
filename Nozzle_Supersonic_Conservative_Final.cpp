#include<stdio.h>
#include<math.h>
//#include<conio.h>
#define nx 31
double V[nx],T[nx],rho[nx],A[nx],M[nx],p[nx];
double Van[nx],Tan[nx],rhoan[nx],Man[nx],pan[nx];
double F1[nx],F2[nx],F3[nx];
double U1[nx],U2[nx],U3[nx];
double U1o[nx],U2o[nx],U3o[nx];
double dU1[nx],dU2[nx],dU3[nx];
double dU1b[nx],dU2b[nx],dU3b[nx];
double U1b[nx],U2b[nx],U3b[nx];
double F1b[nx],F2b[nx],F3b[nx];
double Vo[nx],To[nx],rhoo[nx];
double Vean[nx],Tean[nx],rhoean[nx],pean[nx],Mean[nx];
int s[nx];
double Vb[nx],Tb[nx],rhob[nx];
double drho[nx],dV[nx],dT[nx];
double drhob[nx],dVb[nx],dTb[nx];
double dx,dt,As;
const double C = 0.5;
const double gamma=1.4;
int i;
int displayuf()
{
	for(i=0;i<nx;i++)
	{
		printf("\nU1 %f U2 %f U3 %f F1 %f F2 %f F3 %f",U1[i],U2[i],U3[i],F1[i],F2[i],F3[i]);
	}
	return 0;

}

int areainit()
{
    As=0;
	printf("Area Initialization");
    for(i=0;i<nx;i++)
    {
        A[i]=1+2.2*(i*dx-1.5)*(i*dx-1.5);
		s[i]=0;
    }
	for(i=0;i<nx;i++)
	{
		if((A[i]<A[i-1])&&(A[i+1]>A[i]))
        {
            As=A[i];
        }
		if(As!=0)
		{
			s[i]=1;
		}
		printf("\n%f %f %d",A[i],As,s[i]);
	}
	return 0;
}
double bisection(double ar,int s1)
{
	//s1=0 subsonic
	//s1=1 supersonic
	double c1,c2;
	c1=gamma+1;
	c2=gamma-1;
	double Mh, Ml, Mm, fh, fl, fm;
    if(s1==1)
	{
		Mh=5.0;
		Ml=1.00;
	}
	if(s1==0)
	{
		Mh=0.999;
		Ml=0.001;
	}
	fh=(1/(Mh*Mh))*pow((2/c1)*(1+0.5*c2*Mh*Mh),c1/c2)-ar*ar;
    fl=(1/(Ml*Ml))*pow((2/c1)*(1+0.5*c2*Ml*Ml),c1/c2)-ar*ar;
	//printf("\n%f",fh*fl);
	Mm=8.0;
	fm=8.0;
	//printf("\n%f",fm);
	while(fabs(fm)>1e-12)
	{
		//printf("\nHello");
		Mm=(Mh+Ml)/2.0;
		fm=(1/(Mm*Mm))*pow((2/c1)*(1+0.5*c2*Mm*Mm),c1/c2)-ar*ar;
		if(fm*fh<0)
		{
			fl=fm;
			Ml=Mm;
		}
		else
		{
			fh=fm;
			Mh=Mm;
		}
	}
	//printf("\n%f",fh*fl);
	return Mm;

}

int analytical()
{
	printf("\nAnalytical Solution");
    for(i=0;i<nx;i++)
    {
		Man[i]=bisection(A[i]/As,s[i]);
        //M[i]=newtonraphson(A[i]/As);
        pan[i]=pow((1+(gamma-1)*0.5*Man[i]*Man[i]),-gamma/(gamma-1));
        rhoan[i]=pow((1+(gamma-1)*0.5*Man[i]*Man[i]),-1.0/(gamma-1));
        Tan[i]=pow((1+(gamma-1)*0.5*Man[i]*Man[i]),-1.0);
        Van[i]=Man[i]*sqrt(Tan[i]);
		printf("\nM %f p %f rho %f T %f V %f A %f",Man[i],pan[i],rhoan[i],Tan[i],Van[i],A[i]);
    }
	return 0;
}
int varinit()
{
	printf("\nInitial Conditions");
    for(i=0;i<nx;i++)
    {

        if(i*dx<=0.51)
		{
			rho[i]=1.0;
			T[i]=1.0;
		}
		else if(i*dx<=1.51)
		{
			rho[i]=1.0-0.366*(i*dx-0.5);
			T[i]=1.0-0.167*(i*dx-0.5);
		}
		else if(i*dx<=3.01)
		{
			rho[i]=0.634-0.3879*(i*dx-1.5);
			T[i]=0.833-0.3507*(i*dx-1.5);
		}
		V[i]=0.59/(rho[i]*A[i]);

		//V[i]=Van[i];
		//rho[i]=rhoan[i];
		//T[i]=Tan[i];
		M[i]=V[i]/sqrt(T[i]);
        p[i]=rho[i]*T[i];
		U1[i]=rho[i]*A[i];
		U2[i]=rho[i]*A[i]*V[i];
		U3[i]=rho[i]*(T[i]/(gamma-1)+(gamma*0.5)*V[i]*V[i])*A[i];
		F1[i]=U2[i];
		F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
		F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);
		printf("\nM %f p %f rho %f T %f V %f A %f U1 %f U2 %f U3 %f",M[i],p[i],rho[i],T[i],V[i],A[i],U1[i],U2[i],U3[i]);

	}
	//displayuf();
	//getch();
	return 0;
}
double courant()
{
    double dtt;
    dtt=1000;
    for(i=0;i<nx;i++)
    {
        if(C*dx/(sqrt(T[i])+V[i])<dtt)
        {
           dtt=C*dx/(sqrt(T[i])+V[i]);
        }
		//printf("\ndtt is %f",dtt);
    }
    return dtt;
}
int mac()
{
    //Predictor Step
    for(i=1;i<nx-1;i++)
    {
        dU1[i]=-(F1[i+1]-F1[i])/dx;
		dU2[i]=-(F2[i+1]-F2[i])/dx+(1/gamma)*rho[i]*T[i]*((A[i+1])-(A[i]))/dx;
		dU3[i]=-(F3[i+1]-F3[i])/dx;
	}
    for(i=1;i<nx-1;i++)
    {
        U1b[i]=U1[i]+dU1[i]*dt;
		U2b[i]=U2[i]+dU2[i]*dt;
		U3b[i]=U3[i]+dU3[i]*dt;
		F1b[i]=U2b[i];
		F2b[i]=(U2b[i]*U2b[i])/U1b[i]+((gamma-1)/gamma)*(U3b[i]-(gamma*0.5)*U2b[i]*U2b[i]/U1b[i]);
		F3b[i]=gamma*U3b[i]*U2b[i]/U1b[i]-gamma*(gamma-1)*0.5*pow(U2b[i],3.0)/pow(U1b[i],2.0);
		rhob[i]=U1b[i]/A[i];
		Tb[i]=(gamma-1)*((U3b[i]/U1b[i])-(gamma/2.0)*(U2b[i]/U1b[i])*(U2b[i]/U1b[i]));
	}
	F1b[0]=F1[0];
	F2b[0]=F2[0];
	F3b[0]=F3[0];
	//Corrector Step
    for(i=1;i<nx-1;i++)
    {
        dU1b[i]=-(F1b[i]-F1b[i-1])/dx;
		dU2b[i]=-(F2b[i]-F2b[i-1])/dx+(1/gamma)*rhob[i]*Tb[i]*((A[i])-(A[i-1]))/dx;
		dU3b[i]=-(F3b[i]-F3b[i-1])/dx;
    }
    for(i=1;i<nx-1;i++)
    {
        U1[i]=U1[i]+(dU1b[i]+dU1[i])*0.5*dt;
		U2[i]=U2[i]+(dU2b[i]+dU2[i])*0.5*dt;
		U3[i]=U3[i]+(dU3b[i]+dU3[i])*0.5*dt;
		//printf("\nU1 %f U2 %f U3 %f",U1[i],U2[i],U3[i]);
		F1[i]=U2[i];
		F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
		F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);
		rho[i]=U1[i]/A[i];
		V[i]=U2[i]/U1[i];
		T[i]=(gamma-1)*((U3[i]/U1[i])-(gamma/2.0)*(U2[i]/U1[i])*(U2[i]/U1[i]));
    }
	//displayuf();
	return 0;

}
int boundary()
{
    //Left Boundary
    U1[0]=A[0];
    U2[0]=2*U2[1]-U2[2];
    V[0]=U2[0]/U1[0];
	U3[0]=U1[0]*(T[0]/(gamma-1)+(gamma/2.0)*V[0]*V[0]);
    i=0;
	F1[i]=U2[i];
	F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
	F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);

	//Right Boundary
	U1[nx-1]=2.0*U1[nx-2]-U1[nx-3];
    U2[nx-1]=2.0*U2[nx-2]-U2[nx-3];
    U3[nx-1]=2.0*U3[nx-2]-U3[nx-3];
	i=nx-1;
	F1[i]=U2[i];
	F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
	F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);
	rho[i]=U1[i]/A[i];
	V[i]=U2[i]/U1[i];
	T[i]=(gamma-1)*((U3[i]/U1[i])-(gamma/2.0)*(U2[i]/U1[i])*(U2[i]/U1[i]));
	return 0;
}
double para()
{
    for(i=0;i<nx;i++)
    {
        p[i]=rho[i]*T[i];
        M[i]=V[i]/sqrt(T[i]);
    }
	return 0;
}

double errorrho()
{
	double errrho;
	errrho=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(rho[i]-rhoo[i])>errrho)
		{
			errrho=fabs(rho[i]-rhoo[i]);
		}
	}
	return errrho;
}
double errorv()
{
	double errv;
	errv=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(V[i]-Vo[i])>errv)
		{
			errv=fabs(V[i]-Vo[i]);
		}
	}
	return errv;
}
double errort()
{
	double errt;
	errt=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(T[i]-To[i])>errt)
		{
			errt=fabs(T[i]-To[i]);
		}
	}
	return errt;
}

double errorU1()
{
	double errU1;
	errU1=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(U1[i]-U1o[i])>errU1)
		{
			errU1=fabs(U1[i]-U1o[i]);
		}
	}
	return errU1;
}
double errorU2()
{
	double errU2;
	errU2=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(U2[i]-U2o[i])>errU2)
		{
			errU2=fabs(U2[i]-U2o[i]);
		}
	}
	return errU2;
}
double errorU3()
{
	double errU3;
	errU3=0;
	for(i=1;i<nx-1;i++)
	{
		if(fabs(U3[i]-U3o[i])>errU3)
		{
			errU3=fabs(U3[i]-U3o[i]);
		}
	}
	return errU3;
}

int erroranalytical()
{
	for(i=0;i<nx;i++)
	{
		Tean[i]=100.0*fabs(T[i]-Tan[i])/T[i];
		Mean[i]=100.0*fabs(M[i]-Man[i])/M[i];
		Vean[i]=100.0*fabs(V[i]-Van[i])/V[i];
		pean[i]=100.0*fabs(p[i]-pan[i])/p[i];
		rhoean[i]=100.0*fabs(rho[i]-rhoan[i])/rho[i];
	}
	return 0;
}
int test()
{
	double dU1[nx],dU2[nx],dU3[nx];
	double dU1b[nx],dU2b[nx],dU3b[nx];
	double U1b[nx],U2b[nx],U3b[nx];
	double F1b[nx],F2b[nx],F3b[nx];
	i=5;
	/*
	rho[5]=5.0;
	A[5]=6.0;
	V[5]=7.0;
	T[5]=8.0;
	U1[i]=rho[i]*A[i];
	U2[i]=rho[i]*A[i]*V[i];
	U3[i]=rho[i]*(T[i]/(gamma-1)+(gamma*0.5)*V[i]*V[i])*A[i];
	//U1[i]=1.5;
	//U2[i]=2.5;
	//U3[i]=3.5
	F1[i]=U2[i];
	F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
	F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);
	printf("\nU1 %f U2 %f U3 %f F1 %f F2 %f F3 %f",U1[i],U2[i],U3[i],F1[i],F2[i],F3[i]);
	*/
	dt=0.5*dx;
	printf("\ndt %f dx %f",dt,dx);
	U1[i]=5.0;
	U2[i]=6.0;
	U3[i]=7.0;
	F1[i]=9.0;
	F2[i]=10.0;
	F3[i]=11.0;
	F1[i+1]=7.0;
	F2[i+1]=14.0;
	F3[i+1]=18.0;

	rho[i]=4.0;
	T[i]=2.0;
	A[i]=3.0;
	A[i+1]=5.0;
	dU1[i]=-(F1[i+1]-F1[i])/dx;
	dU2[i]=-(F2[i+1]-F2[i])/dx+(1/gamma)*rho[i]*T[i]*(log(A[i+1])-log(A[i]))/dx;
	dU3[i]=-(F3[i+1]-F3[i])/dx;
	printf("\ndU1 %f dU2 %f dU3 %f",dU1[i],dU2[i],dU3[i]);
	U1b[i]=U1[i]+dU1[i]*dt;
	U2b[i]=U2[i]+dU2[i]*dt;
	U3b[i]=U3[i]+dU3[i]*dt;
	printf("\nU1b %f U2b %f U3b %f",U1b[i],U2b[i],U3b[i]);
	F1b[i]=U2b[i];
	F2b[i]=(U2b[i]*U2b[i])/U1b[i]+((gamma-1)/gamma)*(U3b[i]-(gamma*0.5)*U2b[i]*U2b[i]/U1b[i]);
	F3b[i]=gamma*U3b[i]*U2b[i]/U1b[i]-gamma*(gamma-1)*0.5*pow(U2b[i],3.0)/pow(U1b[i],2.0);
	rhob[i]=U1b[i]/A[i];
	Tb[i]=(gamma-1)*((U3b[i]/U1b[i])-(gamma/2.0)*(U2b[i]/U1b[i])*(U2b[i]/U1b[i]));
	printf("\nF1b %f F2b %f F3b %f rhob %f Tb %.9f",F1b[i],F2b[i],F3b[i],rhob[i],Tb[i]);
	F1b[i]=5;
	F2b[i]=3;
	F3b[i]=6;
	F1b[i-1]=4;
	F2b[i-1]=5;
	F3b[i-1]=1;
	rhob[i]=1.5;
	Tb[i]=2;
	A[i]=3;
	A[i-1]=2;
	dU1b[i]=-(F1b[i]-F1b[i-1])/dx;
	dU2b[i]=-(F2b[i]-F2b[i-1])/dx+(1/gamma)*rhob[i]*Tb[i]*(log(A[i])-log(A[i-1]))/dx;
	dU3b[i]=-(F3b[i]-F3b[i-1])/dx;
	printf("\ndU1b %f dU2b %f dU3b %f",dU1b[i],dU2b[i],dU3b[i]);
	U1[i]=U1[i]+(dU1b[i]+dU1[i])*0.5*dt;
	U2[i]=U2[i]+(dU2b[i]+dU2[i])*0.5*dt;
	U3[i]=U3[i]+(dU3b[i]+dU3[i])*0.5*dt;
	printf("\nU1 %f U2 %f U3 %f",U1[i],U2[i],U3[i]);
	F1[i]=U2[i];
	F2[i]=(U2[i]*U2[i])/U1[i]+((gamma-1)/gamma)*(U3[i]-(gamma*0.5)*U2[i]*U2[i]/U1[i]);
	F3[i]=gamma*U3[i]*U2[i]/U1[i]-gamma*(gamma-1)*0.5*pow(U2[i],3.0)/pow(U1[i],2.0);
	rho[i]=U1[i]/A[i];
	V[i]=U2[i]/U1[i];
	T[i]=(gamma-1)*((U3[i]/U1[i])-(gamma/2.0)*(U2[i]/U1[i])*(U2[i]/U1[i]));
	//getch();
	return 0;
}
int main()
{
    int it=0;
    double errrho,errt,errv;
    double errU1,errU2,errU3;
	dx=3.0/(nx-1);
    //test
	//test();
	//Area initialization
    areainit(); //checked
    //Analytical Solution
    analytical(); //checked
    //Variable initialization
    varinit(); //checked
	errU1=1;
	errU2=2;
	errU3=3;
    errrho=1;
	errv=1;
	errt=1;
	while(errrho>1e-6||errv>1e-6||errt>1e-6)
	//while(errU1>1e-6||errU2>1e-6||errU3>1e-6)
    {
		//old = new
		for(i=0;i<nx;i++)
		{
			rhoo[i]=rho[i];
			Vo[i]=V[i];
			To[i]=T[i];
			//U1o[i]=U1[i];
			//U2o[i]=U2[i];
			//U3o[i]=U3[i];
		}
        //Selection of dt
        dt=courant(); //checked
        //MacCormack
        mac(); //checked
        //Boundary
        boundary(); //looks correct
        //Residual Parameters
        para(); //looks correct
        //Iteration Upgrade
        errrho=errorrho();
		errt=errort();
		errv=errorv();
		//errU1=errorU1();
		//errU2=errorU2();
		//errU3=errorU3();
		it++;
		//printf("\nIteration number %d Error U1 %f Error U2 %f Error U3 %f",it,errU1,errU2,errU3);
		printf("\nIteration number %d Error rho %f Error V %f Error T %f",it,errrho,errv,errt);
		/*
		for(i=0;i<nx;i++)
			printf("\nM %f p %f rho %f T %f V %f A %f",M[i],p[i],rho[i],T[i],V[i],A[i]);
		*/
    }
	erroranalytical();
	printf("\nFinal Solution (E denotes percentage error in variable when compared with analytical solution)");
	for(i=0;i<nx;i++)
		printf("\nM %f E %f p %f E %f rho %f E %f T %f E %f V %f E %f A %f",M[i],Mean[i],p[i],pean[i],rho[i],rhoean[i],T[i],Tean[i],V[i],Vean[i],A[i]);

	//getch();
	return 0;
}
