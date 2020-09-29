#include<stdio.h>
#include<math.h>
#define nx 31
double V[nx],T[nx],rho[nx],A[nx],M[nx],p[nx];
double Van[nx],Tan[nx],rhoan[nx],Man[nx],pan[nx];
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
        rho[i]=1-0.3146*i*dx;
        T[i]=1-0.2314*i*dx;
        V[i]=(0.1+1.09*i*dx)*pow(T[i],0.5);
        M[i]=V[i]/sqrt(T[i]);
        p[i]=rho[i]*T[i];
		printf("\nM %f p %f rho %f T %f V %f A %f",M[i],p[i],rho[i],T[i],V[i],A[i]);
    }
	return 0;
}
double courant()
{
    double dtt;
    dtt=1000;
    for(i=1;i<nx-1;i++)
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
        //Rho
        drho[i]=-rho[i]*(V[i+1]-V[i])/dx-rho[i]*V[i]*(log(A[i+1])-log(A[i]))/dx-V[i]*(rho[i+1]-rho[i])/dx;
        //Velocity
        dV[i]=-V[i]*(V[i+1]-V[i])/dx-(1/gamma)*((T[i+1]-T[i])/dx+(T[i]/rho[i])*(rho[i+1]-rho[i])/dx);
        //Temperature
        dT[i]=-V[i]*(T[i+1]-T[i])/dx-(gamma-1)*T[i]*((V[i+1]-V[i])/dx+V[i]*(log(A[i+1])-log(A[i]))/dx);
    }
    for(i=1;i<nx-1;i++)
    {
        rhob[i]=rho[i]+drho[i]*dt;
        Vb[i]=V[i]+dV[i]*dt;
        Tb[i]=T[i]+dT[i]*dt;
    }
	rhob[0]=rho[0];
	Vb[0]=V[0];
	Tb[0]=T[0];
    //Corrector Step
    for(i=1;i<nx-1;i++)
    {
        //Rho
        drhob[i]=-rhob[i]*(Vb[i]-Vb[i-1])/dx-rhob[i]*Vb[i]*(log(A[i])-log(A[i-1]))/dx-Vb[i]*(rhob[i]-rhob[i-1])/dx;
        //Velocity
        dVb[i]=-Vb[i]*(Vb[i]-Vb[i-1])/dx-(1/gamma)*((Tb[i]-Tb[i-1])/dx+(Tb[i]/rhob[i])*(rhob[i]-rhob[i-1])/dx);
        //Temperature
        dTb[i]=-Vb[i]*(Tb[i]-Tb[i-1])/dx-(gamma-1)*Tb[i]*((Vb[i]-Vb[i-1])/dx+Vb[i]*(log(A[i])-log(A[i-1]))/dx);
    }
    for(i=1;i<nx-1;i++)
    {
        drho[i]=0.5*(drho[i]+drhob[i]);
        dV[i]=0.5*(dV[i]+dVb[i]);
        dT[i]=0.5*(dT[i]+dTb[i]);
        rho[i]=rho[i]+drho[i]*dt;
        V[i]=V[i]+dV[i]*dt;
        T[i]=T[i]+dT[i]*dt;
    }
	return 0;

}
int boundary()
{
    //Left Boundary
    rho[0]=1.0;
    T[0]=1.0;
    V[0]=2*V[1]-V[2];
    //Right Boundary
    rho[nx-1]=2.0*rho[nx-2]-rho[nx-3];
    T[nx-1]=2.0*T[nx-2]-T[nx-3];
    V[nx-1]=2.0*V[nx-2]-V[nx-3];
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
int main()
{
    int it=0;
    double errrho,errt,errv;
    dx=3.0/(nx-1);
    //Area initialization
    areainit(); //checked
    //Analytical Solution
    analytical(); //checked
    //Variable initialization
    varinit(); //checked
	errrho=1;
	errv=2;
	errt=3;
    while(errrho>1e-6||errv>1e-6||errt>1e-6)
    {
		//old = new
		for(i=0;i<nx;i++)
		{
			rhoo[i]=rho[i];
			Vo[i]=V[i];
			To[i]=T[i];
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
		it++;
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


	return 0;
}
