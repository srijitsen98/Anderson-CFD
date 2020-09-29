#include<stdio.h>
#include<math.h>
#define dy 1/20.0
#define ny 21
const int n=ny-1;

double x[ny-1];
int tdma(double a, double b, double c, double B[n])
{
    int i,j,k;
    double a1[n],b1[n],c1[n],d1[n],r;
    for(i=1;i<n;i++)
    {
        a1[i]=a;
        b1[i]=b;
        c1[i]=c;
        d1[i]=B[i];
    }
    c1[1]=c1[1]/a1[1];
    d1[1]=d1[1]/a1[1];
    a1[1]=1;

    for(i=2;i<n;i++)
    {
        c1[i]=c1[i]/(a1[i]-c1[i-1]*b1[i]);
        d1[i]=(d1[i]-d1[i-1]*b1[i])/(a1[i]-c1[i-1]*b1[i]);
        a1[i]=1;
    }
    x[n-1]=d1[n-1];
    for(i=n-2;i>=1;i--)
    {
        x[i]=d1[i]-c1[i]*x[i+1];
    }
    return 0;
}
int main()
{
    double a,b,c,A[ny-1][ny-1],B[ny-1],E,dt,Re,u[ny],uan[ny],err;
    int cnt=0;
    int j;
    Re=5000.0;
    E=10.0;
    dt=E*Re*dy*dy;
    a=E+1;
    b=-0.5*E;
    c=-0.5*E;
    //Initial Conditions
    for(j=0;j<ny-1;j++)
    {
        u[j]=0;
    }
    u[ny-1]=1.0;
    //TIme lop
    while(1)
    {
        B[1]=(1-E)*u[1]+0.5*E*(u[0]+u[2])+0.5*E*u[0];
        for(j=2;j<ny-2;j++)
        {
           B[j]=(1-E)*u[j]+0.5*E*(u[j-1]+u[j+1]);
        }
        B[ny-2]=(1-E)*u[ny-2]+0.5*E*(u[ny-1]+u[ny-3])+0.5*E*u[ny-1];
        tdma(a,b,c,B);
        for(j=1;j<ny-1;j++)
        {
            u[j]=x[j];
        }
        err=0;
        for(j=0;j<ny;j++)
        {
            uan[j]=j*dy;
            if(fabs(uan[j]-u[j])/uan[j]>err)
            {
                err=fabs(uan[j]-u[j])/uan[j];
            }
        }
        cnt++;
        printf("\nStep no. %d Time %f Max Error %f",cnt,cnt*dt,err);
        if(err<1e-6)
        {
            break;
        }
        //printf("\nX is \n");
    }
    printf("\nPrint Velocities");
    for(j=0;j<ny;j++)
    {
        printf("\nPosition %f Velocity %f",j*dy,u[j]);
    }
}

