#include<stdio.h>
#include<math.h>
#define ydiv 41
#define pi 3.14159
#define c 0.40825

class var
{
    public:
    double F[50][ydiv][5],G[50][ydiv][5],dFdxi[50][ydiv][5],p[50][ydiv],rho[50][ydiv],u[50][ydiv],v[50][ydiv],M[50][ydiv],T[50][ydiv],xi[50],eta[50][ydiv];
}fg;

double newtonraphson(double f1)
{
    double f,M,err,fdash,C,dy;
    err=1;
    M=1.5;
    f=(1/c)*atan(c*sqrt(pow(M,2.0)-1));
    f=f-atan(sqrt(pow(M,2.0)-1));
    f=f-f1;
    fdash=M/((c*c*(M*M-1)+1)*sqrt((M*M-1)))-1/(M*sqrt(M*M-1));
    while(err>1e-6)
    {
        err=M;
        M=M-f/fdash;
        err-=M;
        err=fabs(err);
        f=(1/c)*atan(c*sqrt(pow(M,2.0)-1));
        f=f-atan(sqrt(pow(M,2.0)-1));
        f=f-f1;
        fdash=M/((c*c*(M*M-1)+1)*sqrt((M*M-1)))-1/(M*sqrt(M*M-1));
    }
    return M;
}

int main()
{
    int i,j,k;
    double xi[50],dxi,h,deta,detadx[ydiv],eta[ydiv],Fpred[ydiv][5],Gpred[ydiv][5];
    double A1,B1,C1,gamma;
    double rhopred[ydiv],dFdxipred[ydiv][5];
    double fi,fact,fcal,psi,ucal,vcal,Mcal,rhocal,Tcal,pcal;
    //double psi,fi,fcal;

    for(i=0;i<50;i++)
    {
        fg.xi[i]=0;
    }

    gamma=1.4;

    i=0;
    fg.xi[0]=12.11;
    h=40+(fg.xi[i]-10)*tan(5.352*pi/180);
    printf("\n%f %f",fg.xi[0],h);
    deta=1.0/(ydiv-1);
    j=1;
    eta[j]=j*deta;
    detadx[j]=(1-eta[j])*tan(5.352*pi/180)/h;
    printf("\n%f %f",eta[j],detadx[j]);
    h=40.2;
    dxi=0.818;
    fg.F[i][j][1]=0.744E3;
    fg.F[i][j+1][1]=0.798E3;
    fg.G[i][j][1]=-0.435E2;
    fg.G[i][j+1][1]=-0.193E2;
    k=1;
    //fg.dFdxi[i][j][k]=detadx[j]*(fg.F[i][j][k]-fg.F[i][j+1][k])/deta;
    //fg.dFdxi[i][j][k]=(1/h)*(fg.G[i][j][k]-fg.G[i][j+1][k])/deta;
    fg.dFdxi[i][j][k]=detadx[j]*(fg.F[i][j][k]-fg.F[i][j+1][k])/deta+(1/h)*(fg.G[i][j][k]-fg.G[i][j+1][k])/deta;
    printf("\n%f",fg.dFdxi[i][j][k]);
    fg.dFdxi[i][j][k]=-28.99;
    Fpred[j][k]=fg.F[i][j][k]+fg.dFdxi[i][j][k]*dxi;
    printf("\n%f",Fpred[j][k]);
    Fpred[j][k]=721.0;
    Fpred[j][2]=585E3;
    Fpred[j][3]=-0.388E5;
    Fpred[j][4]=0.372E9;
    A1=pow(Fpred[j][3],2.0)/(2.0*Fpred[j][1])-Fpred[j][4];
    B1=((gamma)/(gamma-1))*Fpred[j][1]*Fpred[j][2];
    C1=-0.5*((gamma+1)/(gamma-1))*pow(Fpred[j][1],3.0);
    rhopred[j]=(-B1+sqrt(pow(B1,2.0)-4.0*A1*C1))/(2.0*A1);
    printf("\n%f %f %f %f",A1,B1,C1,rhopred[j]);
    rhopred[j]=1.02;
    Gpred[j][1]=rhopred[j]*Fpred[j][3]/Fpred[j][1];
    Gpred[j][2]=Fpred[j][3];
    Gpred[j][3]=rhopred[j]*pow((Fpred[j][3]/Fpred[j][1]),2.0)+Fpred[j][2]-pow(Fpred[j][1],2.0)/rhopred[j];
    Gpred[j][4]=(gamma/(gamma-1))*(Fpred[j][2]-pow(Fpred[j][1],2.0)/rhopred[j])*(Fpred[j][3]/Fpred[j][1])+rhopred[j]*0.5*(Fpred[j][3]/Fpred[j][1])*(pow(Fpred[j][1]/rhopred[j],2.0)+pow(Fpred[j][3]/Fpred[j][1],2.0));
    printf("\n%f",Gpred[j][1]);
    Gpred[j][1]=-0.552E2;
    Gpred[j-1][1]=-0.658E2;
    Fpred[j-1][1]=703;
    dFdxipred[j][k]=detadx[j]*(Fpred[j-1][k]-Fpred[j][k])/deta+(1/h)*(Gpred[j-1][k]-Gpred[j][k])/deta;
    printf("\n%f",dFdxipred[j][k]);
    dFdxipred[j][k]=-0.122E2;
    fg.dFdxi[i][j][k]=(fg.dFdxi[i][j][k]+dFdxipred[j][k])*0.5;
    printf("\n%f",fg.dFdxi[i][j][k]);
    fg.dFdxi[i][j][k]=-20.5;
    fg.F[i+1][j][k]=fg.F[i][j][k]+fg.dFdxi[i][j][k]*dxi;
    printf("\n%f",fg.F[i+1][j][k]);
    fg.F[i+1][j][k]=0.728E3;//after visc;
    fg.F[i+1][j][2]=0.590E6;
    fg.F[i+1][j][3]=-0.36E5;
    fg.F[i+1][j][4]=0.375E9;
    i=i+1;
    A1=pow(fg.F[i][j][3],2.0)/(2.0*fg.F[i][j][1])-fg.F[i][j][4];
    B1=((gamma)/(gamma-1))*fg.F[i][j][1]*fg.F[i][j][2];
    C1=-0.5*((gamma+1)/(gamma-1))*pow(fg.F[i][j][1],3.0);
    fg.rho[i][j]=(-B1+sqrt(pow(B1,2.0)-4.0*A1*C1))/(2.0*A1);
    printf("\n%f",fg.rho[i][j]);
    fg.rho[i][j]=1.04;
    fg.u[i][j]=fg.F[i][j][1]/fg.rho[i][j];
    fg.v[i][j]=fg.F[i][j][3]/fg.F[i][j][1];
    printf("\n%f",fg.u[i][j]);
    printf("\n%f",fg.v[i][j]);
    fg.u[i][j]=701;
    fg.v[i][j]=-49.4;
    fg.p[i][j]=fg.F[i][j][2]-fg.F[i][j][1]*fg.u[i][j];
    fg.T[i][j]=fg.p[i][j]/(287.0*fg.rho[i][j]);
    printf("\n%f",fg.p[i][j]);
    printf("\n%f",fg.T[i][j]);
    i=i-1;
    j=0;
    fg.dFdxi[i][j][k]=-26.1;
    eta[j]=j*deta;
    detadx[j]=(1-eta[j])*tan(5.352*pi/180)/h;
    dFdxipred[0][k]=detadx[1]*(Fpred[0][k]-Fpred[1][k])/deta+(1/h)*(Gpred[0][k]-Gpred[1][k])/deta;
    printf("\n%f",dFdxipred[0][k]);
    dFdxipred[j][k]=-0.1218E2;
    fg.dFdxi[i][j][k]=0.5*(fg.dFdxi[i][j][k]+dFdxipred[j][k]);
    printf("\n%f",fg.dFdxi[i][j][k]);
    //fg.dFdxi[i][j][k]=-20.5;
    fg.F[i][j][k]=696;
    fg.F[i+1][j][k]=fg.F[i][j][k]+fg.dFdxi[i][j][k]*dxi;
    printf("\n%f",fg.F[i+1][j][k]);
    //*/

    i=i+1;
    printf("\n%d",i);
    fg.xi[1]=12.928;
    vcal=-74.6;
    ucal=707.0;

    if(fg.xi[i]<10.0)
    {
        fi=atan(vcal/ucal);
    }

    else
    {
        //psi=atan2(fabs(vcal),ucal);
        //psi=atan(74.6/ucal);
        psi=atan2(fabs(vcal),ucal);
        printf("\n%f %f",vcal,ucal);
        //psi=0;
        fi=5.352*pi/180.0-psi;
    }



    //psi=atan2(fabs(vcal),ucal);
    //fi=5.352*pi/180.0-psi;
    //printf("\n%f %f",psi,fi);
    Mcal=2.22;//fg.M[i][j];
    pcal=0.705E5;//fg.p[i][j];
    Tcal=255;//fg.T[i][j];
    rhocal=0.963;
    fcal=sqrt((gamma-1)*(pow(2.22,2.0)-1)/(gamma+1));
    fcal=atan(fcal);
    fcal=sqrt((gamma+1)/(gamma-1))*fcal;
    fcal=fcal-atan(sqrt(pow(2.22,2.0)-1.0));
    //fcal=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)*(pow(2.22,2.0)-1)/(gamma+1)))-atan(sqrt(pow(2.22,2.0)-1.0));
    //fcal=fcal+fi;
    printf("\n%f %f",fcal,fcal*180/pi);
    fact=fcal+fi;
    printf("\n%f %f",fact,fact*180/pi);

    fg.M[i][j]=newtonraphson(fact);
    printf("\n%f",fg.M[i][j]);
    /*
    fg.p[i][j]=pcal*pow((1+(gamma-1)*0.5*Mcal*Mcal)/(1+(gamma-1)*0.5*fg.M[i][j]*fg.M[i][j]),gamma/(gamma-1));
    fg.T[i][j]=Tcal*(1+(gamma-1)*0.5*Mcal*Mcal)/(1+(gamma-1)*0.5*fg.M[i][j]*fg.M[i][j]);
    fg.rho[i][j]=fg.p[i][j]/(287.0*fg.T[i][j]);
    */
}
