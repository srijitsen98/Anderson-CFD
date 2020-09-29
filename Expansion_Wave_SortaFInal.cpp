#include<stdio.h>
#include<math.h>
#define ydiv 41
#define pi 3.14159
#define theta_1 5.352*pi/180.0
#define c 0.40825

class var
{
    public:
    double F[100][ydiv][5],G[100][ydiv][5],dFdxi[100][ydiv][5],p[100][ydiv],rho[100][ydiv],u[100][ydiv],v[100][ydiv],M[100][ydiv],T[100][ydiv],xi[100],eta[100][ydiv];
}fg;

double newtonraphson(double f1) //updated and working
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
double selectdx(double dy, int i) //working when last checked on this code
{
    double val,dx,theta,mu;
    int j;
    val=0;
    for(j=1;j<ydiv-1;j++)
    {
        theta=atan(fg.v[i][j]/fg.u[i][j]);
        mu=asin(1/fg.M[i][j]);
        if(fabs(tan(theta+mu))>val)
        {
            val=fabs(tan(theta+mu));
        }
        if(tan(theta-mu)>val)
        {
            val=fabs(tan(theta-mu));
        }
    }
    dx=0.5*dy/val;
    return dx;
}
int main()
{
    double gamma,a,dxi,dy,deta,h,detadx[ydiv],eta[ydiv];
    double Fpred[ydiv][5],upred[ydiv],ppred[ydiv],rhopred[ydiv],Gpred[ydiv][5],dFdxipred[ydiv][5];
    double A1, B1, C1;
    double psi,Mcal,ucal,vcal,fi,pcal,Tcal,fcal,fact;
    int i,j,k;
    double C, SF, Cy; //Courant Number
    Cy=0.6;
    C=0.5;
    gamma=1.4;
    dy=40.0/(ydiv-1);
    //Initial Conditions
    //M_in = 2
    //p_in = 1.01*10^5
    //T_in = 286.1
    //rho_in = 1.23;
    i=0;
    printf("\nInitial Conditions");
    fg.xi[i]=0;
    for(j=0;j<ydiv;j++)
    {
        fg.eta[i][j]=j*dy;
        fg.p[i][j]=1.01E5;
        fg.T[i][j]=286.1;
        fg.M[i][j]=2;
        fg.rho[i][j]=1.23;
        a=sqrt(gamma*fg.p[i][j]/fg.rho[i][j]);
        fg.u[i][j]=fg.M[i][j]*a;
        fg.v[i][j]=0;
        fg.F[i][j][1]=fg.rho[i][j]*fg.u[i][j];
        fg.F[i][j][2]=fg.rho[i][j]*pow(fg.u[i][j],2.0)+fg.p[i][j];
        fg.F[i][j][3]=fg.rho[i][j]*fg.u[i][j]*fg.v[i][j];
        fg.F[i][j][4]=(gamma/(gamma-1))*fg.p[i][j]*fg.u[i][j]+fg.rho[i][j]*fg.u[i][j]*(pow(fg.u[i][j],2.0)+pow(fg.v[i][j],2.0))*0.5;
        fg.G[i][j][1]=fg.rho[i][j]*fg.v[i][j];
        fg.G[i][j][2]=fg.rho[i][j]*fg.u[i][j]*fg.v[i][j];
        fg.G[i][j][3]=fg.rho[i][j]*pow(fg.v[i][j],2.0)+fg.p[i][j];
        fg.G[i][j][4]=(gamma/(gamma-1))*fg.p[i][j]*fg.v[i][j]+fg.rho[i][j]*fg.v[i][j]*(pow(fg.u[i][j],2.0)+pow(fg.v[i][j],2.0))*0.5;
        printf("\n%f %d %d %f %f %f %f %f %f %f %f %f %f",fg.xi[i],i,j,fg.u[i][0],fg.v[i][0],fg.rho[i][0],fg.p[i][0],fg.T[i][0],fg.M[i][0],fg.F[i][0][1],fg.F[i][0][2],fg.F[i][0][3],fg.F[i][0][4]);
//printf("\n%d %d %f %f %f %f %f %f %f %f",i,j,fg.F[i][j][1],fg.F[i][j][2],fg.F[i][j][3],fg.F[i][j][4],fg.G[i][j][1],fg.G[i][j][2],fg.G[i][j][3],fg.G[i][j][4]);
    }

    //function to select dx in terms of dy
    dxi=selectdx(1.0,0);
    deta=1.0/(ydiv-1);
    printf("\n%f %f",dxi,dy);
    //xi is modified x
    //eta is modified y
    while(1)
    {
        if(fg.xi[i]<10)
        {
            h=40.0;
            for(j=0;j<ydiv;j++)
            {
                detadx[j]=0;
            }
        }
        if(fg.xi[i]>10)
        {
            h=40+(fg.xi[i]-10.0)*tan(5.352*pi/180);
            for(j=0;j<ydiv;j++)
            {
                eta[j]=j*deta;
                detadx[j]=(1-eta[j])*tan(theta_1)/h;
                //printf("\ndetadx[%d] %f",j,detadx[j]);
            }
        }


        //dy=h/(ydiv-1);
        //dxi=selectdx(dy,0);
        //predictor step
        //F differentials
        for(j=0;j<ydiv-1;j++)
        {
            for(k=1;k<=4;k++)
            {
                //double temp=predictorf(i,j,k,h,detadx,deta);
                fg.dFdxi[i][j][k]=detadx[j]*(fg.F[i][j][k]-fg.F[i][j+1][k])/deta+(1/h)*(fg.G[i][j][k]-fg.G[i][j+1][k])/deta;
            }
            //printf("\nFdifferential1 %d %d %f %f %f %f",i,j,fg.dFdxi[i][j][1],fg.dFdxi[i][j][2],fg.dFdxi[i][j][3],fg.dFdxi[i][j][4]);
        }

        //F predicted values
        for(j=0;j<ydiv-1;j++)
        {
            for(k=1;k<=4;k++)
            {
                Fpred[j][k]=fg.F[i][j][k]+fg.dFdxi[i][j][k]*dxi;
                if(j>0)
                {
                    SF=Cy*(fg.F[i][j+1][k]-2*fg.F[i][j][k]+fg.F[i][j-1][k])*fabs(fg.p[i][j+1]-2*fg.p[i][j]+fg.p[i][j-1])/(fg.p[i][j+1]+2*fg.p[i][j]+fg.p[i][j-1]);
                    Fpred[j][k]=Fpred[j][k]+SF;
                }
            }
            //printf("\nFpred %d %d %f %f %f %f %f",i,j,Fpred[j][1],Fpred[j][2],Fpred[j][3],Fpred[j][4],SF);
            A1=pow(Fpred[j][3],2.0)/(2.0*Fpred[j][1])-Fpred[j][4];
            B1=((gamma)/(gamma-1))*Fpred[j][1]*Fpred[j][2];
            C1=-0.5*((gamma+1)/(gamma-1))*pow(Fpred[j][1],3.0);
            rhopred[j]=(-B1+sqrt(pow(B1,2.0)-4.0*A1*C1))/(2.0*A1);
            Gpred[j][1]=rhopred[j]*Fpred[j][3]/Fpred[j][1];
            Gpred[j][2]=Fpred[j][3];
            Gpred[j][3]=rhopred[j]*pow((Fpred[j][3]/Fpred[j][1]),2.0)+Fpred[j][2]-pow(Fpred[j][1],2.0)/rhopred[j];
            Gpred[j][4]=(gamma/(gamma-1))*(Fpred[j][2]-pow(Fpred[j][1],2.0)/rhopred[j])*(Fpred[j][3]/Fpred[j][1])+rhopred[j]*0.5*(Fpred[j][3]/Fpred[j][1])*(pow(Fpred[j][1]/rhopred[j],2.0)+pow(Fpred[j][3]/Fpred[j][1],2.0));
            upred[j]=Fpred[j][1]/rhopred[j];
            ppred[j]=Fpred[j][2]-Fpred[j][1]*upred[j];
            //printf("\nFGpred %d %d %f %f %f %f %f %f %f %f %f",i,j,Fpred[j][1],Fpred[j][2],Fpred[j][3],Fpred[j][4],Gpred[j][1],Gpred[j][2],Gpred[j][3],Gpred[j][4]);
        }
        for(k=1;k<=4;k++)
        {
            Fpred[ydiv-1][k]=fg.F[i][ydiv-1][k];
            Gpred[ydiv-1][k]=fg.G[i][ydiv-1][k];
        }
        //printf("\nFpred %d 0 %f %f %f %f %f",i,Fpred[0][1],Fpred[0][2],Fpred[0][3],Fpred[0][4]);
        //corrector step
        //F differentials
        for(k=1;k<=4;k++)
        {
            dFdxipred[0][k]=detadx[0]*(Fpred[0][k]-Fpred[1][k])/deta+(1/h)*(Gpred[0][k]-Gpred[1][k])/deta;
        }
        for(j=1;j<ydiv-1;j++)
        {
            for(k=1;k<=4;k++)
            {
                //double temp=predictorf(i,j,k,h,detadx,deta);
                dFdxipred[j][k]=detadx[j]*(Fpred[j-1][k]-Fpred[j][k])/deta+(1/h)*(Gpred[j-1][k]-Gpred[j][k])/deta;
            }
            //printf("\nFdiffcorr%d %f %f %f %f",j,dFdxipred[j][1],dFdxipred[j][2],dFdxipred[j][3],dFdxipred[j][4]);
        }
        for(j=0;j<ydiv-1;j++)
        {
            for(k=1;k<=4;k++)
            {
                fg.dFdxi[i][j][k]=(fg.dFdxi[i][j][k]+dFdxipred[j][k])*0.5;
                //printf("\nFdiffcorr%d %f",j,fg.dFdxi[i][j][k]);
                fg.F[i+1][j][k]=fg.F[i][j][k]+fg.dFdxi[i][j][k]*dxi;
                //printf("\nFdiffcorr%d %f",j,fg.F[i+1][j][k]);

                if(j>0)
                {
                    SF=Cy*(Fpred[j+1][k]-2*Fpred[j][k]+Fpred[j-1][k])*fabs(ppred[j+1]-2*ppred[j]+ppred[j-1])/(ppred[j+1]+2*ppred[j]+ppred[j-1]);
                    fg.F[i+1][j][k]=fg.F[i+1][j][k]+SF;
                }
            }
        }
        i=i+1;

        //Other variables
        for(j=0;j<ydiv-1;j++)
        {
            A1=pow(fg.F[i][j][3],2.0)/(2.0*fg.F[i][j][1])-fg.F[i][j][4];
            B1=((gamma)/(gamma-1))*fg.F[i][j][1]*fg.F[i][j][2];
            C1=-0.5*((gamma+1)/(gamma-1))*pow(fg.F[i][j][1],3.0);
            fg.rho[i][j]=(-B1+sqrt(pow(B1,2.0)-4.0*A1*C1))/(2.0*A1);
            fg.u[i][j]=fg.F[i][j][1]/fg.rho[i][j];
            fg.v[i][j]=fg.F[i][j][3]/fg.F[i][j][1];
            fg.p[i][j]=fg.F[i][j][2]-fg.F[i][j][1]*fg.u[i][j];
            fg.T[i][j]=fg.p[i][j]/(287.0*fg.rho[i][j]);
            fg.M[i][j]=sqrt(pow(fg.u[i][j],2.0)+pow(fg.v[i][j],2.0));
            a=sqrt(gamma*fg.p[i][j]/fg.rho[i][j]);
            fg.M[i][j]=fg.M[i][j]/a;

            fg.G[i][j][1]=fg.rho[i][j]*fg.F[i][j][3]/fg.F[i][j][1];
            fg.G[i][j][2]=fg.F[i][j][3];
            fg.G[i][j][3]=fg.rho[i][j]*pow((fg.F[i][j][3]/fg.F[i][j][1]),2.0)+fg.F[i][j][2]-pow(fg.F[i][j][1],2.0)/fg.rho[i][j];
            fg.G[i][j][4]=(gamma/(gamma-1))*(fg.F[i][j][2]-pow(fg.F[i][j][1],2.0)/fg.rho[i][j])*(fg.F[i][j][3]/fg.F[i][j][1])+fg.rho[i][j]*0.5*(fg.F[i][j][3]/fg.F[i][j][1])*(pow(fg.F[i][j][1]/fg.rho[i][j],2.0)+pow(fg.F[i][j][3]/fg.F[i][j][1],2.0));
            /*
            fg.G[i][j][1]=fg.rho[i][j]*fg.v[i][j];
            fg.G[i][j][2]=fg.rho[i][j]*fg.u[i][j]*fg.v[i][j];
            fg.G[i][j][3]=fg.rho[i][j]*pow(fg.v[i][j],2.0)+fg.p[i][j];
            fg.G[i][j][4]=(gamma/(gamma-1))*fg.p[i][j]*fg.v[i][j]+fg.rho[i][j]*fg.v[i][j]*(pow(fg.u[i][j],2.0)+pow(fg.v[i][j],2.0))*0.5;
            */
        }
        fg.xi[i]=fg.xi[i-1]+dxi;
        printf("\n%d %f",i,fg.xi[i]);
        printf("\nElement 0");
        j=0;
        //printf("\n%f %d %d %f %f %f %f %f %f %f %f %f %f",fg.xi[i],i,j,fg.u[i][0],fg.v[i][0],fg.rho[i][0],fg.p[i][0],fg.T[i][0],fg.M[i][0],fg.F[i][0][1],fg.F[i][0][2],fg.F[i][0][3],fg.F[i][0][4]);
        //printf("\n%f %d %d %f %f %f %f %f",
        //last element values
        for(k=1;k<=4;k++)
        {
            fg.F[i][ydiv-1][k]=fg.F[i-1][ydiv-1][k];
            fg.G[i][ydiv-1][k]=fg.G[i-1][ydiv-1][k];
        }
        fg.M[i][ydiv-1]=fg.M[i-1][ydiv-1];
        fg.u[i][ydiv-1]=fg.u[i-1][ydiv-1];
        fg.v[i][ydiv-1]=fg.v[i-1][ydiv-1];
        fg.p[i][ydiv-1]=fg.p[i-1][ydiv-1];
        fg.T[i][ydiv-1]=fg.T[i-1][ydiv-1];
        fg.rho[i][ydiv-1]=fg.rho[i-1][ydiv-1];
        if(fg.xi[i]>68.0)
        {
            break;
        }
        j=0;
        if(fg.xi[i]<10.0)
        {
            fi=atan(fg.v[i][j]/fg.u[i][j]);
        }
        else
        {
            psi=atan(fabs(fg.v[i][j])/fg.u[i][j]);
            fi=theta_1-psi;
        }
        //fcal=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)*(pow(fg.M[i][j],2.0)-1)/(gamma+1)))-atan(sqrt(pow(fg.M[i][j],2.0)-1.0));
        fcal=sqrt((gamma-1)*(pow(fg.M[i][j],2.0)-1)/(gamma+1));
        fcal=atan(fcal);
        fcal=sqrt((gamma+1)/(gamma-1))*fcal;
        fcal=fcal-atan(sqrt(pow(fg.M[i][j],2.0)-1.0));
        fact=fcal+fi;
        Mcal=fg.M[i][j];
        pcal=fg.p[i][j];
        Tcal=fg.T[i][j];
        fg.M[i][j]=newtonraphson(fact);
        fg.p[i][j]=pcal*pow((1+(gamma-1)*0.5*Mcal*Mcal)/(1+(gamma-1)*0.5*fg.M[i][j]*fg.M[i][j]),gamma/(gamma-1));
        fg.T[i][j]=Tcal*(1+(gamma-1)*0.5*Mcal*Mcal)/(1+(gamma-1)*0.5*fg.M[i][j]*fg.M[i][j]);
        fg.rho[i][j]=fg.p[i][j]/(287.0*fg.T[i][j]);

        if(fg.xi[i]<10.0)
        {
            fg.v[i][j]=0;
        }
        else
        {
            fg.v[i][j]=-fg.u[i][j]*tan(theta_1);
        }
        fg.G[i][j][1]=fg.rho[i][j]*fg.v[i][j];
        fg.G[i][j][2]=fg.rho[i][j]*fg.u[i][j]*fg.v[i][j];
        fg.G[i][j][3]=fg.rho[i][j]*pow(fg.v[i][j],2.0)+fg.p[i][j];
        fg.G[i][j][4]=(gamma/(gamma-1))*fg.p[i][j]*fg.v[i][j]+fg.rho[i][j]*fg.v[i][j]*(pow(fg.u[i][j],2.0)+pow(fg.v[i][j],2.0))*0.5;

        for(j=0;j<ydiv;j++)
        {
            //printf("\n%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f",i,j,fg.u[i][j],fg.v[i][j],fg.rho[i][j],fg.p[i][j],fg.T[i][j],fg.M[i][j],fg.F[i][j][1],fg.F[i][j][2],fg.F[i][j][3],fg.F[i][j][4],fg.G[i][j][1],fg.G[i][j][2],fg.G[i][j][3],fg.G[i][j][4]);
            printf("\n%d %d %f %f %f %f %f %f",i,j,fg.u[i][j],fg.v[i][j],fg.rho[i][j],fg.p[i][j],fg.T[i][j],fg.M[i][j]);
        }
        dxi=selectdx(1.0,i); //based on 40 grid points

        //Step Operating Paramters
        //fg.dxi[i]=

    }
}
