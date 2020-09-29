#include<stdio.h>
#include<math.h>
#include<conio.h>
#define L 0.5
#define H 0.01
#define nx 21
#define ny 11
int main()
{
    double p[nx][ny], u[nx+1][ny], v[nx+2][ny+1];
    double ps[nx][ny], us[nx+1][ny], vs[nx+2][ny+1];
    double pd[nx][ny], ud[nx+1][ny], vd[nx+2][ny+1];
    double dx,dy,dt;
    double a,b,c,d[nx][ny];
    double ue,mu,vbar,v2bar,ubar,u2bar,rho,w,err,err1;
    int i,j,loci,locj,it,ito,exit,temp;
    double As,Bs;
    printf("\nHello");
    rho=0.002377;
    ue=1;
    mu=3.737E-7;
    //Operating parameters
    dx=L/(nx-1);
    dy=H/(ny-1);
    dt=0.001;
    a=2.0*((dt/(dx*dx))+(dt/(dy*dy)));
    b=-dt/(dx*dx);
    c=-dt/(dy*dy);
    //Initial Boundary Conditions
    //Top wall
    for(i=0;i<=nx;i++)
    {
        u[i][ny-1]=ue;
    }
    for(i=0;i<nx;i++)
    {
        pd[i][ny-1]=0;
    }
    //Bottom Wall
    for(i=0;i<=nx;i++)
    {
        u[i][0]=0;
    }
    for(i=0;i<nx;i++)
    {
        pd[i][0]=0;
    }
    //Left Wall
    for(j=0;j<=ny;j++)
    {
        v[1][j]=0;
    }
    for(j=0;j<ny;j++)
    {
        pd[0][j]=0;
    }
    //Right Wall
    for(j=0;j<ny;j++)
    {
        pd[nx-1][j]=0;
    }
    //Initial Conditions
    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            p[i][j]=0;
			pd[i][j]=0;
        }
    }
    for(i=0;i<nx+1;i++)
    {
        for(j=0;j<ny;j++)
        {
            u[i][j]=0;
			us[i][j]=0;
        }
    }
    for(i=0;i<nx+2;i++)
    {
        for(j=0;j<ny+1;j++)
        {
            v[i][j]=0;
			vs[i][j]=0;
        }
    }
    //Top wall
    for(i=0;i<=nx;i++)
    {
        u[i][ny]=ue;
    }
    v[14][4]=ue/2.0;
	vs[14][4]=ue/2.0;
    //Loop
    ito=0;
	exit=0;
	temp=1;
    while(temp==1)
    {

        //Assign Pstar?
        for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                ps[i][j]=p[i][j];
            }
        }
		for(i=0;i<nx+1;i++)
		{
			for(j=0;j<ny;j++)
			{
				us[i][j]=u[i][j];
			}
		}
		for(i=0;i<nx+2;i++)
		{
			for(j=0;j<ny+1;j++)
			{
				vs[i][j]=v[i][j];
			}
		}
        //Calculate ustar
        for(i=1;i<nx;i++)
        {
            for(j=1;j<ny-1;j++)
            {
                vbar=0.5*(v[i][j+1]+v[i+1][j+1]);
                v2bar=0.5*(v[i][j]+v[i+1][j]);
                As=-((rho*u[i+1][j]*u[i+1][j]-rho*u[i-1][j]*u[i-1][j])/(2.0*dx)+(rho*u[i][j+1]*vbar-rho*u[i][j-1]*v2bar)/(2.0*dy))+mu*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/(dy*dy));
                us[i][j]=(rho*us[i][j]+As*dt-(dt/dx)*(ps[i][j]-ps[i-1][j]))/rho;
            }
        }
        //Calculate vstar
        for(i=1;i<nx+1;i++)
        {
            for(j=1;j<ny;j++)
            {
                ubar=0.5*(u[i][j-1]+u[i][j]);
                u2bar=0.5*(u[i-1][j-1]+u[i-1][j]);
                Bs=-((rho*v[i+1][j]*ubar-rho*v[i-1][j]*u2bar)/(2.0*dx)+(rho*v[i][j+1]*v[i][j+1]-rho*v[i][j-1]*v[i][j-1])/(2.0*dy))+mu*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy));
				vs[i][j]=(rho*vs[i][j]+Bs*dt-(dt/dy)*(ps[i-1][j]-ps[i-1][j-1]))/rho;
				if(i==14&&j==4)
				{
					//printf("\nubar %f",ubar);
					//printf("\nu2bar %f",u2bar);
					//printf("\nv[%d][%d] %f",i+1,j,v[i+1][j]);
					//printf("\nv[%d][%d] %f",i-1,j,v[i-1][j]);
					//printf("\nv[%d][%d] %f",i,j+1,v[i][j+1]);
					//printf("\nv[%d][%d] %f",i,j-1,v[i][j-1]);
					//printf("\nv[%d][%d] %f",i,j,v[i][j]);
					//printf("\nps[%d][%d] %f",i-1,j+1,ps[i-1][j+1]);
					//printf("\nps[%d][%d] %f",i-1,j,ps[i-1][j]);
					//printf("\nBs %f",Bs);
					//printf("\nBs part 1 %f",-((rho*v[i+1][j]*ubar-rho*v[i-1][j]*u2bar)/(2.0*dx)+(rho*v[i][j+1]*v[i][j+1]-rho*v[i][j-1]*v[i][j-1])/(2.0*dy)));
					//printf("\nBs part 2 %f",mu*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy)));
					//printf("\nvs[%d][%d] %f",i,j,vs[i][j]);
				}
            }
        }
        //Pressure Correction Jacobian
        it=0;

		for(i=0;i<nx;i++)
		{
			for(j=0;j<ny;j++)
			{
				pd[i][j]=0;
			}
		}

        while(1)
        {
            w=1.5;
			err=0;
            for(i=1;i<nx-1;i++)
            {
                for(j=1;j<ny-1;j++)
                {
					err1=pd[i][j];
					//printf("\nOld pd[%d][%d] %f",i,j,pd[i][j]);
					//printf("\nus[%d][%d] %f",i+1,j,us[i+1][j]);
                    //printf("\nus[%d][%d] %f",i,j,us[i][j]);
                    //printf("\nvs[%d][%d] %f",i+1,j+1,vs[i+1][j+1]);
                    //printf("\nvs[%d][%d] %f",i+1,j,vs[i+1][j]);
					//printf("\na %f b %f c %f",a,b,c);
					d[i][j]=(rho/dx)*(us[i+1][j]-us[i][j])+(rho/dy)*(vs[i+1][j+1]-vs[i+1][j]);
					pd[i][j]=(1-w)*pd[i][j]+w*(-(1.0/a)*(b*pd[i+1][j]+b*pd[i-1][j]+c*pd[i][j+1]+c*pd[i][j-1]+d[i][j]));
					err1=err1-pd[i][j];
					//printf("\nNew pd[%d][%d] %f",i,j,pd[i][j]);
					if(fabs(err1)>err)
                    {
                        err=fabs(err1);
						loci=i;
                        locj=j;
					}
				}
            }
			/*
            for(i=1;i<nx-1;i++)
            {
                for(j=1;j<ny-1;j++)
                {
                    if(err<fabs(pd[i][j]-(1.0/a)*(b*pd[i+1][j]+b*pd[i-1][j]+c*pd[i][j+1]+c*pd[i][j-1]+d[i][j])))
                    {
                        loci=i;
                        locj=j;
                        err=fabs(pd[i][j]-(1.0/a)*(b*pd[i+1][j]+b*pd[i-1][j]+c*pd[i][j+1]+c*pd[i][j-1]+d[i][j]));

                    }
					//printf("\n i %d j %d err %f",i,j,pd[i][j]-(1.0/a)*(b*pd[i+1][j]+b*pd[i-1][j]+c*pd[i][j+1]+c*pd[i][j-1]+d[i][j]));
                }
            }
			*/
            printf("\nIteration number is %d loci is %d locj is %d Error is %.8f",it,loci,locj,err);
            it++;
            if(err<1e-8)
            {
                break;
            }

        }
        //Pressure Change
        for(i=1;i<nx-1;i++)
        {
            for(j=1;j<ny-1;j++)
            {
                p[i][j]=ps[i][j]+0.1*pd[i][j];
            }
        }
        //U Velocity Correction
        for(i=1;i<nx;i++)
        {
            for(j=1;j<ny-1;j++)
            {
                ud[i][j]=-(dt/(rho*dx))*(pd[i][j]-pd[i-1][j]);
                u[i][j]=us[i][j]+ud[i][j];
            }
        }
        //V Velocity Correction
        for(i=2;i<nx+1;i++)
        {
            for(j=1;j<ny;j++)
            {
                vd[i][j]=-(dt/(rho*dy))*(pd[i-1][j]-pd[i-1][j-1]);
                v[i][j]=vs[i][j]+vd[i][j];
            }
        }
        //Some Boundary Conditions
        //v bottom & top
        for(i=1;i<nx+1;i++)
        {
            v[i][0]=-v[i][1];
            v[i][ny]=-v[i][ny-1];
        }
        //v left & right
        for(j=0;j<ny+1;j++)
        {
            v[0][j]=0;
            v[1][j]=0;
            v[nx+1][j]=v[nx][j];
        }
        //u bottom & top
        for(i=0;i<nx+1;i++)
        {
            u[i][0]=0;
            u[i][ny-1]=ue;
        }
        //u right and left
        for(j=0;j<ny;j++)
        {
            u[0][j]=u[1][j];
            u[nx][j]=u[nx-1][j];
        }
        printf("\nIteration number [outer loop] is %d, v[14][4] is %.8f",ito,v[14][4]);
		temp=0;
		for(i=2;i<nx+1;i++)
        {
            for(j=1;j<ny;j++)
            {
                if(fabs(v[i][j])>1e-6)
				{
					temp=1;
					continue;
				}
            }

        }
		if(temp==0)
		{
			i=14;
			for(j=0;j<ny;j++)
			{
				printf("\nu[%d][%d] is %f",i,j,u[i][j]);
			}

		}
        ito++;


    }
}
