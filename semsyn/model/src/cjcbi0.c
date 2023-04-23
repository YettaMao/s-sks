/*
  Solve for eigen value by Jacobi method
  Refer to Xu Shiliang, P74
*/
#include  <stdio.h>
#include  <stdlib.h>
#include  <fcntl.h>
#include  <math.h>
#include  <malloc.h>

main(int argc,char **argv)
{ 
    int cjcbi();
    void inv_index();

    int axisnum; /* in which direction the plane wave propagate*/
    float density; 
    axisnum = 3;
    char elafile[80];  /* file for storing the initial elastic constant tensor */

    int i,j;
    int op,qr;

    printf("argc=%d\n",argc);
    for(i=0;i<argc;i++){
       printf("%s ",argv[i]);
    }
    printf("\n");
    if(argc != 4) 
    {
       printf("Elastic system rotation \n\n");
       printf("Usage: elastic_rotate elafile propagating_axis_num density\n");
    }
    sscanf(&argv[1][0], "%s",elafile);
    sscanf(&argv[2][0], "%d",&axisnum);
    sscanf(&argv[3][0], "%f",&density);
    axisnum --;

    FILE * fp, *fop;
    float c_opqr[6][6]; /* storing initial matrix of elastic constants */
    float c_ij[3][3]; /* T matrix */

    if ((fp = fopen( elafile, "rt")) == NULL){
       printf( "cannot open this file\n");
       exit(0);
    } 
    printf(" ---------- Old system -------------\n");
    for(op=0;op<6;op++){
       for(qr=0;qr<6;qr++){
           fscanf(fp,"%f",&c_opqr[op][qr]);
           printf("%-8.2f",c_opqr[op][qr]);
       }
       printf("\n");
    }
    fclose(fp);

/*    if ((fop = fopen( rotfile, "wt")) == NULL){
       printf( "cannot open this file\n");
       exit(0);
    } 
    printf("%-8.2f",c_jkmn[jk][mn]);
    fprintf(fop,"%-8.2f",c_jkmn[jk][mn]);
    fclose(fop);
*/
    
    float velocity[3] ;
    printf("Plane wave Propagating along axis-%d, its T matrix is\n",axisnum+1);
    for(i=0;i<3;i++){
       inv_index(i,axisnum,&op);
       for(j=0;j<3;j++){
           inv_index(j,axisnum,&qr);
           printf("(%d %d) ", op+1,qr+1);
           c_ij[i][j] = c_opqr[op][qr];
       }
       printf("\n");
    }

    float eps, aniso;
    float v[3][3];
    eps=0.000001;
    i=cjcbi(c_ij,3,v,eps,100);
    printf("velocity is:\n");
    if (i>0) { 
       for (i=0; i<=2; i++){
          velocity[i] = sqrt(c_ij[i][i]/density);
          printf("%-8.3f ",velocity[i]);
       }
       printf("\n\n");
       for (i=0; i<=2; i++){ 
            for (j=0; j<=2; j++) printf("%13.7e  ",v[i][j]);
              printf("\n");
       }
       printf("\n");
    }
    aniso = fabs((velocity[0]-velocity[1])/(velocity[0]+velocity[1])*2*100);
    printf("Shear wave anisotropy percent= %6.2f %\n\n", aniso);

}

/* Musgrave index */
/* in c array
0  1  2  3  4  5
00 11 22 12 02 01
or 
1  2  3  4  5  6
11 22 33 23 13 12
*/
void inv_index(i,j,op)
int i,j;
int * op;
{
   if(i==0 && j==0) *op = 0;
   if(i==0 && j==1) *op = 5;
   if(i==0 && j==2) *op = 4;
   if(i==1 && j==0) *op = 5;
   if(i==1 && j==1) *op = 1;
   if(i==1 && j==2) *op = 3;
   if(i==2 && j==0) *op = 4;
   if(i==2 && j==1) *op = 3;
   if(i==2 && j==2) *op = 2;
}

  int cjcbi(a,n,v,eps,jt)
  int n,jt;
  float a[],v[],eps;
  { 
    int i,j,p,q,u,w,t,s,l;
    float fm,cn,sn,omega,x,y,d;
    l=1;
    for (i=0; i<=n-1; i++){ 
        v[i*n+i]=1.0;
        for (j=0; j<=n-1; j++)
          if (i!=j) v[i*n+j]=0.0;
    }
    while (1==1)
      { fm=0.0;
        for (i=1; i<=n-1; i++)
        for (j=0; j<=i-1; j++)
          { d=fabs(a[i*n+j]);
            if ((i!=j)&&(d>fm))
              { fm=d; p=i; q=j;}
          }
        if (fm<eps)  return(1);
        if (l>jt)  return(-1);
        l=l+1;
        u=p*n+q; w=p*n+p; t=q*n+p; s=q*n+q;
        x=-a[u]; y=(a[s]-a[w])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0) omega=-omega;
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a[w];
        a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
        a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
        a[u]=0.0; a[t]=0.0;
        for (j=0; j<=n-1; j++)
        if ((j!=p)&&(j!=q))
          { u=p*n+j; w=q*n+j;
            fm=a[u];
            a[u]=fm*cn+a[w]*sn;
            a[w]=-fm*sn+a[w]*cn;
          }
        for (i=0; i<=n-1; i++)
          if ((i!=p)&&(i!=q))
            { u=i*n+p; w=i*n+q;
              fm=a[u];
              a[u]=fm*cn+a[w]*sn;
              a[w]=-fm*sn+a[w]*cn;
            }
        for (i=0; i<=n-1; i++)
          { u=i*n+p; w=i*n+q;
            fm=v[u];
            v[u]=fm*cn+v[w]*sn;
            v[w]=-fm*sn+v[w]*cn;
          }
      }
    return(1);
  }

