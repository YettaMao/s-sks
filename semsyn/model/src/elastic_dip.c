/*************************************************************************/
/*   Program for dipping elastic tensor                                  */
/*   Written by: Liang Zhao, IGG, CAS                                    */
/*   Based on: elastic_rotate.c                                          */
/*   Last modified: Feb 2, 2009                                          */
/*************************************************************************/

/* rotating elastic tensor in x1-x2 plane
   Usage:  ./elastic_rotate par=*.par
   Input:  get par from *.par
           *.ela
   Output: new.ela 
*/
/*
Refer to: Crampin, S.,1977, A review of the effects of anisotropic layering on the propagation of seismic waves, Geophys. J. R. astr. Soc., 49, 181-208.
Cjkmn=g(j,o)g(k,p)g(m,q)g(n,r)Copqr
here g(j,o) represents a partial derivative of xj wrt xo
*/
#include  <stdio.h>
#include  <stdlib.h>
#include  <sys/file.h>
#include  <fcntl.h>
#include  <math.h>
#include  <malloc.h>

main(int argc,char **argv)
{
   void get_index();
   float rotate();

   float angle; /* rotating angle, unit in degree*/
   char elafile[80];  /* file for storing the initial elastic constant tensor */
   char rotfile[80];  /* file for storing elastic constant tensor after rotation */

   int ijk;
   printf("argc=%d\n",argc);
   for(ijk=0;ijk<argc;ijk++){
      printf("%s ",argv[ijk]);
   }
   printf("\n");
   if(argc != 4) 
   {
      printf("Elastic system rotation \n\n");
      printf("Usage: elastic_rotate elafile rotfile angle \n");
   }
   sscanf(&argv[1][0], "%s",elafile);
   sscanf(&argv[2][0], "%s",rotfile);
   sscanf(&argv[3][0], "%f",&angle);

//   setpar(ac,av);
//   mstpar("elafile",     "s", elafile); /* file for storing the initial elastic constant tensor */
//   mstpar("rotfile",     "s", rotfile); /* file for storing elastic constant tensor after rotation*/
//   mstpar("angle",       "f", &angle); /* rotating angle, unit in degree*/
//   endpar();   

   FILE * fp, *fop;
   float c_opqr[6][6]; /* storing initial matrix of elastic constants */
   float c_jkmn[6][6]; /* matrix of elastic constants after rotation */

   int op,qr;
   int jk,mn;
   
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

   int j, k, m, n;

   printf(" ---------- New system -------------\n"); 
   if ((fop = fopen( rotfile, "wt")) == NULL){
      printf( "cannot open this file\n");
      exit(0);
   } 

   for(jk=0;jk<6;jk++){
      get_index(jk,&j,&k);
      for(mn=0;mn<6;mn++){
          get_index(mn, &m, &n);
//          printf("j=%-3d k=%-3d m=%-3d n=%-3d\n",j,k,m,n);
          c_jkmn[jk][mn] = rotate(j,k,m,n,c_opqr,angle);
          printf("%-8.2f",c_jkmn[jk][mn]);
          fprintf(fop,"%-8.2f",c_jkmn[jk][mn]);
      }
      printf("\n");
      fprintf(fop,"\n");
   }
   fclose(fop);
}

/* Musgrave index */
/* in c array
0  1  2  3  4  5
00 11 22 12 02 01
or 
1  2  3  4  5  6
11 22 33 23 13 12
*/

void get_index(op,i,j)
int op;
int * i, * j;
{
   if(op == 0) {*i = 0; *j = 0;}
   if(op == 1) {*i = 1; *j = 1;}
   if(op == 2) {*i = 2; *j = 2;}
   if(op == 3) {*i = 1; *j = 2;}
   if(op == 4) {*i = 0; *j = 2;}
   if(op == 5) {*i = 0; *j = 1;}
}

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

float axis_rot(j,o,angle)
int j,o; /* j: new ; o: old */
float angle; /* in degree, angle between new x1 - old x1, counter-clockwise */
{
   float axis, ang;
   float pi = 3.1415926536;

   if( j==0 && o == 0) ang = angle;
   if( j==0 && o == 1) ang = angle - 90;
   if( j==0 && o == 2) ang = 90;
   if( j==1 && o == 0) ang = angle + 90;
   if( j==1 && o == 1) ang = angle;
   if( j==1 && o == 2) ang = 90;
   if( j==2 && o == 0) ang = 90;
   if( j==2 && o == 1) ang = 90;
   if( j==2 && o == 2) ang = 0;

   ang = ang/180 * pi;
   axis = cos(ang);
   return axis;
}

/*axis_dip within x1-x3 plane, x3: vertical*/
float axis_dip(j,o,angle)
int j,o; /* j: new ; o: old */
float angle; /* in degree, angle between new x1 - old x1, counter-clockwise */
{
   float axis, ang;
   float pi = 3.1415926536;

   if( j==0 && o == 0) ang = angle;
   if( j==0 && o == 1) ang = 90;
   if( j==0 && o == 2) ang = angle - 90;
   if( j==1 && o == 0) ang = 90;
   if( j==1 && o == 1) ang = 0;
   if( j==1 && o == 2) ang = 90;
   if( j==2 && o == 0) ang = angle + 90;
   if( j==2 && o == 1) ang = 90;
   if( j==2 && o == 2) ang = angle;

   ang = ang/180 * pi;
   axis = cos(ang);
   return axis;
}

/*new x1 rotates wrt x1*/
float rotate(j,k,m,n,c_opqr,angle)
int j,k,m,n;
float angle;
float c_opqr[6][6];
{
    int o,p,q,r;
    int op,qr;
    float jo, kp, mq,nr;
    float value = 0.;
    void inv_index();
    float axis_rot();
    float temp;

    for(o=0; o< 3; o++){
       jo = axis_dip(j,o,angle);
       for(p=0; p< 3; p ++){
          kp = axis_dip(k,p,angle);
          inv_index(o,p,&op);
          for(q=0;q<3; q++){
             mq = axis_dip(m,q,angle);
             for(r=0;r<3;r++){
                 nr = axis_dip(n,r,angle);
                 inv_index(q,r,&qr);
                 temp = jo * kp * mq * nr * c_opqr[op][qr];
                 value += temp;
//if(j==0 && k==0 && m == 0 && n ==0) printf("%-3d %-3d %-3d %-3d value=%-8.2f\n",o,p,q,r,temp); 
             }
          }
       }
    }
    return value;
}
