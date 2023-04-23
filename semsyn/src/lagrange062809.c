/*
   files for storing subroutines for Lagrange interpolants
   Written by Liang Zhao, IGG, CAS.
   Begin time: Nov 27, 2007
   References: Looking for WiKiPedia
*/

#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <stdlib.h>
#include "sem2d.h"

/* function calculating Lagrange function of degree Na for variat at seta
   Na: high-polynomials number 
   ith: the ith interpolant
   gll: 1-D float value array storing Gauss-Lobatto Legendre points
   seta: input value
*/
float Lag_N_i(Na,ith,gll,seta)
int Na, ith;
float *gll;
float seta;
{
    float LNi;
    LNi = 1;
    int j;
    for(j=0; j<=Na; j++){
         if(j != ith)    LNi = LNi* (seta-gll[j])/(gll[ith]-gll[j]);
    }
    return LNi;
}

/* function calculating derivative of Lagrange interpolants of degree Na
   Na: high-polynomials number 
   alpha: the No. alpha interpolant
   gll: 1-D float value array storing Gauss-Lobatto Legendre points
   seta: input value
*/

float dLag_N_i(Na,alpha,gll,seta)
int Na, alpha;
float *gll;
float seta;
{
    int i,k;
    float value,LNi;
    value = 0;
    float down= 1.0;
    for(i=0;i<=Na;i++){
       if(i != alpha) down = down * (gll[alpha]-gll[i]);
    }
//printf("down=%f ",down);

    for(i=0;i<=Na;i++){
       LNi = 1;
       if(i!= alpha){
          for(k=0; k<=Na; k++){
             if(k != alpha && k != i) LNi = LNi* (seta-gll[k]);
          }
          value = value + LNi;
       }
    }
    value = value/down;
//printf(" value=%f\n",value);
    return value;
}

/* function 2 calculating derivative of Lagrange interpolants of degree Na
   Na: high-polynomials number 
   i: the No. i interpolant
   gll: 1-D float value array storing Gauss-Lobatto Legendre points
   k: input reference point
   see Wang TK's PHd thesis, P17
*/

float dLag_N2(Na,i,k,gll)
int Na, i;
float *gll;
int k;
{
    float value,rk,ri;
    float Legengre();
    if(i==0 && k ==0) value = -Na*(Na+1)/4.;
    else if(i==Na && k== Na) value = Na*(Na+1)/4.;
    else if(i==k && i<Na && i>0 )  value = 0.;
    else {
        rk = gll[k]; ri= gll[i];
        value = Legengre(Na,rk)/Legengre(Na,ri)/(rk-ri);
    }
    return value;
}

/*
   Function calculating Legendre polynomial of degree N
   Nl: degree 
   x: input 
*/
float Legengre(Nl,x)
int Nl;
float x;
{
    float pn;
    int i;
    float p0, p1, pn1,pn2;

    if(Nl == 0) { pn = 1; }
    else if(Nl==1) {pn = x;}
    else{ 
       p0= 1; pn2 = p0;
       p1= x; pn1 = p1;
       for(i=2; i<=Nl;i++){
          pn = ((2*i-1)*x*pn1 - (i-1)*pn2)/i;
          pn2 = pn1;
          pn1 = pn;
       }
    }
    return pn;
}

/*
   Function calculating the derivative of the Legendre polynomial of degree N
   Nl: degree 
   x: input 
*/
float d_Legengre(Nl,x)
int Nl;
float x;
{
    float dpn;
    int i;
    float dp0, dp1, dpn1,dpn2;
    float pn1;

    if(Nl == 0) {dpn = 0;}
    else if(Nl == 1) {dpn = 1;}
    else{
       dp0= 0; dpn2 = dp0;
       dp1= 1; dpn1 = dp1;
       for(i=2; i<=Nl;i++){
          pn1 = Legengre(i-1,x);
          dpn = ((2*i-1)*(x*dpn1 + pn1) - (i-1)*dpn2)/i;
          dpn2 = dpn1;
          dpn1 = dpn;
       }
    }
    return dpn;
}

/*
  Functions calculating the shape functions
  input: rx,rz: reference coordinate
         ith: the ith of 4 shape functions
  output: value of the ith shape function
  reference: see my .doc file for equations : shape functions
*/
float shape_4(rx,rz,ith)
float rx,rz;
int ith;
{
    float nv; 
    if(ith == 1) { nv = (1-rx-rz+rx*rz)/4; }
    if(ith == 2) { nv = (1+rx-rz-rx*rz)/4; }
    if(ith == 3) { nv = (1-rx+rz-rx*rz)/4; }
    if(ith == 4) { nv = (1+rx+rz+rx*rz)/4; }
    return nv;
}

/*
   subroutine for calculating the Jacobian and its inverse Matrix of degree 1 for 4 control anchors
   input: rx,rz: reference coordinate
          xa,za, a= 1,2,3,4, x and z of the four control points
   output: Je[2*2], Jacobian matrix 2*2
   reference: see my .doc file for equations : shape functions
*/
float Jacobian_4(rx,rz,xa,inv_m)
float rx,rz;
float inv_m[4];
struct Point2d xa[4];
{
    int i;
    float Nrx[4], Nrz[4]; /* the derivative of shape function WRT rx or rz by analytically inducing*/

    struct Je4_one Je4;
    float jev;

    Nrx[0] = -1 + rz; Nrz[0] = -1 + rx;
    Nrx[1] =  1 - rz; Nrz[1] = -1 - rx;
    Nrx[2] = -1 - rz; Nrz[2] =  1 - rx;
    Nrx[3] =  1 + rz; Nrz[3] =  1 + rx;

    for(i=0;i<4;i++) Je4.je[i] = 0;
    for(i=0;i<4;i++){
       Je4.je[0] += xa[i].x*Nrx[i]/4.;
       Je4.je[1] += xa[i].x*Nrz[i]/4.;
       Je4.je[2] += xa[i].z*Nrx[i]/4.;
       Je4.je[3] += xa[i].z*Nrz[i]/4.;
    }
    jev= Je4.je[0]*Je4.je[3] - Je4.je[1]*Je4.je[2];
//    printf("%f %f %f %f je=%f\n", Je4.je[0], Je4.je[1],Je4.je[2],Je4.je[3],jev);
    /*Get inverse matrix*/
    inv_m[0] =  Je4.je[3]/jev;
    inv_m[1] = -Je4.je[2]/jev;
    inv_m[2] = -Je4.je[1]/jev;
    inv_m[3] =  Je4.je[0]/jev;
//    printf("%f\n", inv_m[0]);
    return jev;
}

/*
Weights of GLL integration (Abrmowttz, 1984)
   Input: theta: the value of coordinate of the GLL points, -1<= theta <= 1
              N: degree of the Legendre polynomials
   Output: weight of the GLL integrants
*/
void GLLWeight(gll, np, gllw)
float * gll, *gllw;
int np;
{
   float gglv,theta,tempf;
   int i;
   
   for(i=0;i<=np;i++){
      theta = gll[i];
      if(theta>1 || theta< -1) {
         printf("theta is out of range (-1,1) in GLLWeitht(), lagrange.h\n");
         exit(-1);
      }
      if(i ==0 || i == np){
         gglv= 2./(np*(np+1));
      }
      else {
         tempf = Legengre(np,theta);
         gglv= 2./(np*(np+1)*tempf*tempf);
      }
      gllw[i] = gglv;
   }
}

/*
  Gll5() Calculating GLL points for degree = 5
  For Pn'(x) np=5
  get the n+ 1 roots of : ax4+cx2+e = 0
  np: degree of Polynomials,
*/
void GLL5(np,gllp)
int np;
float gllp[6];
{
    if(np != 5) {
       printf("Error invoke: degree is not equal to 5 !, in GLL5() \n");
       exit(-1);
    }
    float a,c,e;
    a = 1575./40.; c = -1050./40.; e= 85./40.;
 
    gllp[0] = -1.;
    gllp[1] = -sqrt((-c + sqrt(c*c-4*a*e))/2/a);
    gllp[2] = -sqrt((-c - sqrt(c*c-4*a*e))/2/a);
    gllp[3] =  sqrt((-c - sqrt(c*c-4*a*e))/2/a);
    gllp[4] =  sqrt((-c + sqrt(c*c-4*a*e))/2/a);
    gllp[5] = 1.;
}

/*
  References: Wiki Pedia: Legendre polynomias
-------------------
  n	Pn(x)
  0	1
  1	x
  2     0.5*(3x**2-1)
  3     1/2*(5x**3-3*x)
  4     1/8*(35x**4-30x**2+3)
  5     1/8*(63x**5-70x**3+15x)
  6     1/16*(231x**6-315x**4+105x**2-5)
  7     1/16*(429x**7-693x**5+315x**3-35x)
  8	1/128*(6345x**8-12012x**6+6930x**4-1260x**2+35)
-------------------
  Gll() Calculating GLL points for degree <=5
  for degree 4: (35x*x-15)x(1-x*x)/2=0
             3: 3(5x*x-1)(1-x*x)/2=0
  np: degree of Polynomials,
*/
void GLL(np,gllp)
int np;
float * gllp;
{
    if(np > 5 || np<3) {
       printf("Error invoke: degree is out of 3,4 or 5 !, in GLL() \n");
       exit(-1);
    }

    if(np==3){
        gllp[0] = -1.;
        gllp[1] = -sqrt(1/5.);
        gllp[2] =  sqrt(1/5.);
        gllp[3] = 1.;
    }
    if(np==4){
        gllp[0] = -1.;
        gllp[1] = -sqrt(15./35.);
        gllp[2] = 0.;
        gllp[3] = sqrt(15./35.);
        gllp[4] = 1.;
    }
    if(np==5){
        float a,c,e;
        a = 1575./40.; c = -1050./40.; e= 75./40.;
        gllp[0] = -1.;
        gllp[1] = -sqrt((-c + sqrt(c*c-4*a*e))/2/a);
        gllp[2] = -sqrt((-c - sqrt(c*c-4*a*e))/2/a);
        gllp[3] =  sqrt((-c - sqrt(c*c-4*a*e))/2/a);
        gllp[4] =  sqrt((-c + sqrt(c*c-4*a*e))/2/a);
        gllp[5] = 1.;
    }
}

/*
   Calculate Ricker wavelet time-sequences
   npts: time step number
   rp: time sequence
   f0: main frequence
   dt: time step
   A0: amplitude
   t0: time delay
*/
void Ricker_Wavelet(rp,npts,dt,f0,t0,a0)
int npts;
float dt,f0,t0,a0;
float * rp;
{
    int i, nt0;
    float pi_f0, t_t0;
    float pi = 3.14159265;

    pi_f0 = (pi/f0) * (pi/f0);
    nt0 = (int) t0/dt;
    printf("f0=%f pi_f0=%f nt0=%d\n",f0, pi_f0,nt0);
    for(i=0;i<npts;i++){
       t_t0 = (i - nt0) * (i - nt0) * dt *dt;
       rp[i] = a0*(1 - 2*pi_f0*t_t0) * exp(-pi_f0*t_t0);
    }
}
