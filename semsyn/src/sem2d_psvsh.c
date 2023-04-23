
/*********************************************************************/
/*                                                                   */
/*          Codes for Spectral-element method                        */
/*          2-D version                                              */
/*  Written by: Liang Zhao, Institute of Geology and Geophysics, CAS */
/*  Begin time: Thusday, Nov. 22, 2007.                              */
/*                                                                   */
/*********************************************************************/

/*
  In this code, Element means a element unit, while Grid means a interpolation point of mesh. 
  Each element contains (hpn+1)*(hpn+1) grid points.
  Cartesian (x,y,z) with   x = east,
                           y = north,
                           z = upward;
  Note that the z direction here is different from those in GRT calculations, see GRT_interface()
*/

/* modified history:
  Dec 10, 07: storing c[i][j][k][l] instead of calculating in each step to save time, 
              but at t he expense of increasing memory.
  Dec 17, 07: storing Fik[alpha][beta] instead of caculating in each step to save time
  Dec 27, 07: modifying source term: the source sharing by four adjacent elements
              adding attenuation codes
  Jan 11, 08: Applying absorbing boundary conditions used in FD 
              Absorbs2() instead of Absorbs()
  Jan 15, 08: Absorbs() is ok now, it seems that something is wrong with the equation of absorbing boundary 
  May 25, 11: Begin to reduce the high-frequency noise in the waveform  
              Note that the signs of vertical component in GRT and SEM are different 
  July 1, 11: V1.0, finished 
  20140327: using openMP    
  20141028: extended for P-SV & SH imping at the same time      
*/

#include   <stdio.h>
#include   <stdlib.h>
#include   <sys/file.h>
#include   <fcntl.h>
#include   <math.h>
#include   <malloc.h>
#include   <omp.h>
#include   "nrutil.h"
#include   "sem2d.h"

/* Some parameters for checking */
int  check = 1;    /*check=1: output file for debug*/
/* global variants */
int   nex, nez;  /* element number horizontal and vertical */
int   nx, nz;    /* grid points number horizontal and vertical */
int   neh = 1;         /* number of elements for absorbing and attenuation*/
int   nel;       /* element number */
float eh ;       /* spatial step of the spectral element*/
float h;         /* grid point spacing = eh/hpn */
int   nfil = 1;  /* the layer number where SEM ends */
int   hpn = 5;   /* high-ploynomials-degree of Lagrange interpolants, number= hpn + 1*/
float Uni_Je,Uni_JeI; /*Je and its inverse matrix for each grid point*/
                      /*If partition uniformly, Je and JeI is very simple*/
float Uni_Jb,Uni_JbI; /*Jb and its inverse matrix for each grid point for boundary element*/
                      /*If partition uniformly, Jb and JbI is very simple*/
float * JeArray; /*If partition not regularly*/
float * gllp, *gllw; /* GLL points and their weights for integration*/
float * dlagP, *lagP; /*value of Lagrange function at the [-1,-1][1,1] reference grid points*/
float * dlagP_JeI, *dlagP_Gllw ; /*value of Lagrange function at the [-1,-1][1,1] reference grid points * Uni_JeI */

float * grtuw;  /*grt excitation*/
float * grtuw0; /*grt excitation at previous time step*/
float * grtv;  /*grt SH excitation*/
float * grtv0; /*grt SH excitation at previous time step*/

/* filenames */
char greenfile[80];  /* file for GRT output */
char shgreenfile[80];  /* file for GRT SH output */
char semmodel[80];   /* file for model of sem region */
char raymodel[80];   /* file for model in GRT calculation */

/*input parameters*/
int  readpars = 1; /*readpar = 1: input model parameters*/
float xmin;    /*xmin: the horizontal epicentral distance of zero point of SEM */
int   output_vd = 0; /*in vel.xy, output c44-c55*/

main(ac,av)
int ac; 
char **av;
{
   int   i, j ,len1;
   long int len;
   int   flat =1;       /* flag for earth: 1-> flat layers 0->spherical */
   FILE *infile, *shinfile, *open_file();

/* ----- Part 1: definitions and preparation ----- */
/* elemental discretization */

   int uniform=1; /*partition uniformly or not: =1, uniformly; = 0, unevenly */
   int nh=6;     /*number of meshes for absorbing region*/
   int atten =0;  /**/
   float ass=0.93;/*attenuation coefficient*/
   float* taper;  /*for attenuation*/
   int absorb=1;  /*1: using absorbing boundary condition else 0: not */
/* GRT-SEM interface description*/
   int mtotal;    /* GRT-SEM interfaces grid */
   long int rec_len;
   float **sy, **sy_sh, **matrix(); /*For GRT-SEM interfaces*/
   int * ns;
   float rdt; /* = dst/dp */
   void interpl();
/* source definition*/
   int msource;  /*=1: using moment source; else using GRT excitation at the left and bottom interface*/
   int mex, mez; /* source position at the (mex,mez) element */
   float a0, f0, t0;/*Ricker wavelet: amplitude, major frequency and delaytime*/
   float *Ricker;/*Ricker wavelet*/
/* output parameters */
   int   get_rec=0; /*if using receiver position input*/
   int   plot_trace=0, itrace;
   int   xtrace=20,ztrace=1, idxplot=20,idzplot=1;
   int   ixbeg,izbeg;
   float xbeg=0.,zbeg=0.0;
   FILE  *fopu1, *fopv1, *fopw1, *rfp;
   char  rec_file[80]; /* assgin receiver's position*/
   float rec_x;
   int   ixrec,ix1[500],iz1[500]; 

/* time distretization */
   float dp;  /*time step for GRT*/
   int it,nt; /*time step counter; number of time steps for GRT*/
   int nst;   /*number of time steps for SEM*/
   float sdt; /*time step for SEM, second */

/* subroutine definition, all subroutine need to be declare here before use*/
   int   discret(); /* discretization */
   void  Mass_Matrix(); /* Mass_Matrix for all the elements*/
   void  Absorb_Matrix(); /*absorbing boundary for the boundary elements*/
   float Lag_N_i(), GLL5(), GLL(), GLLWeight(), dLag_N2(),dLag_N_i();
   void  atten_taper(), atten_wave(), atten_wave2();/*attenuation*/
   void  All_Je(); /*Jacobian matrix array for all the elements */
   void  Initial_uvw();
   void  Copy_Displacement();
   void  Stiff_Matrix();
   void  Time_diff(), Time_diff_hybrid();
   void  GRT_Interface(); /*GRT excitation*/
   void  GRT_Interface2();/*GRT excitation*/
   void  GRT_Interface3();/*GRT excitation*/
   void  GRT_Interface4();/*GRT excitation*/
   void  Explosive(); /*Using explosive source*/
   void  Ricker_Wavelet();
   void  kirrecord_zhl();/* For output waveform at the receivers */

/* For check only */
   int ist=100, incre=50;
   int plot_snap= 0 ;
   void snapshot();

/* assign parameters from the parameter file using getpar tool */
   setpar(ac,av);
   getpar("flat",          "d", &flat);
   mstpar("greenfile",     "s", greenfile);
   mstpar("shgreenfile",   "s", shgreenfile);
   mstpar("raymodel",      "s", raymodel);
   mstpar("semmodel",      "s", semmodel);

   mstpar("nt",            "d", &nt);
   mstpar("dp",            "f", &dp);

   mstpar("nex",           "d", &nex);
   mstpar("eh",            "f", &eh);
   mstpar("xmin",          "f", &xmin);
   mstpar("nst",           "d", &nst);
   mstpar("sdt",           "f", &sdt);
   getpar("nfil",          "d", &nfil);
   getpar("uniform",       "d", &uniform);
   uniform = 1; /*now, we can only deal with uniform grid, ^_^, update come soon*/
   getpar("hpn",           "d", &hpn);
   getpar("absorb",        "d", &absorb);
   getpar("readpars",      "d", &readpars);
   getpar("neh",           "d", &neh);
   nh = neh * hpn ;
   getpar("atten",         "d", &atten);
   getpar("ass",           "f", &ass);
   getpar("check",         "d", &check);

   getpar("output_vd",     "d", &output_vd);
   getpar("xtrace",        "d", &xtrace);
   mstpar("xbeg",          "f", &xbeg);
   getpar("zbeg",          "f", &zbeg);
   getpar("get_rec",       "d", &get_rec);
   getpar("plot_trace",    "d", &plot_trace);
   if(get_rec == 1){
      mstpar("rec_file",   "s", rec_file);
   }
   mstpar("idxplot",       "d", &idxplot);
 
   mstpar("msource",	   "d", &msource);
   if(msource == 1){
       mstpar("mex",	   "d", &mex);
       mstpar("mez",	   "d", &mez);
       mstpar("a0",        "f", &a0);
       mstpar("f0",	   "f", &f0);
       mstpar("t0",	   "f", &t0);
   }

/*For test*/
   getpar("startsnap",     "d", &ist);
   getpar("incresnap",     "d", &incre);
   getpar("plot_snap",     "d", &plot_snap);
   endpar();

   rec_len = (long) (nt*sizeof(float));

/*Note: hpn is fixed to be 3-6 now, since we only use analytical determination of GLL points for degree of 3-6 */

    /*Part 1.2 Get GLL points and its weights for integration*/
   gllp = (float *) malloc((hpn+1)*sizeof(float));
   gllw = (float *) malloc((hpn+1)*sizeof(float));
   len1  = (hpn+1)*(hpn+1);
   dlagP = (float *) malloc(len1*sizeof(float));
   dlagP_JeI  = (float *) malloc(len1*sizeof(float));
   dlagP_Gllw = (float *) malloc(len1*sizeof(float));
   lagP  = (float *) malloc(len1*sizeof(float));
   GLL(hpn, gllp);
   GLLWeight(gllp, hpn, gllw);
   printf("Lagrange function value, degree=%d\n",hpn);
/* For saving time, calculate Lagrange function and its derivative in 2-D*/
   for(i=0;i<=hpn;i++){
      printf("gllp[%d]=%f gllw[%d]=%f\n",i,gllp[i],i,gllw[i]);
      for(j=0;j<=hpn;j++){
//         dlagP[j*(hpn+1)+i] = dLag_N2(hpn,i,j,gllp);
           /* dLag_N_i() calculate l'i[j] stored as dlagP[j][i], lagP[j][i]  */
           dlagP[j*(hpn+1)+i] = dLag_N_i(hpn,i,gllp,gllp[j]);
           lagP[j*(hpn+1)+i]  = Lag_N_i(hpn,i,gllp,gllp[j]);
//           printf(" %8.5f | ", lagP[j*(hpn+1)+i]);
           printf(" %8.5f | ", dlagP[j*(hpn+1)+i]);
      }
      printf("\n");
   }

/*If using attenuation absorbing boundary*/
   if(atten){
       taper = (float *) malloc (nh*sizeof(float));
       atten_taper(nh,ass,taper);
//       for(i=0;i<nh;i++) printf("taper[%d]=%f \n",i,taper[i]);
   }
/* ----- Part 2: discretization ----- */
   int meshsize;
   struct cmatrix * meshs; /*store matrix of the stiffness tensor, including rho and Cijkl*/
   meshsize = discret(flat,&meshs);
   printf("eh=%f h=%f\n",eh,h);
   /*check the validation of meshing by output the model */
   FILE * fop;
   if(check == 1) {
        float gx,gz;
        fop = fopen("grid.xz","wt");
        for(i=0; i<nx; i=i+2){
            for(j=0;j<nz;j=j+2){
                gx= i* h; gz= (nz-1-j)*h;
                fprintf(fop,"%-8.2f %-9.2f %9.4f\n",gx,gz, meshs[j*nx+i].rho);
            }
            fprintf(fop,"%-8.2f %-9.2f %9.4f\n",gx,gz, meshs[(nz-1)*nx+i].rho);
        }
        fclose(fop);

        float lx,leg,temp;
        fop = fopen("legendre.xz","wt");
        for(i=0; i<=hpn; i++){
            fprintf(fop,"> %-9.5f\n",gllp[i]);
            for(lx=-1;lx <= 1.0;lx=lx+0.05){
                leg = Lag_N_i(hpn,i,gllp,lx);
                fprintf(fop,"%-8.2f %-9.3f\n",lx,leg);
            }
        }
        fclose(fop);
   }
/* ----- Part 3: Jacobian matrix -----
   construct Jacobian array for all the elements
   If the domain is partitioned uniformly, then Je= eh*eh/4
                               inverse[]=[2/eh,0,0,2/eh];
*/
      float je0, jmax=0.; 
      int are_same = 1;
      All_Je(&JeArray);
      je0 = JeArray[0];
      for(i = 1; i < nx*nz; i++){
          if(fabs(je0 - JeArray[i]) > 1.e-2) {
             are_same = 0;
             if(fabs(je0 - JeArray[i]) > jmax) jmax = fabs(je0 - JeArray[i]);
          }
      }

      if(are_same == 1) {
          Uni_Je  = (eh*eh/4.); // Physical meaning: Area ratio between two coordinate systems
          Uni_JeI = 2./eh;
          Uni_Jb  = eh/2.;
          Uni_JbI = 2./eh;
      }
      printf("jmax= %f Uniform=%d meshs Je=%f Je_Inv=%f\n",jmax, are_same,Uni_Je,Uni_JeI);

/* For saving time, storing Lagrange function and its derivative multiple Uni_Jei */
    for(i=0;i<=hpn;i++){
      for(j=0;j<=hpn;j++){
           dlagP_JeI[j*(hpn+1)+i] = dlagP[j*(hpn+1)+i] * Uni_JeI;
           dlagP_Gllw[j*(hpn+1)+i] = dlagP[j*(hpn+1)+i] * gllw[j];
      }
    }
/* ----- Part 4: Mass Matrix -----
   Construct Mass Matrix for all the elements
*/
   float * Mass, *Mass_dt2;
   len = (long) (nx * nz +1)* sizeof(float);
   Mass     = (float *) malloc(len);
   Mass_dt2 = (float *) malloc(len);
   Mass_Matrix(Mass, uniform, meshs);

   for(i=0;i<nx*nz;i++)  Mass_dt2[i] = sdt*sdt/Mass[i];

/* ----- Part 5: Absorbing boundary -----*/
   struct vector * Ab_matrix;
   if(absorb >0){
       printf("Using the absorbing boudary condition\n");
       Ab_matrix = (struct vector *) malloc((nx*nz+1)*sizeof(struct vector));
       for(i=0; i<nx*nz; i++)
          for(j=0;j<3;j++) Ab_matrix[i].u[j] = 0.;
       Absorb_Matrix(Ab_matrix,meshs,uniform,sdt);
       if(check ==1){
           for(i=0;i<nx;i=i+200)
              for(j=0;j<nz; j=j+200) printf("mass[%d][%d]=%f ab= %f\n", i,j, Mass[j*nx+i],Ab_matrix[j*nx+i].u[0]);
              printf("mass[283][0]= %9.4f Ab=%9.4f \n", Mass[283],Ab_matrix[283].u[0]);
       }
   }
   

/* Part 6.0 For output at the receivers*/
   if(plot_trace){
       /*  Output for test only */
       fopu1 = open_file("seis.u","wb");
       fopw1 = open_file("seis.w","wb");
       fopv1 = open_file("seis.v","wb"); /* zhl, 7/1/2005 */
       ixbeg = (int) (xbeg/h - 1);
       izbeg = (int) (zbeg/h);
       izbeg += nz -1;
       if(zbeg == 0.) izbeg = nz -1;
       if(izbeg > nz-1) izbeg= nz-1;
       printf("ixbeg=%d, izbeg=%d, idxplot=%d, idzplot=%d\n",ixbeg,izbeg,idxplot,idzplot);
       if(get_rec == 0){ /* if not assign receiver's position by raw data */
          for(itrace=0; itrace<xtrace; itrace++){
              for(j=0; j<ztrace; j++){
                  ix1[itrace+j*xtrace] = ixbeg + idxplot*itrace;
                  iz1[itrace+j*xtrace] = izbeg - idzplot*j;
                  if(itrace==0)printf("%d  %d  iz1 = %d\n",j,j*xtrace,iz1[itrace+j*xtrace]);
              }
          }
       }
       else{ /* zhl, 10/27/2006 */
          ztrace = 1;
          rfp  = open_file(rec_file,"rt");
          fop = fopen("receiver.xz","wt");
          fscanf(rfp,  "%d\n", &xtrace);
          fprintf(fop,"%d\n\n",xtrace);
          for(itrace=0; itrace<xtrace; itrace++){
              fscanf(rfp,"%f", &rec_x);
/*              ixrec = (int) ((rec_x*111.2 - xmin)/h); 
              fprintf(fop,"%-9.4f %9.4f\n",(rec_x*111.2 - xmin),0);
  20140108 commented by zhl*/
              ixrec = (int) (rec_x/h);
              fprintf(fop,"%-9.4f %9.4f\n",rec_x,0);
              ix1[itrace] = ixrec;
              iz1[itrace] = izbeg;
          }
          fclose(rfp);
          fclose(fop);
       }
   }

/* ----- Part 6: Time proceeding, second-order difference ----- */
   int shtotal; /*like mtotal, but for SH incidenting */
   if(msource == 0){
      infile=open_file(greenfile,"rb");
      shinfile=open_file(shgreenfile,"rb");
      mtotal = 2*(nx + nz);
      shtotal = nx + nz;
      grtuw  = (float *) malloc(mtotal*sizeof(float));
      grtuw0 = (float *) malloc(mtotal*sizeof(float));

      grtv  = (float *) malloc(shtotal*sizeof(float));
      grtv0 = (float *) malloc(shtotal*sizeof(float));

      sy = matrix(0,(mtotal-1),0,3);
      sy_sh = matrix(0,(shtotal-1),0,3);
      ns = (int *) malloc(4*sizeof(int));
      rdt = sdt/dp;
      printf("mtotal= %d shtotal=%d\n", mtotal, shtotal);
      nt = (int) ((float)((nt-3)*dp/sdt));

      for(j=0;j<4;j++){
          ns[j] = j;
          for(i=0;i<mtotal;i++) sy[i][j]= 0.0;
          for(i=0;i<shtotal;i++) sy_sh[i][j]= 0.0;
      }
   } else{
      Ricker = (float *) malloc((nst+1)*sizeof(float));
      Ricker_Wavelet(Ricker,nst,sdt,f0,t0,a0);
      fop = fopen("ricker.xz","wt");
         for(i=100;i<500; i++) fprintf(fop,"%f %e\n",i*sdt,Ricker[i]/a0);
      fclose(fop);
   }

   len = (long) (nx * nz + 1) * sizeof(struct vector);
   struct vector * uvw, *uvw1, *uvw0;/*displacement field at time point t, t+dt, and t-dt */
   struct vector * GStiff; /*stiffness matrix for all elements*/
   struct vector Force;    /*Force item at time t */
   for(i=0;i<3;i++) Force.u[i] = 0.;
   uvw0   = (struct vector *) malloc(len);
   uvw    = (struct vector *) malloc(len);
   uvw1   = (struct vector *) malloc(len);
   GStiff = (struct vector *) malloc(len);

   float grtmax,umax;
   int gl,ggg,uuu;

   for(it = 0;it < nst; it++){
      if(it>0) {

           /*20140329, using openMP works*/
//           #pragma omp parallel for 
           for(gl=0;gl<mtotal;gl++){ grtuw0[gl] = grtuw[gl];}
           for(gl=0;gl<shtotal;gl++){ grtv0[gl] = grtv[gl];}
      }

      if(!(it%20)) printf("--- Run %d ---\n",it);
      if(msource == 0) {
            if(it < nt){
                interpl(infile,  rec_len,it,rdt, mtotal,ns,sy,   grtuw,nx,nz);
                interpl(shinfile,rec_len,it,rdt,shtotal,ns,sy_sh,grtv, nx,nz);        
// To avoid singularity in GRT resolution
                for(gl=0;gl<mtotal;gl++){
                    if(!(grtuw[gl]<10 && grtuw[gl]>-10)) {
                        // printf("gl=%d grtuw=%e\n",gl,grtuw[gl]);
                        grtuw[gl] = (grtuw[gl-1]+grtuw[gl+1])/2;
                    }
                }
                for(gl=0;gl<shtotal;gl++){
                    if(!(grtv[gl]<10 && grtv[gl]>-10)) {
                        grtv[gl] = (grtv[gl-1]+grtv[gl+1])/2;
                    }
                }

            }
      }

      if(it==0) { 
          Initial_uvw(uvw);
          Initial_uvw(uvw0);
          if(msource == 0)  GRT_Interface(uvw,grtuw,grtv);
      }
      else{
          if(msource == 1) Explosive(mex,mez,Ricker[it],&Force);
//          if(msource == 0) { if(it < nt) GRT_Interface2(uvw,grtuw,grtuw0);}
/*zhl: 20131212 change to using following subroutine*/
//          if(msource == 0) { if(it < nt) GRT_Interface3(uvw,uvw0,grtuw,grtuw0,grtv,grtv0);}
          if(msource == 0) { if(it < nt) GRT_Interface4(uvw,uvw0,grtuw,grtuw0,grtv,grtv0);}  

          Stiff_Matrix(GStiff,uvw,meshs,uniform);

          if(msource == 1) Time_diff(it,uvw1,uvw,uvw0,GStiff,Mass,Mass_dt2, meshs,sdt,msource,mex,mez,Force,absorb,Ab_matrix); 
          else   Time_diff_hybrid(it,uvw1,uvw,uvw0,GStiff,Mass, Mass_dt2, meshs,sdt,absorb,Ab_matrix);
/*zhl: 20131210 remove the comment of the following line*/
/*zhl: 20140220 modify attenuation actions to avoid multiple-time attenuation at several corners*/
//          if(atten>0) atten_wave(uvw1,nh,atten,taper);
          if(atten>0) atten_wave2(uvw1,nh,atten,taper);
          Copy_Displacement(uvw, uvw0);
          Copy_Displacement(uvw1,uvw);
      }

/*
if((!(it%25) && it < nt) || it<20){
   grtmax = 0.; 
   umax=0.;
   for(gl=0;gl<2*(nx+nz);gl++) if(fabs(grtuw[gl])>grtmax) {ggg=gl; grtmax=fabs(grtuw[gl]);}
   for(gl=0;gl<nx;gl++) if(fabs(uvw[gl].u[0])>umax) {uuu=gl; umax=fabs(uvw[gl].u[0]);}
   printf("%d num=%d grtmax=%e unum=%d umax=%e grt_v=%e\n",it,ggg,grtmax,uuu,umax,grtuw[2*nx+30]);
}
*/

      if(plot_trace) kirrecord_zhl(uvw,iz1,ix1,xtrace,ztrace,fopu1,fopw1,fopv1,0);
      if(plot_snap){
	  /* Output snapshots */
          if(it==ist)snapshot("out1",uvw,nx,nz,2,2);
          if(it==ist+incre)snapshot("out2",uvw,nx,nz,2,2);
          if(it==ist+2*incre)snapshot("out3",uvw,nx,nz,2,2);
          if(it==ist+3*incre)snapshot("out4",uvw,nx,nz,2,2);
          if(it==ist+4*incre)snapshot("out5",uvw,nx,nz,2,2);
          if(it==ist+5*incre)snapshot("out6",uvw,nx,nz,2,2);
          if(it==ist+6*incre)snapshot("out7",uvw,nx,nz,2,2);
          if(it==ist+7*incre)snapshot("out8",uvw,nx,nz,2,2);
          if(it==ist+8*incre)snapshot("out9",uvw,nx,nz,2,2);
          if(it==ist+9*incre)snapshot("out10",uvw,nx,nz,2,2);
          if(it==ist+10*incre)snapshot("out11",uvw,nx,nz,2,2);
          if(it==ist+11*incre)snapshot("out12",uvw,nx,nz,2,2);
          if(it==ist+12*incre)snapshot("out13",uvw,nx,nz,2,2);
          if(it==ist+13*incre)snapshot("out14",uvw,nx,nz,2,2);
          if(it==ist+14*incre)snapshot("out15",uvw,nx,nz,2,2);
          if(it==ist+15*incre)snapshot("out16",uvw,nx,nz,2,2);
          if(it==ist+16*incre)snapshot("out17",uvw,nx,nz,2,2);
          if(it==ist+17*incre)snapshot("out18",uvw,nx,nz,2,2);
          if(it==ist+18*incre)snapshot("out19",uvw,nx,nz,2,2);
          if(it==ist+19*incre)snapshot("out20",uvw,nx,nz,2,2);
          if(it==ist+20*incre)snapshot("out21",uvw,nx,nz,2,2);
          if(it==ist+21*incre)snapshot("out22",uvw,nx,nz,2,2);
          if(it==ist+22*incre)snapshot("out23",uvw,nx,nz,2,2);
          if(it==ist+23*incre)snapshot("out24",uvw,nx,nz,2,2);
          if(it==ist+24*incre)snapshot("out25",uvw,nx,nz,2,2);
          if(it==ist+25*incre)snapshot("out26",uvw,nx,nz,2,2);
      }
   }
   printf("SEM calculation finished!\n");
   if(plot_trace) { fclose(fopw1); fclose(fopu1); fclose(fopv1);} 
   if(msource == 0) {
      fclose(infile);
      fclose(shinfile);
   }
}

/* ===========================================================================
   ===========================================================================
   ===========================================================================
*/
/* subroutines */
/*
   subroutine mesh2d() discretize the SEM region
   get the node parameters
 ---- global variants
   nex, nez,

*/

int discret(flat, mesh)
int flat;
struct cmatrix ** mesh; /*store matrix of the stiffness tensor*/ 
{
   int  tablesize;
   int  i,j,k;
   FILE *fop,*open_file();
   
   struct cmatrix * c_matrix; /*store matrix of the stiffness tensor*/ 
   struct cmatrix * ppc;
   void   read_model();
   float  Cijkl(), Cijkl_3D();
   float  thk2;

/* about read SEM model: np c11 c13 c16 c33 c36 c44 c45 c55 c66 rho
   vp ==> c11
   vs ==> c44
   add c13 c16 c33 c36 c45 c55 c66, 7 parameters
*/
   int   *np, nll; /*nll: layer number of the SEM model*/
   float *vp, *vs, *rho, **an, **bn; 
   float *mc12,*mc13,*mc16,*mc22,*mc23,*mc26,*mc33,*mc36,*mc45,*mc55,*mc66;
   float lip(); 
   float f0,x;
   int   *ratio;

   float fvp,fvs,fden; /*vp vs and rho at free surface */
   float ivp,ivs,iden; /*vp vs and rho at the incident boundary */

/* step 1: first, assign with a 1-D model such as the PREM 
           reading the medium parameters from the raymodel used in the GRT 
           calculation;  
   jo : layer number of the PREM     nb : source layer
   nen: layer number of the ith ray  nna: the layer number of the ith ray
   nray: property of ray             nrec: the last layer of GRT
*/
   int   jo, nb; 
   int   lfinal, *nen, *ncoun, *nna, *nray;
   float *cc,*ss,*dd,*tth; /*vp vs density and thickness of 1-D model*/
   int   firstlayer, nrec; /*nrec: the last layer of GRT calculation */
   float thickness;
/* temp variants using in loop of assign parameters*/
   int   nlay1,nlay2,nlaynum, inlay;
   int   iz1, iz2;
   float t_c11, t_c13, t_c44;
//   struct cmatrix2d apoint;
   struct cmatrix3d apoint3d;
   struct cmatrix3d apoint;
   struct cmatrix temp_c;

   firstlayer = 0;
   read_model(flat,raymodel,&jo,&nb,&lfinal,&cc,&ss,&dd,&tth,
		       &nen,&ncoun,&nna,&nray);
   nrec = nna[nen[0]-2];
   printf("nrec = %d,  nfil = %d\n",nrec,nfil);

   thickness=0.0;
   if(nfil >nrec){
       for(i=nrec; i<nfil;i++)
           thickness += tth[i];
   } else {
       for(i=nfil-1; i<nrec-1;i++) {
           thickness += tth[i];
	   printf("i = %d, thick= %e, depth= %f vp = %f\n",i+1,tth[i], thickness,cc[i]);
       }
   }
   nez = (int) ((thickness-1.e-3)/eh) ;
   printf("Warning: eh was changed from %f --->",eh);
   eh = (thickness-1.e-3)/nez;
   printf(" %f h= %f\n",eh,eh/hpn);
   nel = nex * nez;  /* element number */
   printf("thickness = %f, nez= %d\n", thickness, nez);
   
   ivp = cc[nrec-1]; ivs=ss[nrec-1]; iden= dd[nrec-1]; /*the last layer of GRT calculation*/

   long int len;
   nx = nex * hpn + 1;  /* horizontal grid points*/
   nz = nez * hpn + 1;  /* vertical grid points  */
   h = eh/hpn;    /* grid point spacing */
   len = (long) (nx * nz) * sizeof(struct cmatrix);
   *mesh = (struct cmatrix *) malloc(len);
   c_matrix= *mesh;
   ppc= c_matrix;
   if(*mesh == NULL){
      printf("cannot allocate meshs memory in subroutine discret()\n");
      exit(-1);
   }

   iz1 = 0;
/*20110602: test changing nrec-2 to nrec -1*/
/*20140117: test changing nrec-1 to nrec -2*/
   nlay1 = ((nfil>nrec)? nrec:nrec-2);
   nlay2 = ((nfil>nrec)? nfil+200:nfil-200);
   if(nlay2 < 0) nlay2=0;
   nlaynum = nlay2 - nlay1;
   if(nlaynum <0) nlaynum *= -1;
   nlaynum += 1;
   inlay = ((nfil>nrec)? 1:-1);
   printf("nlay1=%d nlay2=%d nlaynum=%d\n",nlay1, nlay2,nlaynum);

/*20140116: thickness -> thk2 */
   thk2 = 0.0;
/*From down to up */
/*  
 |  nz-1 --------------------
 |  nz-2 --------------------  
 |  ...
 |  0    --------------------
*/
   int cnt,i1,j1,k1,l1;
   for(k=0; k < nlaynum; k++){
      i = nlay1 + k*inlay;
      thk2 += tth[i];
      iz2 = (int) ((thk2 + 1.e-3)/h) + 1;
      if(iz2 > nz-1) iz2 = nz-1;
      if(iz2 < iz1 ) break;
/* vp = sqrt(c11/rho) vs= sqrt(c44/rho): getting c11 and c44 */
      ivp = cc[i]* cc[i]; ivs= ss[i]*ss[i]; iden= dd[i];
      t_c11 = ivp*iden; t_c44 = ivs*iden;
      t_c13 = t_c11 - 2.0* t_c44;
//      printf("iz1=%4d iz2=%4d vp=%f vs=%f den=%f c11=%f c13=%f c44=%f\n",iz1,iz2,cc[i],ss[i],dd[i],t_c11,t_c13,t_c44);
      apoint3d.c11= t_c11; apoint3d.c12= t_c13;  apoint3d.c13= t_c13;  apoint3d.c16= 0.;
      apoint3d.c22= t_c11; apoint3d.c23= t_c13;  apoint3d.c26= 0.0;
      apoint3d.c33= t_c11; apoint3d.c36= 0.0;    apoint3d.c44= t_c44;
      apoint3d.c45= 0.0 ;  apoint3d.c55= t_c44;  apoint3d.c66= t_c44;
      for(i1=0;i1<3;i1++){
           for(j1=0;j1<3;j1++){
              for(k1=0;k1<3;k1++){
                 for(l1=0;l1<3;l1++) temp_c.c[i1][j1][k1][l1] = Cijkl_3D(i1,j1,k1,l1,apoint3d);
              }
           }
      }

      for(j=iz1;j<=iz2;j++){
          cnt = nx;
          for(ppc= c_matrix + j*nx; cnt--; ppc++){
               ppc[0].rho = iden;
               for(i1=0;i1<3;i1++){
                   for(j1=0;j1<3;j1++){
                      for(k1=0;k1<3;k1++){
                         for(l1=0;l1<3;l1++) ppc[0].c[i1][j1][k1][l1] = temp_c.c[i1][j1][k1][l1];
                      }
                   }
               }
          }
      }
      iz1 = iz2 + 1;
      if(iz2 == nz -1) break;
   }
   free(nen); free(nna);free(ncoun); free(nray);

   /*assign with the input model*/
   ppc= c_matrix;
   if(readpars == 1){
        fop = open_file(semmodel,"r");
        fscanf(fop,"%d\n",&nll);
        printf("nll= %d\n",nll);
        if(nll<0){
           firstlayer= 1; nll *=-1;
        }
        ratio   = (int   *) malloc ((nll+2)*sizeof(int)); 
        np   = (int   *) malloc ((nll+2)*sizeof(int)); 
        vp   = (float *) malloc ((nll+2)*sizeof(float)); 
        vs   = (float *) malloc ((nll+2)*sizeof(float)); 
        /* vp==>c11 vs==>c44 and mc13 mc16 mc33 mc36 mc45 mc55 mc66 in order */
        mc12 = (float *) malloc ((nll+2)*sizeof(float));
        mc13 = (float *) malloc ((nll+2)*sizeof(float));  
        mc16 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc22 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc23 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc26 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc33 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc36 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc45 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc55 = (float *) malloc ((nll+2)*sizeof(float)); 
        mc66 = (float *) malloc ((nll+2)*sizeof(float)); 
        rho  = (float *) malloc ((nll+2)*sizeof(float)); 
        /*05252011, changed to 1000 from 6000*/
        an   = matrix(0,nll-1,0,1000);
        bn   = matrix(0,nll-1,0,1000);

        for(j=0; j<nll; j++){
            ratio[j] = 0;
            fscanf(fop, "%d %f %f %f %f %f %f %f",   &np[j],&vp[j],   &mc12[j],&mc13[j],&mc16[j],&mc22[j],&mc23[j],&mc26[j]);
            fscanf(fop, "%f %f %f %f %f %f %f\n",           &mc33[j], &mc36[j],&vs[j],  &mc45[j],&mc55[j],&mc66[j],&rho[j]);
            printf("%-4d   %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f",np[j],vp[j],  mc12[j],mc13[j],mc16[j],mc22[j],mc23[j],mc26[j]);
            printf("%-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f\n",           mc33[j],mc36[j],vs[j],  mc45[j],mc55[j],mc66[j],rho[j]);
            if(np[j] < 0){
                ratio[j] = 1;
                np[j]    = - np[j];
            }
            if(np[j] > 1000){
                printf("an bn allocation error, np[j] >6000\n");
                exit(-2);
            }
            for(i=0; i<np[j]; i++){
                fscanf(fop,"%f ",&f0);
                an[j][i] = f0;
            } 
            fscanf(fop,"\n");
            for(i=0; i<np[j]; i++){
                fscanf(fop,"%f ",&f0);
                bn[j][i] = f0;
            } 
            fscanf(fop,"\n");
        }
        fclose(fop);

        printf("assign parameters for the SEM region: 0->nx=%d, 0->nz=%d\n", nx, nz);
        int iz;
        /*from down to up*/
        for(i=0; i< nx; i++){  
            x= i* h;
            iz1 = 0;
            if(firstlayer ==1) 
                iz1 = (int) ((lip(0,np,an,bn,x)+1.e-3)/h);
            if(iz1 < 0) iz1 = 0;
            for(j=firstlayer; j<nll; j++){
                iz2 = (int) ((lip(j,np,an,bn,x)+1.e-3)/h) ;
	        if(iz2 > nz-1) iz2 = nz-1;
//		if(!(i % 100)) printf("x=%f, j = %d, iz1 = %d,  iz2 = %d\n",x,j,iz1,iz2);
                t_c11 = vp[j]; t_c44 = vs[j]; /* zhl, 2/27/2006 */
                t_c13 = mc13[j];              /* zhl, 2/27/2006 */
                if(ratio[j] == 1){ /*if assign by multiplying with coefficients*/
                    apoint.c11 *= t_c11;   apoint.c12 *= mc12[j]; apoint.c13 *= mc13[j]; apoint.c16 *= mc16[j];
                    apoint.c22 *= mc22[j]; apoint.c23 *= mc23[j]; apoint.c26 *= mc26[j]; apoint.c33 *= mc33[j]; apoint.c36 *= mc36[j];  
                    apoint.c44 *= t_c44;   apoint.c45 *= mc45[j]; apoint.c55 *= mc55[j]; apoint.c66 *= mc66[j];
                    apoint.rho *= rho[j];
                } else {  /* real value */
                    apoint.c11 = t_c11;   apoint.c12 = mc12[j]; apoint.c13 = mc13[j]; apoint.c16 = mc16[j];
                    apoint.c22 = mc22[j]; apoint.c23 = mc23[j]; apoint.c26 = mc26[j]; apoint.c33 = mc33[j]; apoint.c36 = mc36[j];  
                    apoint.c44 = t_c44;   apoint.c45 = mc45[j]; apoint.c55 = mc55[j]; apoint.c66 = mc66[j];
                    apoint.rho = rho[j];
                }
                for(i1=0;i1<3;i1++){
                    for(j1=0;j1<3;j1++){
                       for(k1=0;k1<3;k1++){
                          for(l1=0;l1<3;l1++) {
                              temp_c.c[i1][j1][k1][l1] = Cijkl_3D(i1,j1,k1,l1,apoint);
                          }
                       }
                    }
                }
                temp_c.rho= apoint.rho;
                /*For a layer */
                for(iz = iz1; iz<=iz2; iz++){
                    for(i1=0;i1<3;i1++){
                        for(j1=0;j1<3;j1++){
                           for(k1=0;k1<3;k1++){
                              for(l1=0;l1<3;l1++) ppc[iz*nx+i].c[i1][j1][k1][l1] = temp_c.c[i1][j1][k1][l1];
                           }
                        }
                    }
                    ppc[iz*nx+i].rho= temp_c.rho;
                }
               iz1 = iz2 + 1;
            }
        }

        free(np); free(vp); free(vs); free(rho);
        free(mc13);free(mc16);free(mc33);free(mc36);free(mc45);free(mc55);free(mc66);
        free(ratio);
        /*05252011,ZHL, changed from 6000 to 1000 */
        free_dmatrix(an,0,nll-1,0,1000);
        free_dmatrix(bn,0,nll-1,0,1000);
   }
   free(cc); free(ss); free(dd); free(tth);
   printf("thk2=%-8.2f\n",thk2);
   mesh= &c_matrix;

/*20140115 output model to file for check*/
   if(check == 1){
      FILE * fom;
      float tempv,tv44,tv55;
      fom = open_file("vel.xy","w");
      for(i=0;i<nx;i++){
          for(j=0;j<nz;j++){
              if(output_vd ==0){
                 tempv = sqrt(ppc[j*nx+i].c[0][0][0][0]/ppc[j*nx+i].rho);
                 fprintf(fom,"%-8.3f %-8.2f %-8.2f\n", (xmin+i*h)/(6371*3.14159265/180.),(nz-1-j)*h,tempv);
              }
              else{
                 tv44 = sqrt(ppc[j*nx+i].c[1][2][1][2]/ppc[j*nx+i].rho);
                 tv55 = sqrt(ppc[j*nx+i].c[0][2][0][2]/ppc[j*nx+i].rho);
                 tempv = (tv44-tv55)/(tv44+tv55)/2.;
                 fprintf(fom,"%-8.3f %-8.2f %-8.2f\n", (xmin+i*h)/(6371*3.14159265/180.),(nz-1-j)*h,tempv);
              }
          }
      }
      fclose(fom);
   }

   tablesize = nx*nz;
   return(tablesize);
}

/*subroutine for interpolating the SEM model by Control Depth Points*/
float lip(layer,np,an,bn,x)
int layer, *np; float **an, **bn; float x;
{
    float dep; int i;
    if( x < an[layer][0] || x > an[layer][np[layer]-1] ){
        printf("Model out of range, x=%9.3f an= %9.3f or %9.3f \n", x, an[layer][0], an[layer][np[layer]-1]);
        printf("layer= %d \n", layer);
        exit(-3);
    }
    for(i=0; i<np[layer]-1; i++)
        if(x >=an[layer][i]  && x <an[layer][i+1]){
            dep  = (bn[layer][i+1]-bn[layer][i])/(an[layer][i+1]-an[layer][i]);
            dep  = bn[layer][i] + dep*(x-an[layer][i]); 
        }
    return dep;
}

/* All_Je() construct Jacobian array for each element
   Calculating one element by element
   But assemble to Grid-Array for storing
Requires
   global: nex, nez : elemental number, horizontal and vertical 
   eh: element spacing
   gllp: GLL integration points 
   Output: Je4_one array
*/
void All_Je(JeArray)
float ** JeArray;
{
   int i, j;/*element loop: x, z in order*/
   int m, n;/*grid loop in each element: x,z in order*/
   long int len;
   int index;
   float * Je4, rx,rz; /*rx,rz: GLL points */
   struct Point2d points[4];
   float  invm[4];
   float Jacobian_4();

   len = (long) (nx * nz+1) * sizeof(float) ;
   *JeArray = (float *) malloc(len);
   Je4 = *JeArray;
   printf("len=%d\n",nx*nz);
   for(i=0;i<nex;i++){
      for(j=0;j<nez;j++){
         points[0].x = i*eh;     points[0].z = j*eh;
         points[1].x = (i+1)*eh; points[1].z = j*eh;
         points[2].x = i*eh;     points[2].z = (j+1)*eh;
         points[3].x = (i+1)*eh; points[3].z = (j+1)*eh;
         for(m=0;m<=hpn;m++){
            for(n=0;n<=hpn;n++){
                rx = gllp[m]; rz= gllp[n];
                index = (j*hpn+ n) * nx + i*hpn + m;
                Je4[index] = Jacobian_4(rx,rz,points,invm);
//                printf("%d %d %f %f Je[%d]=%f Inv= %f\n", m,n,rx,rz,index, Je4[index],invm[0]);
            }
         }
      }
   }
}

/*
Quadrilateral element with 4 control points
     3(x3,z3)------4(x4,z4)       (r3,t3)-----(r4,t4)
     |             |        --->  |-1,1         |1,1
     |             |              |             | 
     |             |              |-1,-1        |1,-1
     1(x1,z1)------2(x2,z2)       (r1,t1)-----(r2,t2)
     ex: column number of the element
     ez: raw number of the element
Requires:
     Global: nx, nz,hpn,eh
*/

/*
Get idex in mesh arrays: from down to up 0->nez-1, from left to right 0->nex-1
    Input: 
      ex: column number of the element
      ez: raw number of the element
    Output:
      return index of the grid point at the left corner of the element
Requires:
    Global: nx,nz,hpn
*/
int GetIndex(ex,ez)
int ex,ez;
{
    int index;
    index= ez*hpn*nx + ex*hpn ; /*if you check, please be careful*/
    return index;
}

/*
  Calculate system equation by second-difference of Time
Require:
  uvw1,uvw,uvw0: displacement vector at time t+dt, t, and t-dt
  GStiff: global stiffness matrix
  Mass: global mass matrix
  sdt: time step
  Global: nx,nz
  msource: using Moment Source at the (mex,mez) element
  Force: Force item by Moment source
  absorb: 1: using absorbing boundary condition
         note that ab_matrix has multiplied sdt/2
*/
void Time_diff(it,p1,p,p0,GStiff,Mass,Mass_dt2, meshs,sdt,msource,mex,mez,Force,absorb,ab_matrix)
int    it;
struct vector * p1,*p,*p0;
struct vector * GStiff;
float *Mass, *Mass_dt2, sdt;
struct cmatrix * meshs;
int msource,mex,mez;
struct vector Force;
int absorb;
struct vector * ab_matrix;
{
   int i,j,k,gn,index;
   float tempf,coe,coe1,dt2;
   void Absorbs(), Absorbs2();

   dt2 = sdt*sdt;
   if(msource == 1){
      index = (mez*hpn)*nx + mex* hpn ; /*The left-bottom corner of (mex,mez) element */
      if(absorb==1){
          for(i=0;i<index;i++){
//             coe = -dt2/Mass[i];
             for(k=0;k<3;k++){
                p1[i].u[k] = -Mass_dt2[i]*GStiff[i].u[k] + 2.*p[i].u[k] - p0[i].u[k];
             }
          }
          for(k=0;k<3;k++){
                tempf = Force.u[k] - GStiff[index].u[k];
                p1[index].u[k] = tempf*Mass_dt2[index] + 2.*p[index].u[k] - p0[index].u[k];
          }
          for(i=index+1;i<nx*nz;i++){
//             coe = -dt2/Mass[i];
             for(k=0;k<3;k++){
                p1[i].u[k] = -Mass_dt2[i]*GStiff[i].u[k] + 2.*p[i].u[k] - p0[i].u[k];
             }
          }
          Absorbs(p1,p,p0,GStiff,Mass,sdt,ab_matrix);/*Jan 10, 2008*/
//          Absorbs2(p1,meshs,sdt);
      } else{
          for(i=0;i<index;i++){
//             coe = -dt2/Mass[i];
             for(k=0;k<3;k++){
                p1[i].u[k] = -Mass_dt2[i]*GStiff[i].u[k] + 2.*p[i].u[k] - p0[i].u[k];
             }
          }
          for(k=0;k<3;k++){
                tempf = Force.u[k] - GStiff[index].u[k];
                p1[index].u[k] = tempf*Mass_dt2[index] + 2.*p[index].u[k] - p0[index].u[k];
          }
          for(i=index+1;i<nx*nz;i++){
//             coe = -dt2/Mass[i];
             for(k=0;k<3;k++){
                p1[i].u[k] = -Mass_dt2[i]*GStiff[i].u[k] + 2.*p[i].u[k] - p0[i].u[k];
             }
          }
      }
   } 
}

/*
  Calculate system equation by second-difference of Time without internal source
Require:
  uvw1,uvw,uvw0: displacement vector at time t+dt, t, and t-dt
  GStiff: global stiffness matrix
  Mass: global mass matrix
  sdt: time step
  Global: nx,nz
  absorb: 1: using absorbing boundary condition
         note that ab_matrix has multiplied sdt/2
*/
void Time_diff_hybrid(it,p1,p,p0,GStiff,Mass,Mass_dt2,meshs,sdt,absorb,ab_matrix)
int    it;
struct vector * p1,*p,*p0;
struct vector * GStiff;
float *Mass, *Mass_dt2, sdt;
struct cmatrix * meshs;
int absorb;
struct vector * ab_matrix;
{
   int i,k;
   float coe,dt2;
   void Absorbs(), Absorbs2();

   dt2 = sdt*sdt;

   for(k=0;k<3;k++){
    /*20140327, the following line of openMP works */
//   # pragma omp parallel for private (coe)
      for(i=0;i<nx*nz;i++){
//         coe = -dt2/Mass[i];
         p1[i].u[k] = -Mass_dt2[i]*GStiff[i].u[k] + 2.*p[i].u[k] - p0[i].u[k];
      }
   }
   Absorbs(p1,p,p0,GStiff,Mass,sdt,ab_matrix);
//      Absorbs2(p1,meshs,sdt);
}

/*
   Ref: Komatitsch and Tromp, 1999, GJI, Page 812
   Mass_Matrix() calculated mass matrix()
   (1) calculate each element one by one
   (2) Assemble to Grid Array
Requires
   density: density at each grid, stored in mesh
   positions: of the anchors of each element 
   gllw: GLL points' weights for integration
   Je: Jacobian matrix
*/
void Mass_Matrix(mass, uniform, meshs)
float * mass;
int uniform;/*Spatial partition uniformly or not, if uniformly, Jacobian is very simple*/
struct cmatrix * meshs;
{
   int i, j;/*element loop: x, z in order*/
   int m, n;/*grid loop in each element: x,z in order*/
   int index;
   float rx,rz, density,je; /*rx,rz: GLL points' weight for integration*/

   for(i=0;i<nx*nz;i++)  mass[i] = 0.;

   for(i=0;i<nex;i++){
      for(j=0;j<nez;j++){
         for(m=0;m<=hpn;m++){
            rx = gllw[m]; 
            /*20140327, using openMP, works */
//            #pragma omp parallel for private (rz,index,density,je)
            for(n=0;n<=hpn;n++){
                rz= gllw[n];
                index = (j*hpn+ n) * nx + i*hpn + m;
                density = meshs[index].rho;
                if(uniform ==1) je= Uni_Je;
                else je = JeArray[index];
                mass[index] = mass[index] + rx * rz * je * density;
            }
         }
      }
   }
}

/*
  Calculate absorbing boundary condition 
  for saving time in time difference here BE=BE*sdt/.
Requires: 
  Ab_matrix: absorbing coefficient for all the grids
             
  meshs: storing mesh properties
  uniform: partition uniformly or not
  Globals:
      nex, nez, hpn

  nz-1 |                    | 
       |                    |
       |                    |
       |                    |
       |                    |
       |                    | 
    0   ---------------------
         0                nx-1
*/
void Absorb_Matrix(Ab_matrix,meshs,uniform,sdt)
struct vector * Ab_matrix;
struct cmatrix * meshs;
int uniform;
float sdt;
{
  int ex, ez;
  int i,j, al, be;
  int index;
  struct cmatrix  mesh;
  float vp,vs;

  /*left boundary x=0; n-> -x ; t1-> z; t2 -> y;*/
  for(ez=0;ez<nez;ez++){
     for(be=0;be<=hpn;be++){
        index = (ez*hpn + be) * nx ;
        mesh = meshs[index];
        vp = sqrt( mesh.c[0][0][0][0]/mesh.rho);
        vs = sqrt( mesh.c[1][2][1][2]/mesh.rho);
  /*20110622: change to ab += gll ... */
        Ab_matrix[index].u[0] +=   gllw[be] * mesh.rho * vp *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[1] +=   gllw[be] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[2] +=   gllw[be] * mesh.rho * vs *Uni_Jb*sdt/2.;
//if(!(ez%hpn) && be==hpn) printf("ez=%d ab=%f vp=%f vs=%f \n",ez,Ab_matrix[index].u[1],vp,vs);
     }
  }
  /*bottom boundary z=0; n-> -z; t1-> x; t2 -> y;*/
  for(ex=0;ex<nex;ex++){
     for(al=0;al<=hpn;al++){
        index = ex*hpn + al;
        mesh = meshs[index];
        vp = sqrt( mesh.c[0][0][0][0]/mesh.rho);
        vs = sqrt( mesh.c[1][2][1][2]/mesh.rho);
        Ab_matrix[index].u[0] +=   gllw[al] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[1] +=   gllw[al] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[2] +=   gllw[al] * mesh.rho * vp *Uni_Jb*sdt/2.;
//printf("ex=%d ab=%f vp=%f vs=%f \n",ex,Ab_matrix[index].u[1],vp,vs);
     }
  }  
  /*right boundary x=nx-1; n-> x; t1-> -z; t2 -> y;*/
  for(ez=0;ez<nez;ez++){
     for(be=0;be<=hpn;be++){
        index = (ez*hpn + be) * nx + nx -1 ;
        mesh = meshs[index];
        vp = sqrt( mesh.c[0][0][0][0]/mesh.rho);
        vs = sqrt( mesh.c[1][2][1][2]/mesh.rho);
        Ab_matrix[index].u[0] +=   gllw[be] * mesh.rho * vp *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[1] +=   gllw[be] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[2] +=   gllw[be] * mesh.rho * vs *Uni_Jb*sdt/2.;
//if(!(ez%hpn) && be==hpn) printf("ez=%d ab=%f vp=%f vs=%f \n",ez,Ab_matrix[index].u[1],vp,vs);
     }
  }
  /*top boundary z=0; n-> z; t1-> x; t2 -> y;*/
/* July 06, 2011, added top boundary absorbing, the waveform is clearer, and the corner effects are reduced
   but the reflection is weaker!
   but Yann: not necessary, because it is the free surface boundary
  
  for(ex=0;ex<nex;ex++){
     for(al=0;al<=hpn;al++){
        index = (nz-1) * nx + ex*hpn + al;
        mesh = meshs[index];
        vp = sqrt( mesh.c[0][0][0][0]/mesh.rho);
        vs = sqrt( mesh.c[1][2][1][2]/mesh.rho);
        Ab_matrix[index].u[0] += gllw[al] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[1] += gllw[al] * mesh.rho * vs *Uni_Jb*sdt/2.;
        Ab_matrix[index].u[2] += gllw[al] * mesh.rho * vp *Uni_Jb*sdt/2.;
     }
  }  
*/
}

/*
  Calculate Stiffness matrix for all the elements
Input:
  uvw: displacement array
  meshs: storing mesh properties
  uniform: partition uniformly or not
  Global: nex, nez
          hpn
*/
void Stiff_Matrix(Stiffs,p,meshs,uniform)
struct vector * Stiffs;
struct vector * p;
struct cmatrix * meshs;
int uniform;
{
   int ex,ez,i,j,k,l,index;
   int alpha, beta;/*at each element, the reference*/
   void Grid_Stiff();
   float dUj_di(); 
   float Stress_ij();
   float F_ik();
   void dUj_di_2D();

   float  a_Stiff[3];
   struct cmatrix  mesh;
   struct FIK * FikM; /*Matrix storing F_ik for all grid points in an element*/
   float  ukl[3][3];  /*ukl: displacement vector derivative: 0,1,2 component */
   float  Tij[3][3];  /*T[i][j]: stress  vector: 0, 1, 2 components */

   FikM = (struct FIK *) malloc((hpn+1)*(hpn+1)*sizeof(struct FIK));
/* initialize first */
   for(j=0;j<3;j++) {
       /*20140329, using openMP */
//       #pragma omp parallel for 
       for(i=0;i<nx*nz;i++){
            Stiffs[i].u[j] = 0.;
       }
    }

/* calculate all the element one by one, some share points between two elements need to add*/
   for(ex=0;ex<nex;ex++){
       for(ez=0;ez<nez;ez++){
         /*calculate Fik[alpha][beta] in each element */

         for(alpha=0;alpha<=hpn;alpha++){
            for(beta=0;beta<=hpn;beta++){
               index = (ez*hpn+beta)*nx + ex*hpn + alpha;
               mesh = meshs[index];
               /* using the symmetrical property Uk,l = Ul,k*/
/*
               for(k=0;k<3;k++){
                  ukl[k][0]= dUj_di(ex,ez,alpha,beta,p,k,0,uniform);
                  ukl[k][1]= 0.;
                  ukl[k][2]= dUj_di(ex,ez,alpha,beta,p,k,2,uniform);
               }
*/
               dUj_di_2D(ex,ez,alpha,beta,p,ukl,uniform);
               for(i=0;i<3;i++) {
/*071209, why Tij[i][1] = 0: see FikM equation in which Ti1 is not used since d/dx1=0 (Eq.A5, Komatitsch, GJI 1999) */
/*071509, IMPORTANT Questions 
   Tij[i][1] = 0, means T01 T11 and T12 (T12 T22 and T32 ) = 0
   But in My JGR 2008 Paper, only T22 = 0
   Try T32 = T23, T12 = T21
*/
                  Tij[i][0] = Stress_ij(i,0,ukl,mesh,uniform);
//                  Tij[i][1] = Stress_ij(i,1,ukl,mesh,uniform);
                  Tij[i][1] = 0.;
                  Tij[i][2] = Stress_ij(i,2,ukl,mesh,uniform);
               }
/*071509 add 2 lines, anyway, it is right although they are not used in Fik calculation since inverse Jacabian matrix=0 for axis y*/
//               Tij[0][1]= Tij[1][0];
//               Tij[2][1]= Tij[1][2];

//if(ex==0 &&ez==16&&alpha==0&&beta==0) printf("Tij=%e %e %e %e %e %e\n", Tij[0][0],Tij[0][2],Tij[1][0],Tij[1][2],Tij[2][0],Tij[2][2]);
               for(i=0;i<3;i++) {
                   FikM[beta*(hpn+1)+alpha].F0[i]= F_ik(i,0,Tij,uniform);
//                   FikM[beta*(hpn+1)+alpha].F1[i]= F_ik(i,1,Tij,uniform);
                   FikM[beta*(hpn+1)+alpha].F2[i]= F_ik(i,2,Tij,uniform);
               }
            }
         }

         for(alpha=0;alpha<=hpn;alpha++){
            for(beta=0;beta<=hpn;beta++){
               index = (ez*hpn+beta)*nx + ex*hpn + alpha;
               Grid_Stiff(alpha,beta,FikM,a_Stiff,uniform);
               for(j=0;j<3;j++) Stiffs[index].u[j] += a_Stiff[j];
            }
         }
       }
   }
   free(FikM);
}

/*
  Calculate Stiffness Matrix at a Grid point at an elemental level
Input
  alpha,beta: reference number horizontal and vertical
  uvw: displacement 
  stiff: storing mesh properties
Requires:
  Global:
     hpn: grid number in each element
     gllp:
     gllw:
     Uni_Je
Output: stiffness matrix of 3 components
*/

/*20140621: ZHL, modified dlag_P -> dlagP_Gllw*/
void Grid_Stiff(alpha,beta,FikM,aStiffM,uniform)
int alpha,beta;
struct FIK * FikM;
float aStiffM[3]; /*Grid point stiffness Matrix*/
int uniform;
{
    int al,be;
    int i;/*displacement components: i=0,1,2 */
    float je,item0,item1,item2;
    float Fik;/*Fik and derivative of Lagrange polynomials*/

    je = Uni_Je;
    for(i=0;i<3;i++){
        item0 = 0.; item1 = 0.; item2 = 0.;

/*        # pragma omp parallel for private (Fik) reduction (+:item0)
         20140328: if using openMP, it became very slowly */

        for(al=0;al<=hpn;al++){
            Fik = FikM[beta*(hpn+1)+al].F0[i];
//            item0 = item0 + gllw[al]*Fik*dlagP[al*(hpn+1)+alpha];
            item0 = item0 + Fik*dlagP_Gllw[al*(hpn+1)+alpha];
        }

/*2011060: adding item for y component of stiffness tensor, adding the following two lines */
/*Results: almost the same, but Yann said it is not necessary*/
//        Fik = FikM[beta*(hpn+1)+alpha].F1[i];
//        item1 = gllw[alpha]*gllw[beta]*Fik;


        for(be=0;be<=hpn;be++){
            Fik = FikM[be*(hpn+1)+alpha].F2[i];
//            item2 = item2 + gllw[be]*Fik*dlagP[be*(hpn+1)+beta];
            item2 = item2 + Fik*dlagP_Gllw[be*(hpn+1)+beta];
        }
//        aStiffM[i] = (item0 * gllw[beta] + item2 * gllw[alpha]) *je + item1;

        aStiffM[i] = (item0 * gllw[beta] + item2 * gllw[alpha]) *je;
    }
}

/*
   Calculate component Fik for (i,k) i=1,2,3; k=1,3, not 2
   alpha,beta: reference grid points index
   u: displacement vector
   uniform: 1: partition uniformly, 0: not
   stiff: stiffness tensor at the point 
   20110601: add Fik(k=1), test if it is correct?
*/
float F_ik(i,k,Tij,uniform)
int i,k;
float Tij[3][3];
int uniform;
{
   float Fik;
   Fik = 0.;

   if(k == 0) Fik = Tij[i][0] * Uni_JeI;
//   if(k == 1) Fik = Tij[i][1];
   if(k == 2) Fik = Tij[i][2] * Uni_JeI;

   return Fik;
}

/* Getting stress component sigma(i,j)
Require input:
   alpha,beta: reference grid points number
   i,j: component index
   ukl: displacement vector derivative: 0,1,2 component 
   stiff: stiffness tensor at the point 
*/
float Stress_ij(i,j,ukl,stiff,uniform)
int i,j;
float ukl[3][3];
struct cmatrix stiff;
int uniform;
{
   float stress=0.;
   int k,l; 
/* 20140622: Zhl: testing comment the following lines to increase speed
   for(k=0;k<3;k++){
      for(l=0;l<3;l++){
            stress += stiff.c[i][j][k][l]*(ukl[k][l]+ukl[l][k])/2;
      }
   }
*/

/*
      stress  = stiff.c[i][j][0][0]*ukl[0][0];
      stress += stiff.c[i][j][2][2]*ukl[2][2];
      stress += stiff.c[i][j][1][2]*ukl[1][2];
      stress += stiff.c[i][j][0][2]*(ukl[2][0]+ukl[0][2]);
      stress += stiff.c[i][j][0][1]*ukl[1][0];
*/
   if(i == 1 && j == 1) stress = 0.;
   else{
      stress  = stiff.c[i][j][0][0]* ukl[0][0];
      stress += stiff.c[i][j][0][1]*(ukl[0][1]+ukl[1][0]);
      stress += stiff.c[i][j][0][2]*(ukl[2][0]+ukl[0][2]);
      stress += stiff.c[i][j][1][1]* ukl[1][1];
      stress += stiff.c[i][j][1][2]*(ukl[1][2]+ukl[2][1]);
      stress += stiff.c[i][j][2][2]* ukl[2][2];
   }

   return stress;
}

/* Refer to Komatitsch and Tromp, GJI, 1999, P821, equation A2
   Calculate partial derivative of Uj WRT xi in a element (ex,ez)
Require input:
   ex,ez: elemental number horizontal and vertical
   alpha,beta: coordinates in the reference framework using GLL points, in rx, rz directions
   uj: displacement of component j
Globals: 
   gllp: gll points
   hpn
   Inverse Jacabian matrix
*/

/*20140621: ZHL, modified dlag_P -> dlagP_JeI*/
float dUj_di(ex,ez,alpha,beta,p,j,i,uniform)
int ex,ez;
int alpha,beta;
struct vector * p;
int j,i;/*component partial derivative i=0: x; i=1: z*/
int uniform; /*Partition uniformly or not*/
{
    int index,theta;
    float uji=0.;
    float item;
    item = 0.; 

       /*partial derivative wrt x: i=0*/
       if(i==0) {
          /*20140328, using openMP, works */
//          #pragma omp parallel for private (index) reduction (+:item)
          for(theta=0;theta<=hpn;theta++){
              index = (ez*hpn+beta)*nx + ex*hpn + theta;
              item  = item + p[index].u[j] * dlagP_JeI[alpha*(hpn+1)+theta];

          }
//          uji = item * Uni_JeI;
       }
       /*partial derivative wrt z: i=2*/
       if(i==2) {
          /*20140328, using openMP, works*/
//          #pragma omp parallel for private (index) reduction (+:item)

          for(theta=0;theta<=hpn;theta++){
              index = (ez*hpn+theta)*nx + ex*hpn + alpha;
              item  = item + p[index].u[j] * dlagP_JeI[beta*(hpn+1)+theta];
          }
//          uji = item * Uni_JeI;
       }
    return uji;
}

/* Refer to Komatitsch and Tromp, GJI, 1999, P821, equation A2
   Difference between dUj_di: calculate dUj/dxi[][] at the same time
   Calculate partial derivative of Uj WRT xi in a element (ex,ez)
Require input:
   ex,ez: elemental number horizontal and vertical
   alpha,beta: coordinates in the reference framework using GLL points, in rx, rz directions
   uj: displacement of component j
Globals: 
   gllp: gll points
   hpn
   Inverse Jacabian matrix
*/
/*20140621: ZHL, modified dlag_P -> dlagP_JeI*/
void dUj_di_2D(ex,ez,alpha,beta,p,duji,uniform)
int ex,ez;
int alpha,beta;
struct vector * p;
float  duji[3][3];
int uniform; /*Partition uniformly or not*/
{
    int index,theta;
    int i,j;
    for(j=0;j<3;j++){
       for(i=0;i<3;i++) duji[j][i] = 0.;
    }

    /*partial derivative wrt x: i=0*/

    /*20140328, using openMP, worked but used 2 times longer time */
    /* tmp0 =0; tmp1 = 0; tmp2=0;
    #pragma omp parallel for private (index) reduction (+:tmp0) reduction(+:tmp1) reduction(+:tmp2)*/

    for(theta=0;theta<=hpn;theta++){
         index = (ez*hpn+beta)*nx + ex*hpn + theta;
         duji[0][0]  = duji[0][0] + p[index].u[0] * dlagP_JeI[alpha*(hpn+1)+theta];
         duji[1][0]  = duji[1][0] + p[index].u[1] * dlagP_JeI[alpha*(hpn+1)+theta];
         duji[2][0]  = duji[2][0] + p[index].u[2] * dlagP_JeI[alpha*(hpn+1)+theta];
    }
    /*partial derivative wrt z: i=2*/
    for(theta=0;theta<=hpn;theta++){
         index = (ez*hpn+theta)*nx + ex*hpn + alpha;
         duji[0][2]  = duji[0][2] + p[index].u[0] * dlagP_JeI[beta*(hpn+1)+theta];
         duji[1][2]  = duji[1][2] + p[index].u[1] * dlagP_JeI[beta*(hpn+1)+theta];
         duji[2][2]  = duji[2][2] + p[index].u[2] * dlagP_JeI[beta*(hpn+1)+theta];
    }
/*  20140620 used dlagP_JeI 
  for(j=0;j<3;j++){
         for(i=0;i<3;i++) duji[j][i] *= Uni_JeI;
    }
*/

}

/* Get stiffness coefficient from tensor of order 4
Input: i, j,k,l
Require: apoint: storing the stifness property at a grid points with index index
Output: Cijkl
*/
float Cijkl(i,j,k,l,apoint)
int i,j,k,l;
struct cmatrix2d apoint;
{
   float cmn;/*Using Musgrave brevity, 1970*/
   int m,n;
   cmn=0.; 
   m = inv_index(i,j);
   n = inv_index(k,l);
   if(m==0 && n ==0) {cmn=apoint.c11; }
   if(m==0 && n ==2) {cmn=apoint.c13; }
   if(m==0 && n ==5) {cmn=apoint.c16; }
   if(m==2 && n ==0) {cmn=apoint.c13; }
   if(m==2 && n ==2) {cmn=apoint.c33; }
   if(m==2 && n ==5) {cmn=apoint.c36; }
   if(m==3 && n ==3) {cmn=apoint.c44; }
   if(m==3 && n ==4) {cmn=apoint.c45; }
   if(m==4 && n ==3) {cmn=apoint.c45; }
   if(m==4 && n ==4) {cmn=apoint.c55; }
   if(m==5 && n ==0) {cmn=apoint.c16; }
   if(m==5 && n ==2) {cmn=apoint.c36; }
   if(m==5 && n ==5) {cmn=apoint.c66; }
   return cmn;
}

/* Get stiffness coefficient from tensor of order 4
Input: i, j,k,l
Require: apoint: storing the stifness property at a grid points with index index
Output: Cijkl
     0   1   2   3   4   5
    _______________________
0  | 11	 12  13  14  15  16
1  | 21	 22  23  24  25  26
2  | 31	 32  33  34  35  36
3  | 41	 42  43  44  45  46
4  | 51	 52  53  54  55  56
5  | 61	 62  63  64  65  66

*/
float Cijkl_3D(i,j,k,l,apoint)
int i,j,k,l;
struct cmatrix3d apoint;
{
   float cmn;/*Using Musgrave brevity, 1970*/
   int m,n;
   cmn=0.; 
   m = inv_index(i,j);
   n = inv_index(k,l);
   if(m==0 && n ==0) cmn=apoint.c11; 
   if(m==0 && n ==1) cmn=apoint.c12; 
   if(m==0 && n ==2) cmn=apoint.c13; 
   if(m==0 && n ==5) cmn=apoint.c16; 

   if(m==1 && n ==0) cmn=apoint.c12; 
   if(m==1 && n ==1) cmn=apoint.c22; 
   if(m==1 && n ==2) cmn=apoint.c23; 
   if(m==1 && n ==5) cmn=apoint.c26; 

   if(m==2 && n ==0) cmn=apoint.c13; 
   if(m==2 && n ==1) cmn=apoint.c23; 
   if(m==2 && n ==2) cmn=apoint.c33; 
   if(m==2 && n ==5) cmn=apoint.c36; 

   if(m==3 && n ==3) cmn=apoint.c44; 
   if(m==3 && n ==4) cmn=apoint.c45; 

   if(m==4 && n ==3) cmn=apoint.c45; 
   if(m==4 && n ==4) cmn=apoint.c55; 

   if(m==5 && n ==0) cmn=apoint.c16; 
   if(m==5 && n ==1) cmn=apoint.c26; 
   if(m==5 && n ==2) cmn=apoint.c36; 
   if(m==5 && n ==5) cmn=apoint.c66; 
   return cmn;
}

/*
   Get index in Musgrave breviety
 Musgrave index 
 in c language array
 0  1  2  3  4  5
 00 11 22 12 02 01
 or in fortran
 1  2  3  4  5  6
 11 22 33 23 13 12
*/
int inv_index(i,j)
int i,j;
{
   int op;
   if(i==0 && j==0) op = 0;
   if(i==0 && j==1) op = 5;
   if(i==0 && j==2) op = 4;
   if(i==1 && j==0) op = 5;
   if(i==1 && j==1) op = 1;
   if(i==1 && j==2) op = 3;
   if(i==2 && j==0) op = 4;
   if(i==2 && j==1) op = 3;
   if(i==2 && j==2) op = 2;
   return op;
}

/*
  Initialize UVW vectors for all the grids
  uvw: displacement vector
Require:
  Global: nx, nz: grid number horizontal and vertical
*/
void Initial_uvw(p)
struct vector * p;
{
   int i,j;
   for(i=0;i<=nx*nz;i++){
       for(j=0;j<3;j++){
           p[i].u[j] = 0.;
       }
   }
}

/*
  Copy uvw field to uvw_l
  uvw: displacement at time t
  uvwl: displacement at time t-dt
Requires:
  nx,nz
*/
void Copy_Displacement(p,pl)
struct vector *p;
struct vector *pl;
{
   int i,j;
   int iex,iez;

   for(j=0;j<3;j++) {
      /*20140329, using openMP works*/
//      #pragma omp parallel for 
      for(i=0;i< nx*nz;i++){
            pl[i].u[j] = p[i].u[j];
/*
              if(!( pl[i].u[j]<2 && pl[i].u[j] > -2)) {
                   iez = (int) (i/nx);
                    iex = i - iez*nx;
                    printf("%d %d %d %e\n",iex,iez,j,pl[i].u[j]); 
                    break;
                  }
*/
       }
   }
}

/*
Using explosive source at the left-bottom corner of the (mex,mez) element
Input:
   md: moment density at time t using Ricker[it]
   force: force component at the element
Note: souce is deployed on the (0,0) point of the element 
*/
void Explosive(mex,mez,md,Force)
int mex, mez;
float md;
struct vector* Force;
{
    int index,i,j,k,alpha,beta;
    float gik[3][3];
    float mij,je,jei,dLag;
    struct vector force;
    force = * Force;

    jei = Uni_JeI;
    je  = Uni_Je;

    for(i=0;i<3;i++){
       for(k=0;k<3;k++){
           gik[i][k] = 0.0;
           for(j=0;j<3;j++){
               if(j !=1 ){
                   if(i == j) mij = md;
                     else mij = 0;/*explosive source*/
                   if(k == j) gik[i][k] = gik[i][k] + mij*jei;
               }
           }
       }
    } 
    for(i=0;i<3;i++){
       force.u[i] = 0;
       for(k=0;k<1;k++){ /*if k<4, need to note the signal + or - at different directions*/

/*shared by four adjacent elements
|---------|---------|
|         |(mex,mez)|
|   2     |    3    |
|---------|---------|
|         |         |
|   k=0   |    1    |
|---------|---------|
*/
          if(k==0) {alpha = hpn; beta= hpn;}
          if(k==1) {alpha = 0;   beta= hpn;}
          if(k==2) {alpha = hpn; beta= 0;}
          if(k==3) {alpha = 0;   beta= 0;}
          dLag= dlagP[alpha*(hpn+1)+alpha];
          force.u[i] = force.u[i]+ gllw[beta]*gllw[alpha]*je*gik[i][0]*dLag;

          dLag= dlagP[beta*(hpn+1)+beta];
          force.u[i] = force.u[i] + gllw[alpha]*gllw[beta]*je*gik[i][2]*dLag;
       }
    }
    * Force = force;
}

/* Output a snapshot at time t
*/
void snapshot(file,uvw,nx0,nz0,increx,increz)
char *file; 
struct vector * uvw;
int nx0,nz0; 
int increx,increz;
{
    int ix, iz, index;
    FILE *fop, *open_file();

    fop=open_file(file,"wt");
    for(ix=0; ix< nx0; ix += increx)
        for(iz=0; iz< nz0; iz +=increz){
           index = iz*nx0 + ix;
           fprintf(fop,"%d %d %e %e\n",ix,iz,uvw[index].u[0], uvw[index].u[1]);
        }
    fclose(fop);        
}    


/*For attenuating the boundary*/
/*
Comment Mar, 18, 2008: it seems that attenuation is not fit for SEM
*/
void atten_taper(nh0,ass,taper)
int nh0;
float ass, *taper;
/* REQUIRES: ass Positive                                      */
/* MODIFIES: taper[0..nh-1]                                    */
/* ENSURES:  taper[0]=1, taper[nh-1]=ass, taper exp function   */
{
    int n;
    double alpha;

    alpha=1.0/nh0*log10((double) (ass))/log10(2.73);
    
    for(n=0; n< nh0; n++)
        taper[n]=exp((double) (alpha*n));
}

void atten_wave(p,nh0,atten,taper)
struct vector * p;
int nh0,atten;
float *taper;
{
  int ix, iz, i,offset;
  float *ptaper;
  
  for (ix=0; ix<nx; ix++)
      for (iz=0; iz< nz; iz++){ 
           /*bottom boundary*/
           if(atten>0 && iz< nh0){
               ptaper = taper + nh0- iz -1;
               offset = iz*nx+ix;
               for(i=0;i<3;i++) p[offset].u[i] *= ptaper[0];
           }

/* 071409, LZ comment attenuation at the left and top boundary */
           /*bottom + left boundaries*/
           if(atten>1 && ix< nh0){
               ptaper = taper + nh0 - ix -1;
               offset = iz*nx + ix;
              for(i=0;i<3;i++) p[offset].u[i] *= ptaper[0];
         }
           /*bottom left + top boundaries*/
           if(atten>2 && iz>= nz-nh0){
               ptaper = taper + iz - nz + nh0;
               offset = iz*nx+ix;
               for(i=0;i<3;i++) p[offset].u[i] *= ptaper[0];
           }
          /*bottom left top + right boundaries*/
           if(atten>3 && ix >= nx-nh0){
               ptaper = taper + ix - nx + nh0;
               offset = iz*nx+ix;
               for(i=0;i<3;i++) p[offset].u[i] *= ptaper[0];
           }
       }
}

/*20140220. ZHL: added atten_wave2() to avoid multiple-time attenuation at several corners*/
void atten_wave2(p,nh0,atten,taper)
struct vector * p;
int nh0,atten;
float *taper;
{
  int ix, iz, i,offset;
  float *ptaper;

  for(i=0;i<3;i++){  
    for (ix=0; ix<nx; ix++){

        /*20140329, using openMP works*/
//       #pragma omp parallel for private (ptaper,offset)
        for (iz=0; iz< nz; iz++){ 
           /*bottom boundary*/
            if(atten>0 && iz< nh0){
               ptaper = taper + nh0- iz -1;
               offset = iz*nx+ix;
               p[offset].u[i] *= ptaper[0];
            }

           /*bottom + left boundaries*/
            if(atten>1 && ix< nh0){
               ptaper = taper + nh0 - ix -1;
               offset = iz*nx + ix;
               if(iz>=nh0) p[offset].u[i] *= ptaper[0];
            }
           /*bottom left + right boundaries*/
            if(atten>2 && ix >= nx-nh0){
               ptaper = taper + ix - nx + nh0;
               offset = iz*nx+ix;
               if(iz>=nh0) p[offset].u[i] *= ptaper[0];
            }

           /*bottom left + right boundaries*/
            if(atten>3 && iz>= nz-nh0){
               ptaper = taper + iz - nz + nh0;
               offset = iz*nx+ix;
               if(ix<nx-nh0 && ix>=nh0)  p[offset].u[i] *= ptaper[0];
            }
        }
    }
  }
}
/*
   Apply absorbing boundary condition
Require:
  uvw1,uvw,uvw0: displacement vector at time t+dt, t, and t-dt
  GStiff: global stiffness matrix
  Mass: global mass matrix
  sdt: time step
  Global: nx,nz
         note that ab_matrix has multiplied sdt/2
*/
void Absorbs(uvw1,uvw,uvw0,GStiff,Mass,sdt,ab_matrix)
struct vector * uvw1,*uvw,*uvw0;
struct vector * GStiff;
float *Mass, sdt;
struct vector * ab_matrix;
{
    int i,gn,k;
    float dt2,coe,coe1;

    dt2 = sdt *sdt;

    /*left boundary*/
    for(i=0;i<nz;i++){
        gn = i*nx ;
        for(k=0;k<3;k++){
             coe  = Mass[gn] + ab_matrix[gn].u[k];
             coe1 = Mass[gn] - ab_matrix[gn].u[k];
             uvw1[gn].u[k] = (-dt2*GStiff[gn].u[k] + 2.*Mass[gn]*uvw[gn].u[k] - coe1*uvw0[gn].u[k])/coe;
        }
    }
    /*bottom boundary*/
    for(i=0;i<nx;i++){
        gn = i ;
        for(k=0;k<3;k++){
             coe  = Mass[gn] + ab_matrix[gn].u[k];
             coe1 = Mass[gn] - ab_matrix[gn].u[k];
             uvw1[gn].u[k] = (-dt2*GStiff[gn].u[k] + 2.*Mass[gn]*uvw[gn].u[k] - coe1*uvw0[gn].u[k])/coe;
        }
    }
    /*right boundary*/
    for(i=0;i<nz;i++){
        gn = i*nx + nx -1;
        for(k=0;k<3;k++){
             coe  = Mass[gn] + ab_matrix[gn].u[k];
             coe1 = Mass[gn] - ab_matrix[gn].u[k];
             uvw1[gn].u[k] = (-dt2*GStiff[gn].u[k] + 2.*Mass[gn]*uvw[gn].u[k] - coe1*uvw0[gn].u[k])/coe;
        }
    }
    /*top boundary*/
/* July 06, 2011, added top boundary absorbing, the waveform is clearer, and the corner effects are reduced
   but the reflection is weaker!
   but Yann: not necessary, because it is the free surface boundary
   20131227: only absorb several elements
*/  
    int ntop;
    ntop = 3;
    for(i=1;i<ntop*hpn+1;i++){
        gn = (nz-1)*nx + i ;
        for(k=0;k<3;k++){
             coe  = Mass[gn] + ab_matrix[gn].u[k];
             coe1 = Mass[gn] - ab_matrix[gn].u[k];
             uvw1[gn].u[k] = (-dt2*GStiff[gn].u[k] + 2.*Mass[gn]*uvw[gn].u[k] - coe1*uvw0[gn].u[k])/coe;
        }
    }

}

/*
   Apply absorbing boundary condition same as that in FD
Require:
  uvw: displacement vector at time  t
  meshs
  Global: nx,nz

*/
void Absorbs2(p,meshs,dt)
struct vector *p ;
struct cmatrix * meshs;
float dt;
{
    int i,gn,k;
    struct cmatrix mesh;
    float fac;
    int ntop;
    fac = dt/h;

    /*left boundary*/
    for(i=0;i<nz;i++){
        gn = i*nx ;
        mesh = meshs[gn+1];
        for(k=0;k<3;k++){
             p[gn].u[k] -= p[gn + 1].u[k]/mesh.rho*fac; 
        }
    }
    /*bottom boundary*/
    for(i=1;i<nx;i++){
        gn = i ;
        mesh = meshs[gn+nx];
        for(k=0;k<3;k++){
             p[gn].u[k] -= p[gn + nx].u[k]/mesh.rho*fac;
        }
    }
    /*right boundary*/
    for(i=0;i<nz;i++){
        gn = i*nx + nx -1;
        mesh = meshs[gn-1];
        for(k=0;k<3;k++){
             p[gn].u[k] -= p[gn -1].u[k]/mesh.rho*fac;
        }
    }
    /*top boundary*/
    ntop = 3; /**/
    for(i=0;i<ntop*hpn+1;i++){
        gn = (nz-1)*nx + i;
        mesh = meshs[gn-nx];
        for(k=0;k<3;k++){
             p[gn].u[k] -= p[gn - nx].u[k]/mesh.rho*fac;
        }
    }

}

/* subroutine for GRT-SEM interfaces*/
/*
   Requires: uvw: displacement matrix
             grtuw: u and w at the GRT-SEM boundaries
   Method: 
*/

void GRT_Interface(p,grtuw,grtv)
struct vector * p;
float *grtuw,*grtv;
{
   int i,m;
   int index,idz;
/* ----------------------------------------------------------------*/
/*05252011: Tried only couple at the elements' anchor points*/
/*Result: it seems give much clearer waveform */
/*05262011: tried only couple horizotal component, i.e., u[0] but not u[1]*/
/*Result: S-> P phases is much weaker, but still obviously*/
/*05262011: Force the vertical component to be zero */
/*05302011: assemble elements by adding */
/*20110601: change the interface to the second elements to the left and bottom boundaries*/
/*nbl is the number of the grid number of GRT-SEM interface to the boundaries*/
/*20110602: change the sign of u[3] component */
/*20141028: extended to u,v,w imping at the same time*/
   int nbl;
   nbl = neh*hpn;

   for(i=0;i<nx;i++){
      index = nbl*nx;
      p[index+ i].u[0] = 0;
      p[index+ i].u[1] = 0;
      p[index+ i].u[2] = 0;
   }
   for(i=0;i<nz;i++){
      index= i*nx+nbl;
      p[index].u[0] = 0;
      p[index].u[1] = 0;
      p[index].u[2] = 0;
   }
// Bottom boundary
   for(i=neh;i<nex;i++){
       for(m=0;m<=hpn;m++){
           index = nbl*nx+ i*hpn + m;
           p[index].u[0] =  grtuw[i*hpn+m];
           p[index].u[1] =  grtv[i*hpn+m];
           p[index].u[2] = -grtuw[nx+i*hpn+m];
       }
   }
// Left boundary
   for(i=neh; i<nez;i++){
       for(m=0;m<=hpn;m++){
           idz   = i*hpn+m;
           index = idz * nx + nbl;
           p[index].u[0] = grtuw[2*nx+idz];
           p[index].u[1] = grtv[nx+idz];
           p[index].u[2] = -grtuw[2*nx+nz+idz];
       }
   }
/* -----------------------------------------------------------------*/
}

/* subroutine for GRT-SEM interfaces taking account into reflection and scattering*/
/*
   Requires: uvw: displacement matrix
             grtuw: u and w at the GRT-SEM boundaries
             grtuw0: u and w at the GRT-SEM interface, previous time step
   set Total = Incident + refelction
   then T(n+1) = I(n+1) + R(n) and R(n) = T(n) - I(n)
      -> T(n+1) = I(n+1) + {T(n) - I(n)}
   keep in mind: for vertical component u[2]: GRT and SEM has different sign
   Method: 
*/

void GRT_Interface2(p,grtuw,grtuw0,grtv,grtv0)
struct vector * p;
float *grtuw, *grtuw0, *grtv, *grtv0;
{
   int i,m;
   int index,idz;
   int nbl;
   nbl = neh * hpn;

   /*bottom and left boudaries*/
/* for vertical component u[2]: GRT and SEM has different sign */
   for(i=nbl;i<nx;i++){
      index = nbl*nx + i;
      p[index].u[0] = + grtuw[i];
      p[index].u[1] = + grtv[i];
      p[index].u[2] = - grtuw[i+nx];
   }
   for(i=nbl;i<nz;i++){
      index= i*nx+nbl;
      p[index].u[0] = + grtuw[2*nx+i];
      p[index].u[1] = + grtv[nx+i];
      p[index].u[2] = - grtuw[2*nx+nz+i];
   }
}


void GRT_Interface3(p,p0,grtuw,grtuw0,grtv,grtv0)
struct vector * p,*p0;
float *grtuw, *grtuw0, *grtv, *grtv0;
{
   int i,m;
   int index,idz;
   int nbl;
   /*bottom and left boudaries*/
/* 
     u(scatter) = u(t)-u(t-1)
     u(t+1) = u(scatter) + u(incident)
*/
/* for vertical component u[2]: GRT and SEM has different sign */
   int istart;
   istart = neh * hpn;
   nbl = neh * hpn;
   for(i=istart;i<nx;i++){
      index = nbl *nx + i;
      p[index].u[0] += - p0[index].u[0] + grtuw[i];
      p[index].u[1] += - p0[index].u[1] + grtv[i];
      p[index].u[2] += - p0[index].u[2] - grtuw[i+nx];
   }
   istart = neh * hpn + 1;
   for(i=istart;i<nz;i++){
      index= i*nx + nbl;
      p[index].u[0] += - p0[index].u[0] + grtuw[2*nx+i];
      p[index].u[1] += - p0[index].u[1] + grtv[nx+i];
      p[index].u[2] += - p0[index].u[2] - grtuw[2*nx+nz+i];
   }
}

void GRT_Interface4(p,p0,grtuw,grtuw0,grtv,grtv0)
struct vector * p,*p0;
float *grtuw, *grtuw0, *grtv, *grtv0;
{
   int i,m;
   int index,idz;
   int nbl;
   /*bottom and left boudaries*/
/* 
     u(scatter) = u(t)-u(t-1)
     u(t+1) = u(scatter) + u(incident)
*/
/* for vertical component u[2]: GRT and SEM has different sign */
   int istart;
   istart = neh * hpn;
   nbl = neh * hpn;
   for(i=istart;i<nx;i++){
      index = nbl *nx + i;
      p[index].u[0] += - grtuw0[i]    + grtuw[i];
      p[index].u[1] += - grtv0[i]     + grtv[i];
      p[index].u[2] +=   grtuw0[i+nx] - grtuw[i+nx];
   }
   istart = neh * hpn + 1;
   for(i=istart;i<nz;i++){
      index= i*nx + nbl;
      p[index].u[0] += - grtuw0[2*nx+i]  + grtuw[2*nx+i];
      p[index].u[1] += - grtv0[nx+i]     + grtv[nx+i];
      p[index].u[2] += grtuw0[2*nx+nz+i] - grtuw[2*nx+nz+i];
   }
}


/*
  Output time-sequence waveform at the receivers
Requires:
  p: uvw wave field
  ix1,iz1: receiver's coordinate
  xtrace, ztrace: trace number
  fu, fw, fv: file pointer
*/
void kirrecord_zhl(p,iz1,ix1,xtrace,ztrace,fu,fw,fv, source)
struct vector * p;
int  *iz1, *ix1, xtrace, ztrace;
FILE *fu, *fw, *fv;
int source;
/* source = 0   -> velocity Vx, Vy and Vz are recorded */ 
/* source = 1   -> P wave and S wave (grad div u) and (curl curl u) are recorded */ 
/* source = 2   -> both velocity and P and S waves are recorded */ 

{
    register int i;
    register struct vector *pp;
    register float *hold;
    int ntrace, offset;

    ntrace = xtrace*ztrace;
    hold= (float *) malloc(sizeof(float)*(ntrace));

    if(source==0 || source <=-2){
        /*20140329, using openMP works*/
//        #pragma omp parallel for private (offset,pp)
	for(i=0; i< ntrace; i++)		/* u component */
	{
	    offset   = iz1[i]*nx + ix1[i];
	    pp       = p      + offset;
            hold[i]  =  pp[0].u[0];
	}
	fwrite(hold,sizeof(float),ntrace,fu);
        /*20140329, using openMP works*/
//        #pragma omp parallel for private (offset,pp)
	for(i=0; i< ntrace; i++)		/* w component */
/*20110607: change the direction of u[2] */
	{
	    offset   = iz1[i]*nx + ix1[i];
	    pp       = p      + offset;
            hold[i]  = - pp[0].u[2];
	}
	fwrite(hold,sizeof(float),ntrace,fw);
        /*20140329, using openMP works*/
//        #pragma omp parallel for private (offset,pp)

	for(i=0; i< ntrace; i++)		/* v component */
	{
	    offset   = iz1[i]*nx + ix1[i];
	    pp       = p      + offset;
            hold[i]  =  pp[0].u[1];
	}
	fwrite(hold,sizeof(float),ntrace,fv);
    }	
}
