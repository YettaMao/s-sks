/*********************************************************************/
/*                                                                   */
/*          C code for read the results from sem2d                   */
/*                                                                   */
/*  Written by: Liang Zhao, IGG, CAS                                 */
/*  Last modified: Jan. 31, 2008.                                    */
/*                                                                   */
/*********************************************************************/

/* modified history:
*/

/* read the output of FD calculation, refer to sem2d.c
   Usage:  ./read_uvw par=*.par
   Input:  get par from *.par
           seis.u seis.v seis.w
   Output: log.txt u.xy v.xy w.xy 
*/

#include   <stdio.h>
#include   <stdlib.h>
#include   <sys/file.h>
#include   <fcntl.h>
#include   <math.h>
#include   <malloc.h>

main(ac,av)
int ac; 
char **av;
{
   int i,k,j,num;
   int ix1[700], iz1[700], itrace;
   int ixbeg, izbeg,idxplot, idzplot;
   int xtrace = 700, ztrace = 1;
   int nex, nez,nx,nz, nt,nst,hh,source,len;
   int hpn;
   float step = 0.2; /* time step for output */
   int kdt, istep;
   int flat = 1,nfil=1;

   float xmin,xbeg,zbeg, tstart,gcarc,zs ;
   float eh,h,dp,dt,time,ndt;
   FILE * rfp; 
   float rec_x;
   int  ixrec; 
   int get_rec; /* 1: assign receiver's position based on raw data; 0: not ;  */
   get_rec = 0;  /* assgin receiver's position ;*/
   char raymodel[80];   /* file for model in GRT calculation */
   char rec_file[80]; /* assgin receiver's position ; */

   FILE *fopu1, *fopw1, *fopv1; 
   FILE *fopu2, *fopw2, *fopv2, *fp, *fop; 

   kdt = 1;
   setpar(ac,av);
   mstpar("raymodel",    "s", raymodel);
   mstpar("nex",         "d", &nex);
   mstpar("nst",         "d", &nst); /* nst for SEM calculation */
   mstpar("hpn",         "d", &hpn);
   mstpar("eh",          "f", &eh);
   mstpar("dp",          "f", &dp);  /* dt for GRT calculation */
   dt = dp;
   getpar("sdt",         "f", &dt);  /* dt for the SEM */
   getpar("nfil",        "d", &nfil);
   getpar("kdt",         "d", &kdt);
   getpar("source",      "d", &source);

   getpar("xmin",        "f", &xmin);
   getpar("tstart",      "f", &tstart);
   getpar("xtrace",      "d", &xtrace);
   getpar("ztrace",      "d", &ztrace);
   mstpar("zbeg",        "f", &zbeg);
   mstpar("idzplot",     "d", &idzplot);
   mstpar("xbeg",        "f", &xbeg);
   mstpar("idxplot",     "d", &idxplot);
   mstpar("zs",          "f", &zs); /* event depth*/

   getpar("get_rec",     "d", &get_rec);
   printf("get_rec=%d\n", get_rec);
   if(get_rec == 1){
      mstpar("rec_file",   "s", rec_file);
      printf("Assign receiver's positions by %s\n",rec_file);
   }
   endpar();
   
/*read 1-D model */
   int   jo, nb; 
   int   lfinal, *nen, *ncoun, *nna, *nray;
   float *cc, *ss,*dd, *tth; /*vp vs density and thickness of 1-D model*/
   int   firstlayer, nrec; /*nrec: the last layer of GRT calculation */
   float thickness;

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
   h = eh/hpn;
   printf("getrec=%d %f h= %f\n",get_rec, eh,eh/hpn);

   /*zhl: for reading the record */   
   nx = nex*hpn+1;
   nz = nez*hpn+1;
   ixbeg = (int) (xbeg/h);
   izbeg = (int) (zbeg/h);
   izbeg += nz + 1;
   if(izbeg > nz-1) izbeg=nz-1;
   printf("nst= %d ixbeg= %d izbeg= %d \n", nst, ixbeg, izbeg);
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
          fscanf(rfp,  "%d", &xtrace);
          for(itrace=0; itrace<xtrace; itrace++){
              fscanf(rfp,"%f", &rec_x);
              ixrec = (int) ((rec_x*111.2 - xmin)/h);
              ix1[itrace] = ixrec;
              iz1[itrace] = izbeg;
          }
         fclose(rfp);
       }
   
   istep = (int) (step/dt/kdt);
   nt = (int) (nst/kdt);
   printf("nst=%-5d nt=%-5d istep=%-5d\n",nst,nt,istep);
   if(istep < 1) istep = 1;
   fopu1 = fopen("seis.u","rb");     
   fopv1 = fopen("seis.v","rb");     
   fopw1 = fopen("seis.w","rb");     

   int it,ntrace,offset;
   ntrace = xtrace*ztrace;
   len = (int) (nst*ntrace/istep);
   float *outu, *outv, *outw;
   float enlarge;
   enlarge = 2;
   
   float *holdu, *holdv, *holdw;
   holdu = (float *) malloc(sizeof(float)*ntrace);
   holdv = (float *) malloc(sizeof(float)*ntrace);
   holdw = (float *) malloc(sizeof(float)*ntrace);

   outu = (float *) malloc(sizeof(float)*len);
   outv = (float *) malloc(sizeof(float)*len);
   outw = (float *) malloc(sizeof(float)*len);

   num = 0;
   printf("\n Now, read the results ...\n");
   printf("xtrace=%d ztrace=%d\n",xtrace, ztrace);
   printf("nt=%d dt=%-10.5f istep=%d \n",nt, dt, istep);
   float maxv_u,maxv_v,maxv_w;
   maxv_u = -999;
   maxv_v = -999;
   maxv_w = -999;
   for(it=0; it<nt; it++)
   {
       time = it * dt *kdt;
       fread(holdu,sizeof(float),ntrace,fopu1);
       fread(holdv,sizeof(float),ntrace,fopv1);
       fread(holdw,sizeof(float),ntrace,fopw1);
       if(!(it%istep)) 
       {
          for(itrace=0; itrace< xtrace; itrace ++)
          {
              for(j=0; j< ztrace; j++)
              {
                  k = j*xtrace + itrace;
                  offset = num*ntrace + k;
                  outu[offset] = holdu[k];
                  outv[offset] = holdv[k];
                  outw[offset] = holdw[k];
                  if(fabs(holdu[k]) > maxv_u) maxv_u = fabs(holdu[k]);
                  if(fabs(holdv[k]) > maxv_v) maxv_v = fabs(holdv[k]);
                  if(fabs(holdw[k]) > maxv_w) maxv_w = fabs(holdw[k]);
              }
          }
          num ++;
       }
       if(!(it%200)) printf("it= %-5d time= %-8.3f sec: %-10.3e %-10.3e %-10.3e\n", it, time,holdu[0],holdv[0],holdw[0]);
   }
   
   fclose(fopu1);
   fclose(fopv1);
   fclose(fopw1);

   int t_m;
   t_m = (int)(it/istep);
   printf("time number= %d it/istep= %d \n", num, t_m);
   if(num > t_m) num--;
   printf("\n Now, write t & amplitude pairs, tstart=%-11.5f\n", tstart);

   float t_e, x_b, x_e;
   float scale_u, scale_v, scale_w,unit;
   float dis2deg= 111.2;
   t_e = tstart + num*istep*dt*kdt;
   x_b = (xmin + ixbeg )/ dis2deg;
   x_e = (xmin + ixbeg + idxplot*(xtrace-1)*h)/ dis2deg;
   unit = idxplot/dis2deg;
   if(maxv_u == 0) maxv_u = 1e-10;
   if(maxv_v == 0) maxv_v = 1e-10;
   if(maxv_w == 0) maxv_w = 1e-10;
   scale_u = unit/maxv_u*enlarge;
   scale_v = unit/maxv_v*enlarge;
   scale_w = unit/maxv_w*enlarge;
   fp   = fopen("log.txt","wt");
   fprintf(fp,"%-9.3f %-9.3f\n", tstart,t_e);
   fprintf(fp,"%-9.3f %-9.3f\n", x_b,x_e);
   fclose(fp);

   fopu2 = fopen("u.xy","wt");     
   fopv2 = fopen("v.xy","wt");     
   fopw2 = fopen("w.xy","wt");     
   
   for(itrace=0; itrace< xtrace; itrace ++)
   {
//       gcarc = (xmin + ixbeg*h + idxplot*itrace*h)/ dis2deg;
       gcarc = (xmin + xbeg + idxplot*itrace*h)/ dis2deg;
       for(j=0; j< ztrace; j++)
       {
           if(itrace==0 && j==0) printf("No.1 receiver at %-8.3f deg\n",gcarc);
           fprintf(fopu2,">> %d\n", num);
           fprintf(fopv2,">> %d\n", num);
           fprintf(fopw2,">> %d\n", num);
           for(k=0; k<num; k++)
           {
              time = tstart + k*istep*dt*kdt;
              offset = k*ntrace + j*xtrace + itrace;
              fprintf(fopu2,"%-9.3f %-11.4e\n", time, gcarc + outu[offset]*scale_u);
              fprintf(fopv2,"%-9.3f %-11.4e\n", time, gcarc + outv[offset]*scale_v);
              fprintf(fopw2,"%-9.3f %-11.4e\n", time, gcarc + outw[offset]*scale_w);
           }
       }
   }
   fclose(fopu2);
   fclose(fopv2);
   fclose(fopw2);
}


