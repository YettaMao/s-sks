/*********************************************************************/
/*                                                                   */
/*          C code for read the results from sem2d.c                 */
/*                     write to sacfiles                             */
/*                                                                   */
/*  Written by: Liang Zhao, IGG, CAS                                 */
/*  Last modified: 20140514, adding codes for fixing position        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

/* modified history:
*/

/* read the output of SEM calculation, refer to sem2d.c
   Usage:  ./uvw2sac par=*.par
   Input:  get par from *.par
           seis.u seis.v seis.w
   Output: u*.sac v*.sac w*.sac 
*/
/* sac head
 *  ................-------------------------------------------
 *  |quake time     |kztime      |        |
 *  |<-   omarker ->|            |        |
 *  |<--     begin time       -->|        |
 *  |<--           end time            -->|
*/

#include "sachead.h"
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/file.h>
#include <malloc.h>

float * b, *pp;

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
   float step ; /* time step for output */
   int kdt=1, istep;
   int flat = 1,nfil=1;

   FILE *fopu1, *fopw1, *fopv1, *fop, *fppm, *fp; 

   float xmin,xbeg,zbeg, tstart,gcarc,zs ;
   float eh,h,dp,dt,time,ndt;
   FILE * rfp; 
   float rec_x,baz;
   int  ixrec; 
   int get_rec; /* 1: assign receiver's position based on raw data; 0: not ;  */
   int fix_rec; /* 1: fix receiver's longitude and latitude */
   get_rec = 0;  /* assgin receiver's position ;*/
   fix_rec = 0;
   float enlarge;
   enlarge = 4;
   baz = 270; 

   char raymodel[80];   /* file for model in GRT calculation */
   char rec_file[80];   /* assgin receiver's position ; */
   char pos_file[80];   /* fix receiver's (long, lat)*/

   void sstep();

   SACHEAD sac_h= sac_null;
/* source */
   int ite,key,idev;
   float amom,dt1, dt2, dt3, rhos;

   kdt = 1;
   idev = 0;
   setpar(ac,av);
   getpar("flat",        "d", &flat);
   mstpar("raymodel",    "s", raymodel);
   mstpar("nex",         "d", &nex);
   mstpar("nst",         "d", &nst); /* ntt for SEM calculation */
   mstpar("hpn",         "d", &hpn);
   mstpar("eh",          "f", &eh);
   mstpar("dp",          "f", &dp);  /* dt for GRT calculation */
   getpar("sdt",         "f", &dt);  /* dt for the SEM */
   step  = dt * 2;
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
   getpar("enlarge",     "f", &enlarge);
   printf("get_rec=%d\n", get_rec);
   if(get_rec == 1){
      mstpar("rec_file",   "s", rec_file);
      printf("Assign receiver's positions by %s\n",rec_file);
   }
   getpar("fix_rec",     "d", &fix_rec);
   if(fix_rec == 1){
      mstpar("pos_file",   "s", pos_file);
      mstpar("baz",        "f", &baz);
      printf("fix receiver's positions by %s\n",rec_file);
   }

   mstpar("Conv_Source", "d", &ite);
   if(ite == 1){
	  mstpar("Moment",        "f", &amom);
	  mstpar("Dt1",           "f", &dt1);
	  mstpar("Dt2",           "f", &dt2);
	  mstpar("Dt3",           "f", &dt3);
	  mstpar("beta",          "f", &rhos);
	  getpar("point_source",  "d", &key);
	  getpar("finaldev",      "d", &idev);
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
/*              ixrec = (int) ((rec_x*111.2 - xmin)/h);*/
              ixrec = (int) (rec_x/h);
              ix1[itrace] = ixrec;
              iz1[itrace] = izbeg;
          }
         fclose(rfp);
       }

   float longs[500], lats[500];
   int ptrace;
   char names[200][50];
   if(fix_rec == 1){ /* if fix receiver's position */
       rfp  = open_file(pos_file,"rt");
       fscanf(rfp,  "%d", &ptrace);       
       if(ptrace != xtrace) printf("ERROR: inconsistent position file and receiver file\n");
       for(j=0;j<ptrace;j++){
           fscanf(rfp,"%f %f %s", &longs[j],&lats[j],names[j]);
       }
       fclose(rfp);
   }
   
   istep = (int) (step/dt/kdt);
   nt = (int) (nst/kdt);
   printf("xtrace=%d ztrace=%d nst=%-5d nt=%-5d istep=%-5d\n",xtrace,ztrace, nst,nt,istep);
   if(istep < 1) istep = 1;

   float t_e, x_b, x_e;
   float dis2deg= 111.2;
   t_e = tstart + nt*dt*kdt;
   if(get_rec == 1){
      x_b = (xmin + ix1[0]*h )/ dis2deg;
      x_e = (xmin + ix1[xtrace-1]*h)/ dis2deg;
   }
   else{
      x_b = (xmin + xbeg )/ dis2deg;
      x_e = (xmin + xbeg + idxplot*(xtrace-1)*h)/ dis2deg;
   }
   printf("x_b=%f x_e=%f\n",x_b,x_e);
   fp   = fopen("log.txt","wt");
   fprintf(fp,"%-9.3f %-9.3f\n", tstart,t_e);
   fprintf(fp,"%-9.3f %-9.3f\n", x_b,x_e);
   fclose(fp);

   fopu1 = fopen("seis.u","rb");     
   fopv1 = fopen("seis.v","rb");     
   fopw1 = fopen("seis.w","rb");     

   int it,ntrace,offset;
   ntrace = xtrace*ztrace;
   len = (int) (nst*ntrace/istep);
   float *outu, *outv, *outw;
   
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

   float maxv = -999; 
   int uvw,t_n, t_trace;
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
                  if(outu[offset] > maxv) {maxv= outu[offset]; uvw = 1;t_n=num; t_trace=itrace;}
                  if(outv[offset] > maxv) {maxv= outv[offset]; uvw = 2;t_n=num; t_trace=itrace;}
                  if(outw[offset] > maxv) {maxv= outw[offset]; uvw = 3;t_n=num; t_trace=itrace;}
              }
          }
          num ++;
       }
   }

   printf("maxv=%-10.3e belong to uvw=%-3d t_n= %-4d iz=%-4d ix=%-4d\n",maxv,uvw,t_n,iz1[t_trace],ix1[t_trace]);
   
   fclose(fopu1);
   fclose(fopv1);
   fclose(fopw1);

   int t_m;
   t_m = (int)(it/istep);
   printf("time number= %d it/istep= %d \n", num, t_m);
   if(num > t_m) num--;
   printf("\n Now, write t & amplitude pairs, tstart=%-11.5f\n", tstart);

   int index=0;
   char u_sac[80],v_sac[80],w_sac[80];
   char pm[80];
   strcpy(pm,"pm.xy");
   char t_index[20];
   float *us, *vs, *ws;
//   float dis2deg= 111.2;

   sac_h.delta = istep*dt*kdt;
   sac_h.nzyear= 2006 ;
   sac_h.nzjday= 0 ;
   sac_h.nzhour= 0 ;
   sac_h.nzmin = 0 ; 
   sac_h.nzsec = 0 ;
   sac_h.nzmsec= 0 ;

   sac_h.evla  = 0.;

   sac_h.evlo  = 0.;
   if(fix_rec == 1){
      if(get_rec == 1) gcarc= (xmin + ix1[0]*h)/ dis2deg;
      else gcarc = (xmin + xbeg + idxplot*h)/ dis2deg;
      if(baz == 270){ /*source from east*/
             sac_h.evlo  = longs[0]-gcarc;
      }
      if(baz == 90) {
          sac_h.evlo  = longs[0]+ gcarc;
      }
      printf("baz=%-8.3f evlo=%-8.3f\n",baz,sac_h.evlo);
      printf( "fix receiver, set evlo=%f\n",sac_h.evlo);
   }
   
   sac_h.evdp  = zs;

   sac_h.stla  = 0.,
   sac_h.stel  = 0.,

   sac_h.o     = - tstart;
   sac_h.b = 0;
   sac_h.e = num*sac_h.delta;
   sac_h.npts= num;
   sac_h.iftype = 1;
   sac_h.leven = 1;
   sac_h.nvhdr = 6;

   us= (float* ) malloc(sizeof(float)*sac_h.npts);
   vs= (float* ) malloc(sizeof(float)*sac_h.npts);
   ws= (float* ) malloc(sizeof(float)*sac_h.npts);

   amom = (amom * (1.e-20))/(4. *3.1415927*rhos);
   fppm = fopen(pm,"wt");
   ndt = sac_h.delta;
   int malloced =0;

   int nfn, nfa, ifound, nnn,nfad,nnum;
/*062710 ZHL: zs has been already defined, now define again, a bug?*/
/*  float xs, zs,xr,zr,cc1, cc2,cc3, roo; */
   float xs, xr,zr,cc1, cc2,cc3, roo;

   float ne_pole = 1.;
   if(fix_rec ==1 && baz ==90) ne_pole =-1.;

   for(itrace=0; itrace< xtrace; itrace ++)
   {
       if(get_rec == 1) gcarc= (xmin + ix1[itrace]*h)/ dis2deg;
       else gcarc = (xmin + xbeg + idxplot*itrace*h)/ dis2deg;
       for(j=0; j< ztrace; j++)
       {
           if(itrace==0 && j==0) printf("No.1 receiver at %-8.3f deg\n",gcarc);
           if(fix_rec==0) sac_h.gcarc = gcarc;
           sac_h.stlo = gcarc;
           if(fix_rec == 1){
               sac_h.stlo = longs[itrace];
/*               printf( "fix receiver, set stlo=%f\n",sac_h.stlo);*/
               strcpy(t_index,names[itrace]);
           }
           else{
              fop = fopen("tt.txt","wt");
              fprintf(fop,"%d",index);
              fclose(fop);   
              fop = fopen("tt.txt","rt");
              fscanf(fop,"%s",t_index);
              fclose(fop); 
           }
           strcpy(sac_h.kstnm,t_index);
           strcpy(u_sac,t_index);strcat(u_sac,".ce");
           strcpy(v_sac,t_index);strcat(v_sac,".cn");
           strcpy(w_sac,t_index);strcat(w_sac,".cz");
           fprintf(fppm,">> %-8.2f %-8.2f %-8.2f %d\n",gcarc,sac_h.b,sac_h.delta,num);
           for(k=0; k<num; k++)
           {
              offset = k*ntrace + j*xtrace + itrace;
              us[k]= outu[offset]*enlarge;
              vs[k]= outv[offset]*enlarge;
              ws[k]= outw[offset]*enlarge;
              fprintf(fppm,"%-11.4e %-11.4e\n",us[k], ws[k]);
           }

/* call source.f here*/
/* convolve with source parameters */
           if(ite == 1){
               float sdt;
               sdt = ndt;
               xs = 0; // xs and zs are source (x,z)
               xr = gcarc * dis2deg;
//               zr = iz1[itrace+j*xtrace]*h;
               zr = 0;
//               printf("xr=%f zr=%f\n",xr,zr);

//               printf("Convolving with source function\n");
               if(malloced == 0){
                  b = (float* ) malloc(sizeof(float)*num);
                  pp= (float* ) malloc(sizeof(float)*num);
                  malloced = 1;
               }
               for(k=0;k<num;k++) pp[k] =0;
               for(nfn=0;nfn<3;nfn++){
                  if(nfn==0) {for(k=0;k<num;k++) b[k] = us[k];}
                  if(nfn==1) {for(k=0;k<num;k++) b[k] = vs[k];}
                  if(nfn==2) {for(k=0;k<num;k++) b[k] = ws[k];}
                  nnn = num;
                  ifound = 0;
                  for(k=num-1;k>=0;k=k-1){
                     if(fabs(b[k]) <= 1.e-30 && ifound != 1)    nnn = k;
                     else ifound = 1;
                  }
                  if(key == 1) {
                     nfad = 1;
                     sstep(nnn,nfad,sdt,b,pp);
                     nnn = (int) nnn/2;
                  }

                  if(key == 1)  {
                      sdt = ndt *2;
                      sac_h.delta = sdt;
                  }
                  else amom = 1.;

                  stime_(&dt1,&dt2,&dt3,&sdt,pp,&nfa);
                  
                  nnum = nnn;
//                  printf("new point number=%d sdt=%f\n",nnum,sdt); 
 
                  if(key==1 && idev == 1) diff_(&nnum,b,&sdt);

/* modify with 1/sqrt(r) */
                  cc1 = xr - xs;
                  cc2 = zr - zs;
                  roo = sqrt(cc1*cc1 + cc2 *cc2);
                  for(k=0; k< nnum; k++) {
                      cc3 = 1.4142/sqrt(roo);
                      b[k] *= cc3;
                  }
  
                  for(k=0;k<nnum;k++) b[k] = b[k]*amom;
                  convt_(b,&nnum,pp,&nfa,b,&nnum,&sdt);
if(itrace==0 && j==0 && nfn==0) {
 for(k=0;k<nnum;k=k+50)  printf("b[%d]=%e\n",k,b[k]);}
 
//                  printf("output to sac files nfn=%d\n",nfn);
                  if(nfn==0) {for(k=0;k<nnum;k++) us[k] = b[k]*ne_pole;}
                  if(nfn==1) {for(k=0;k<nnum;k++) vs[k] = b[k]*ne_pole;}
/*zhl: 20131120 change sign to - : make upward to be positive*/
/*zhl: 20140215 change to previous version since the sign has been changed in sem2d.c */
                  if(nfn==2) {for(k=0;k<nnum;k++) ws[k] = b[k];}
               }
               sac_h.npts= nnum;
           }

/* output to sac files */
           sac_h.cmpaz = 90.; sac_h.cmpinc = 90.;
           write_sac(u_sac,sac_h,us);
           sac_h.cmpaz = 0. ; sac_h.cmpinc = 90.;
           write_sac(v_sac,sac_h,vs);
           sac_h.cmpaz = 0. ; sac_h.cmpinc =  0.;
           write_sac(w_sac,sac_h,ws);
       }
       index ++;
   }
   fclose(fppm);
}

void sstep(nfa,nfad,dp,rp,pp)
int nfa,nfad;
float dp;
float * rp, *pp;
{
    float convs();
    float f[20000],p[20000];
    int l,j,nff;
    fa_(&dp,&nfa,f);
    nff = (int) (nfa/2);
    l = 0;
    for(j=0; j<nfa; j++) pp[j] = rp[j];
    for(j=0; j<nff;j+=nfad){
       p[l] = convs(dp,nfa,j,pp,f);
       rp[l] = p[l];
       l ++;
    }
}

float convs(del,nf,n,fa,fp)
float del;
int nf,n;
float *fp,*fa;
{
    int i,nn,ndo,ip,np;
    float dn,conv,even,odd,ends;
    nn = n;
    dn = del;
    if(nn < 0) conv = 0.;
    else{
       if(nn>(nf-1)/2) ndo = (nf-1)/2;
       else ndo = nn;
       ip = 1;
       np = 2*nn-1;
       even = fp[ip]*fa[np];
       odd = 0;
       if(ndo > 0){
          for(i=1;i<ndo;i++){
            ip ++;
            np --;
            odd += fp[i]*fa[np];
            ip ++;
            np --;
            even += fp[i]*fa[np];
          }
       }
       ends = fp[0]*fa[2*nn] + fp[ip+1]*fa[np-1];
       conv = dn*(ends+ 4*even + 2*odd)/3;
    }
    return conv;
}


