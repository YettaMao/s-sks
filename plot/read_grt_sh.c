/*********************************************************************/
/*                                                                   */
/*          C code for read the results from asersh_sem              */
/*                                                                   */
/*  Written by: Liang Zhao, IGG, CAS                                 */
/*  Last modified: 20141024                                          */
/*                                                                   */
/*********************************************************************/

/* modified history:
*/

/*
   Usage:  ./read_grt par=*.par
   Input:  get par from green.out

   Output: GRT_Plot_sh.xy 
*/

#include   <stdio.h>
#include   <stdlib.h>
#include   <sys/file.h>
#include   <fcntl.h>
#include   <math.h>
#include   <malloc.h>

float *vt34;
int nx,nz,nex,nez;
int   nfil = 1;  /* the layer number where SEM ends */

main(ac,av)
int ac; 
char **av;
{
   int i,k,j,num,size_f;
   int flat = 1;
   int ix1[700], itrace;
   int iz1[700];
   int ixbeg,izbeg,idxplot,idzplot=15;;
   int xtrace = 700, ztrace=20;
   int n,nt,len,lenz,hpn;
   float dp,dt;
   int  kdx =1;
   float step = 0.2; /* time step for output */
   int istep;

   char raymodel[80];   /* file for model in GRT calculation */
   char greenfile[80];
   int mtotal;
   int  **matrix(); 
   void interpl();

   FILE *infile ; 
   FILE *fop; 
   float xmin,xbeg, zbeg,tstart,gcarc,depth ;
   float h,eh;
   float time;
 
   setpar(ac,av);
   mstpar("shgreenfile",   "s", greenfile);
   mstpar("raymodel",      "s", raymodel);
   getpar("nfil",          "d", &nfil);
   mstpar("eh",            "f", &eh);
   mstpar("nex",           "d", &nex);
   mstpar("hpn",           "d", &hpn);
   mstpar("nt",            "d", &nt);
   mstpar("dp",            "f", &dp);
   getpar("xmin",          "f", &xmin);
   getpar("tstart",        "f", &tstart);
   getpar("xtrace",        "d", &xtrace);
   getpar("ztrace",        "d", &ztrace);
   getpar("kdx",           "d", &kdx);
   mstpar("xbeg",          "f", &xbeg);
   mstpar("zbeg",          "f", &zbeg);
   mstpar("idxplot",       "d", &idxplot);
   mstpar("idzplot",       "d", &idzplot);
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
   printf(" %f h= %f\n",eh,eh/hpn);

   /*zhl: for reading the record */   
   nx = nex*hpn+1;
   nz = nez*hpn+1;
   ixbeg = (int) (xbeg/h);
   printf("nt= %d ixbeg= %d \n", nt, ixbeg);
   for(itrace=0; itrace<xtrace; itrace++)
   {
      ix1[itrace] = ixbeg + idxplot * itrace;
   }

   izbeg = (int) (zbeg/h);
   for(itrace=0; itrace<ztrace; itrace++)
   {
      iz1[itrace] = nz + izbeg - idzplot * itrace;
printf("itrace=%d iz1=%d\n",itrace,iz1[itrace]);
   }
   
   fop = fopen("dimension.txt","rt");
   fscanf(fop,"%d",&nz);
   fclose(fop);   

   mtotal = 1*(nx + nz);
   vt34 = (float *) malloc(mtotal*sizeof(float));
 
   istep = (int) (step/dp);
   if(istep <1) istep = 1;
   infile = fopen(greenfile,"rb");     

   int it,ntrace,offset;
   ntrace = xtrace;
   num = 0;
   printf("\n Now, read the results ...\n");
   printf("xtrace=%d istep=%d\n",xtrace,istep);
   printf("nx=%d nz=%d \n",nx, nz);
   float *v3, *holdv;
   float *zv;
   int ix,iz;
   len = (int) (nt*ntrace/istep);
   lenz = (int) (nt*ztrace/istep);

   size_f = sizeof(float);
   holdv = (float *)malloc(len*sizeof(float));

   zv = (float *)malloc(lenz*sizeof(float));

   num = -1;
   for(it=0; it<nt; it++)
   {
       fread(vt34,size_f,mtotal,infile);
       if(!(it%istep)) {      
            num ++ ;
            /*bottom boundary*/
            v3 = vt34 ;
            for(ix=0; ix<nx; ix++) {
                v3++;
                for(itrace=0; itrace < xtrace; itrace++) {
                    if(ix1[itrace] == ix) {
                       offset = num*xtrace + itrace;
                       holdv[offset] = v3[0];  
                       if(itrace==0 && it==0) {
                            printf("offset=%d v= %-11.3e\n",offset, holdv[offset]);
                       }
                    }
                }
            }
            /*left boundary*/
            v3 = vt34 + nx;
            for(iz=0; iz<nz; iz++) {
                v3++;
                for(itrace=0; itrace < ztrace; itrace++) {
                    if(iz1[itrace] == iz) {
                       offset = num*ztrace + itrace;
                       zv[offset] = v3[0];  
                       if(itrace==1 && it==0) {
                            printf("offset=%d v= %-11.3e\n",offset, zv[offset]);
                       }
                    }
                }
            }
       }
   }
   fclose(infile);

   fop = fopen("GRT_Plot_sh.xy","wt");     
  
   fprintf(fop,"%-5d %-5d\n",num,xtrace);
   for(itrace=0; itrace< xtrace; itrace ++)
   {
       gcarc = xmin + xbeg + idxplot*itrace*h;
       if(itrace==0) printf("No.1 receiver at %-8.3f\n",gcarc/111.2);
       fprintf(fop,">> 1 %-8.3f\n", gcarc);
       for(k=0; k<num; k++)
       {
            time = tstart + k*istep*dp;
            offset = k*xtrace + itrace;
            fprintf(fop,"%-11.5f %-11.5e\n", time, holdv[offset]);
       }
   }
   fclose(fop);

   fop = fopen("GRT_Plot_z.xy","wt");     
  
   fprintf(fop,"%-5d %-5d\n",num,ztrace);
   for(itrace=0; itrace< ztrace; itrace ++)
   {
       depth = -zbeg + idzplot*itrace*h;
       if(itrace==0) printf("No.1 receiver at %-8.3f\n",depth);

       fprintf(fop,">> 1 %-8.3f\n", depth);
       for(k=0; k<num; k++)
       {
            time = tstart + k*istep*dp;
            offset = k*ztrace + itrace;
            fprintf(fop,"%-11.5f %-11.5e\n", time, zv[offset]);
       }
   }
   fclose(fop);
}


