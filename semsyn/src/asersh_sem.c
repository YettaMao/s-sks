/******************************************************************/
/*                                                                */
/* GRT for SH   system                                            */
/* Linked with Spectral element method                            */
/* modified from P-SV system                                      */
/*                                                                */
/*                                                                */
/* Last modified: 20140928 by L Zhao                              */
/* added output green files for plot                              */
/*                Dec. 13, 2007 for SEM                           */
/*                Feb 4, 2008 calculate grid point                */
/* 070209 LZ: output: velocity-> displacement, 201 McCone Hall    */
/*            but there is a question about the displacement !    */
/*            However v-> u is definitely necessary               */
/* 20140928*/
/******************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RAY_LENGTH      1000000 
#define RAY_NUMBER      50000
#define DEGREE_KM       111.2

#define DIRECT          0
#define REFLECTED       1
#define TOTAL           2
#define DOWN           -1
#define UP              1

int     nx;              /* x-dimension of mesh */
int     nz;              /* z-dimension of mesh */
int     nd;              /* z-dimension of mesh */
int     nex;             /* x-dimension of elements */
int     nez;             /* z-dimension of elements */
int     nf         =10;  /* number of fluid meshs */
int     nt;              /* number of time steps */

float   h;               /* spatial grid interval */
float   eh;              /* elemental spacing    */
int     hpn =5;          /* High polynomials number */
int     neh = 1;         /* number of elements for absorbing and attenuation*/
float   dp;              /* time digization interval */
float   xmin;            /* the left edge of the FD region */

int     layerrefl = 0;   /* add the reflection of layers (other than CMB) */
int     convert   = 1;   /* ignore convert phaser or not */
int     layincrust = 2;  /* the layer number for the middle crust (PREM)*/
int     laycrust = 3;    /* the layer number for the moho (PREM) */
int     lay210 = 12;     /* the layer number for the 210  (PREM) */
int     lay390 = 21;     /* the layer number for the 390  (PREM) */
int     lay670 = 37;     /* the layer number for the 370  (PREM) */

/* filenames for rayfile and Green's function */
char rayfile[80];

/* following for the GRT */
float tstart = 0.0;
char raymodel[80], greenfile[80];

/* following for the Gauss sources */
float sgm, ts, *source, *dsource;
float theta, dip, lamda, azmuth;
float azmuth1;

float * b, *pp;

main(ac,av)
int ac; char **av;

{   
    FILE *fop_model, *fop_green, *open_file();
    int n34, istress, nstress, isource; 
    int ispecial;
    int nrec, nfil = 1, zlay, nfil0, nfil1, zlaymax;
    int ix, iz, i, i0, mad1, n, j, k,gl;
    float  xx, thickness, dep;
    int nso, nconv, l2n, nx1, ncent; float fscl;
    float *so0, *so, *dso, *so1;
    float *green;

    FILE * fopu,*fopw;
    int plot_trace=0;
    int itrace,xtrace;  /* trace number in x -direction*/
    float xbeg;  /* first location in x direction(relative to xmin)*/
    int idxplot; /* idxplot*h = x direction interval for plot */
    xtrace = 20; xbeg = 0.0; idxplot = 20; 
    int ixbeg, ix1[200];
    int is_a_rec, npts;
    int kdt = 5; /*output time step*/ 

/* following for the model */
    int jo, nb, jo0, nb0;
    float *cc,  *ss,  *dd, *tth;
    float *cc0, *ss0, *dd0,*tth0;
/* following for the ray descriptions */
    int lfinal, lfinal0; 
    int *nen, *ncoun, *na, *nray; 
    int *nen0,*ncoun0,*na0,*nray0; 
    int bounce = 1;
    int raytype;
    int kdx = 1;
    int mtd = 3;

    int der = 0;
    int iflat = 1;     /* =1 flat 0: no (spherical geometry)*/
    int idirect = 0;   /* index for direct rays (if >0, need to consider converted waves 06/12/04 */

    int grt_(), gauss_(), convt0_(), convt_();
    void getrays();
    void GLL();
    void Grid_Point();/*calculate grid point coordinate */

    void Green_2_Displace();

/* source */
    int ite,key,idev;
    float amom,dt1, dt2, dt3, rhos;

    setpar(ac,av);
    mstpar("shraymodel",  "s", raymodel);
    mstpar("shgreenfile", "s", greenfile);
    getpar("tstart",    "f", &tstart);
    getpar("plot_trace","d", &plot_trace);

    mstpar("eh",        "f", &eh);
    mstpar("dp",        "f", &dp);
    mstpar("nt",        "d", &nt);
    mstpar("nex",       "d", &nex);
    mstpar("hpn",       "d", &hpn);
    getpar("neh",       "d", &neh);
    mstpar("nf",        "d", &nf);
    getpar("kdx",       "d", &kdx);
    getpar("mtd",       "d", &mtd);
    getpar("der",       "d", &der);
    mstpar("nd",        "d", &nd);
    mstpar("xmin",      "f", &xmin);

/* following for output green file for plot 
   By Liang Zhao, Jan 26, 2006
*/
    getpar("xtrace",    "d", &xtrace);
    getpar("xbeg",      "f", &xbeg);
    getpar("idxplot",   "d", &idxplot);
    
    getpar("flat",      "d", &iflat);
    getpar("layerrefl", "d", &layerrefl);
    getpar("convert",   "d", &convert);
    getpar("bounce",    "d", &bounce);

    mstpar("sgm",       "f", &sgm);
    mstpar("ts",        "f", &ts);
    mstpar("theta",     "f", &theta);
    mstpar("dip",       "f", &dip);
    mstpar("lamda",     "f", &lamda);
    mstpar("azmuth",    "f", &azmuth);

    mstpar("Green2Dis", "d", &ite);
    if(ite == 1){
	  mstpar("Moment",        "f", &amom);
	  mstpar("Dt1",           "f", &dt1);
	  mstpar("Dt2",           "f", &dt2);
	  mstpar("Dt3",           "f", &dt3);
	  mstpar("beta",          "f", &rhos);
	  getpar("point_source",  "d", &key);
	  getpar("finaldev",      "d", &idev);
    }
/*set moment = 1 temporarily*/
    amom = 1.0;

/* nfil is the layer number where FD ends: Ling 01/30/04 */
    getpar("nfil",    "d", &nfil);
    endpar();

/* for surface: nfil = 1 */
    if(nfil < 2) { nfil = 1; nf = 2;}

/*  calculate the radiation patterns */
    azmuth1 = azmuth + 180;
    so0 = (float *)malloc(5*sizeof(float));
    fault_(&theta,&dip,&lamda,&azmuth,so0);

    so1 = (float *)malloc(5*sizeof(float));
    fault_(&theta,&dip,&lamda,&azmuth1,so1);
/* FFT the Gauss source function */
    nso = (int) (2.*ts/dp);
    source =(float *) malloc (nso*sizeof(float));
    gauss_(&sgm,&dp,&ts,source,&nso);
/*071809: malloc size of 40200 is consistent with soubroutine convt() in dcsource.f*/
    dso =(float *) malloc (40200*sizeof(float));
    convt0_(&nt,&nso,source,dso,&dp,&nconv,&l2n,&nx1,&ncent,&fscl);
    free(source);

/* Read models from the file model */
    read_model(iflat,raymodel,&jo,&nb,&lfinal,&cc,&ss,&dd,&tth,
	       &nen,&ncoun,&na,&nray);
    printf("jo=%d source at nb=%d\n",jo,nb);
    for(n=0; n<lfinal; n++)	       
        ncoun[n]=((na[nen[0]-2]<=nb) ? UP : DOWN);

/* The layer number where n=3,4 are */
    nrec = na[nen[0]-2];

    fop_green = open_file(greenfile,"wb");

    green =(float *) malloc ((nt+3)*sizeof(float));

/* Calculate the initial responses at m=0 */

/* Evaluate the thickness of the FD region; modify the thickness 
   of the last layer, in order to fit the actual region into the 
   FD grids, nz is the number of total grids in the solid 
*/

    thickness=0.0;
    if(nfil>nrec){
        for(n=nrec; n<nfil; n++)   
	    thickness += tth[n];
    } else {	    
        for(n=nfil-1; n<nrec-1; n++)
        {
	    thickness += tth[n];
	    printf("n = %d, tth= %f, thick = %f\n",n,tth[n], thickness); /*zhl: Feb/18/2006 */
        }
    }	    

    nez = (int) ((thickness-1.e-3)/eh); 
    eh = (thickness-1.e-3)/nez;
    h = eh/hpn;
    printf("xmin=%f nex=%d nez=%d eh=%f\n",xmin,nex,nez,eh);
    nz = nez*hpn + 1;
    nx = nex*hpn + 1;

/* Grid point */
    float *gllp; 
    float *grid_x, *grid_z;/*storing the x and z coordinate for the boundaries*/
    gllp = (float *) malloc((hpn+1)*sizeof(float));
    grid_x  = (float *) malloc(nx*sizeof(float));
    grid_z  = (float *) malloc(nz*sizeof(float));
    GLL(hpn, gllp);
    Grid_Point(gllp,grid_x,grid_z,hpn,eh,nex,nez);
   
    if(plot_trace == 1){
        npts = (int) (nt/kdt);
        ixbeg = (int) (xbeg/h);
        fopu = open_file("grt.v","wt");
        fprintf(fopu,">> %d %d\n", xtrace,npts);
        for(itrace=0; itrace<xtrace;itrace++){
           ix1[itrace] = ixbeg + idxplot*itrace;
        }
    }

    n = nfil-1;
    tth[n]= nz*h- (thickness-tth[n])+0.0001;
    thickness = nz*h ;

    n  = (int) (tth[nrec-1]/h) -3;
    if( nd > n ) nd = n;
    printf("thh=%f, nz=%d nfil=%d nrec=%d nd=%d\n",thickness,nz,nfil,nrec,nd);

    if(nfil < nrec) {
       thickness=0.0;
       for(i=nfil-2; i>=0; i--) thickness += tth[i];
       n = (int) (thickness/h) + 2;
       if(nf > n) {
          printf("nf change from %d to %d\n",nf,n);
	  nf = n;
       }
    }
       
/* This set of points are derived from the initial model and rays */

    cc0   = (float *) malloc ((jo+1)*sizeof(float));
    ss0   = (float *) malloc ((jo+1)*sizeof(float));
    dd0   = (float *) malloc ((jo+1)*sizeof(float));
    tth0  = (float *) malloc ((jo+1)*sizeof(float));

    nen0   = (int *) malloc (RAY_NUMBER * sizeof(int));
    ncoun0 = (int *) malloc (RAY_NUMBER * sizeof(int));
    na0    = (int *) malloc (RAY_LENGTH * sizeof(int));
    nray0  = (int *) malloc (RAY_LENGTH * sizeof(int));

    isource=1;
    fwrite(&nz,sizeof(int),1,fop_green);     
    printf("nz=%d nfil=%d nrec=%d\n",nz,nfil,nrec);

    nfil1 = ((nfil>nrec) ? nfil+1: nfil-1);
    if(nfil1<1) nfil1=1;
    if(nfil1>jo) nfil1=jo+1;
    dep = (nz -1.0)*h;
    zlaymax = modify_model(dep,nrec,nfil1,jo,cc,ss,dd,tth,
                            &jo0,cc0,ss0,dd0,tth0);
    zlaymax -= 1;
    if(zlaymax < 1) zlaymax = 1;     
    printf("zlaymax = %d\n",zlaymax);
    
/* Calculate the initial responses at n=3 and 4 */
    iflat = 1;   /* take care of GRT calculation:model has been flattened*/
    ispecial = 0;
    idirect = 0;   // see the initial of this parameter 06/12/04
    float gcarc,gain,ampmax; /*For output GRT waveform */
    ampmax = -999;
    gain = nx*h/DEGREE_KM/xtrace;

    printf("nx= %d nd= %d nz= %d gain=%f\n",nx,nd,nz,gain);
    int kk,is_nan; /* zhl, 6/16/2006 */

    int malloced =0;
    if(malloced == 0 && ite ==1){
         b = (float* ) malloc(sizeof(float)*nt);
         pp= (float* ) malloc(sizeof(float)*nt);
         malloced = 1;
    }
    float xs,zs,xr,zr;
    xs=0; // x of source 
    zs=0; // z of source, temporarily
    for(n34=3; n34<=3; n34++){ 
	dep = neh * h; /*zhl: 20131216 adjust the position of GRT-SEM interface */
        zlay = modify_model(dep,nrec,nfil1,jo,cc,ss,dd,tth,
                            &jo0,cc0,ss0,dd0,tth0);
                                   
	nfil0 = nfinal(jo0,ss0,nrec);
	if(nfil<nrec) nfil0 = 1;     
        if(layerrefl == 1){nfil0 = zlaymax;}
	printf("nfil0=%d zlay=%d\n",nfil0, zlay);
    
        for(istress=0; istress<1; istress++){
/*          nstress =0        v   velocity 
            nstress =12       T12 stress 
*/
            if(n34==3 && istress==0) nstress=0;

	    printf("Getting the Rays for DIRECT for right REGION nstress=%d\n",nstress);
            getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb, 
         	    &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,DIRECT,idirect);
	    printf("End Getting the Rays for DIRECT for right REGION\n");

            for(i=0; i< nx; i++){
                is_a_rec = 0;
                for(itrace =0; itrace< xtrace;itrace++){
                         if(ix1[itrace] == i) is_a_rec =1;
                }
                xx=xmin + grid_x[i];
                gcarc = xx/DEGREE_KM;
		if(!(i % 200)) printf("DIRECT i=%d xx=%f gcarc=%f\n",i,xx,gcarc);
                so = (xx > 0) ? so0 : so1; 
                grt_(&isource,&iflat,&ispecial,&nb0,&jo0,cc0,ss0,dd0,tth0,
		     &lfinal0,nen0,na0,nray0,ncoun0,
                     &dp,&nt,&xx,&tstart,&nstress,so,&der,green);

/*  P-SV subroutine has mtd              grt_(&isource,&iflat,&ispecial,&mtd,&nb0,&jo0,cc0,ss0,dd0,tth0,
		     &lfinal0,nen0,na0,nray0,ncoun0,
                     &dp,&nt,&xx,&tstart,&nstress,so,&der,green); */
                convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,&ncent,&fscl);  

/* zhl for debug, find nan output, 6/16/2006 */
                is_nan = -1;  
                for(kk=0; kk < nt; kk ++){
                    if(!(green[kk]>-2 && green[kk]< 2))  {
                        is_nan = kk ;
                        if(kk >1 && kk<nt && green[kk+1]>-2 && green[kk+1]< 2 ) {
                            green[kk] = (green[kk-1] + green[kk+1])/2.0;
                            is_nan = -1;
                        }
                    }
                    if(istress==0 && i ==0) {
                          if(fabs(green[kk])>ampmax) ampmax = fabs(green[kk]);
                    }
                } 

/*LZ: change the output to be displacement instead of velocity, 070209, 201 McCone Hall */
/*LZ: Does not work, output is a step-like function */
//                Velocity_2_Displace(green,nt,dp);
/*LZ: output green function to displacement, 20110527, in LPG, Nantes*/
/*result: almost similar to that done after SEM calculation */
                if(ite == 1) {   
                     xr= xx; 
                     zr= nz*h - dep;
                     Green_2_Displace(green,nt,dp,xs,zs,xr,zr,amom,key,dt1,dt2,dt3,idev);  
                }

                fwrite(green,sizeof(float),nt,fop_green);  

                if(is_a_rec ==1 && plot_trace ==1){
                    if(n34 == 3){
                       fprintf(fopu,">> %f\n", xx);
                       for(kk=0; kk < nt; kk ++){
                          if(!(kk%kdt)) fprintf(fopu,"%-9.3f %-12.6e\n",kk*dp,gcarc+green[kk]/ampmax*gain);
//                          if(!(kk%500))  printf("%-9.3f %-12.6e ampmax=%e\n",kk*dp,green[kk],ampmax);
                       }            
                    }          
                }
                if(is_nan > 0) printf("GRT Warning: Bottom ix=%-3d nstress=%-3d green=%-10.e\n",i,nstress, green[0]); 
            }
        }
    }

    if(plot_trace == 1 ){
        fclose(fopu); ;
    }

    idirect = 0;
    for(n34=3; n34<=3; n34++){
        xx = xmin + neh*h; /*zhl: 20131216 adjusting position of GRT-SEM interface */
        so = (xx > 0) ? so0 : so1; 

	for(istress=0; istress<1; istress++){
        
/*          nstress =0        v   velocity 
            nstress =12       T12 stress 
*/
            if(n34==3 && istress==0) nstress=0;

	    ispecial =  0;
            for(iz=0; iz<nz; iz++){
                dep = grid_z[iz];
                if(dep == 0.) dep = 1e-4;
                if(!(iz % 100)||iz==nz-1) printf("Total iz=%d  dep=%f\n", iz,dep);
		if(iz < 0 || (iz == 0 && istress == 0))
		    raytype = REFLECTED;
                else
		    raytype = TOTAL;

                if(xx < 0 && raytype == REFLECTED) 
                    raytype = DIRECT;

                if(xx > 0 || raytype == DIRECT){
                    zlay = modify_model(dep,nrec,nfil1,jo,cc,ss,dd,tth,
                                       &jo0,cc0,ss0,dd0,tth0);


                    nfil0 = nfinal(jo0,ss0,nrec);
  		    if(nfil<nrec) nfil0 = 1;  
                    if(layerrefl == 1){nfil0 = zlaymax;}

                    getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb, 
                 	    &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,raytype,idirect);

                    grt_(&isource,&iflat,&ispecial,&nb0,&jo0,cc0,ss0,dd0,tth0,
		         &lfinal0,nen0,na0,nray0,ncoun0,
                         &dp,&nt,&xx,&tstart,&nstress,so,&der,green);

/*                    grt_(&isource,&iflat,&ispecial,&mtd,&nb0,&jo0,cc0,ss0,dd0,tth0,
		            &lfinal0,nen0,na0,nray0,ncoun0,
                            &dp,&nt,&xx,&tstart,&nstress,so,&der,green); */
                    convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,
                       &ncent,&fscl);  

                    if(raytype == DIRECT){
                        for(i = 0; i<nt; i++)
                            green[i] *= -1.0;
                    }
                } else {
                    for(i = 0; i<nt; i++)
                        green[i] = 0.0;
                }

/* zhl for debug, find nan output 6/16/2006 */
                is_nan = -1;  
                for(kk=0; kk < nt; kk ++){
                     if(!(green[kk]>-2 && green[kk]< 2))  {
                        is_nan = kk ;
                        if(kk >1 && kk<nt && green[kk+1]>-2 && green[kk+1]< 2 ) {
                            green[kk] = (green[kk-1] + green[kk+1])/2.0;
                            is_nan = -1;
                        }
                     }
                } 
                if(ite == 1) {   
                     xr= xx; 
                     zr= nz*h - dep;
                     Green_2_Displace(green,nt,dp,xs,zs,xr,zr,amom,key,dt1,dt2,dt3,idev);  
                }

                fwrite(green,sizeof(float),nt,fop_green);     

                if(is_nan > 0) printf("GRT Warning: Left iz=%-3d nstress=%-3d green=%-10.e\n",iz,nstress, green[0]); 
            }
        }
    }
    
    fclose(fop_green);
}    

/* getrays determins the ray parameters for the responses 
   at zlay layer. The input is the ray group for the n=3,4, 
       which are: lfinal, nen, na, nray 
                  zlay: the layer where the receive is
                  nrec: the layer where n=3,4 lines are
                  nfil: the layer where solid-liquid interface is
                  type =0: direct wavefield
                       =1: reflected wavefield
                       =2: total wavefield
*/                  

void getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil,jo,nb,
        lfinal0,nen0,ncoun0,na0,nray0,jo0,nb0,bounce,raytype,idirect)
int lfinal, *nen, *na, *nray, *ncoun; 
int zlay, nrec, nfil, bounce, raytype, jo, nb, jo0, *nb0;
int *lfinal0, *nen0, *na0, *nray0, *ncoun0;
int idirect;

{
    int ray, k, n, j, rays1, rays2, nenn, rays, layer2, incre=0;
    int P_SV, SV_P;

    int **rtype, *type;
    int nbounce, i, m, num;
    int nref, nref1, nref2;
    int myray;     /* Ling: 02/11/04 */
    int layspecial = nfil;  // Ling: 06/12/04

    *nb0=nb;
    if(nfil<nrec+1 && jo != jo0){
        *nb0=nb+1; nrec++; incre++;
    }

    n=k=0;
    myray=0;  /* Ling: 02/11/04 */
    
    for (ray=0; ray<lfinal; ray++){

//        printf("ray = %d,  myray = %d\n",ray+1,myray);        
//        P_SV = nray[nen[ray]-2];
        P_SV = nray[myray+nen[ray]-2];      /* Ling: 02/11/04 */
	SV_P = ( (P_SV==3) ? 5 : 3);
//	printf("P_SV = %d,  SV_P = %d\n",P_SV,SV_P);

	switch(raytype){

	case DIRECT:
	    rays1=0; rays2=0; break;
	case REFLECTED:     
	    rays1=1; rays2=1; break;
	case TOTAL:     
	    rays1=0; rays2=1; break;

        }

	rtype = imatrix(0,bounce,1,power(2,bounce+2));
	type  = ivector(0,bounce);

        for(rays=rays1; rays<=rays2; rays++){

            if(rays == 0){
//                printf("DIRECT WAVE\n");
                nenn=0;

                /* keep the ray paths to n=3, 4 */
                for(j=0; j<nen[ray]-1; j++){
                    na0[k]   =na[j+myray]+incre;    
                    nray0[k] =nray[j+myray];    
                    k++;
                    nenn++;
                } 

                /* calculate the direct wave from the source */
		nbounce  = 1;
		type[0] = P_SV;
		nen0[n] = nenn;

		k  = MultiBounces(nrec,nfil,zlay,k,nbounce,type,
				   na0,nray0,&ncoun0[n],&nen0[n],layspecial);
                n++;
            }     
            else { 
            
                /* reflected waves from the solid-liquid boundary */
		layspecial = nfil;
                for(nbounce=2; nbounce<=bounce; nbounce++){

		    if((zlay == nfil) && ((nbounce % 2) == 0)) 
			continue;     /* skip the uparriving ray */
		    if((zlay == nrec) && ((nbounce % 2) != 0)) 
			continue;     /* skip the downarriving ray */
		    
                    num = BounceRays(nbounce,rtype,P_SV,SV_P);

		    for(j=1; j<=num; j++){

                      /* modified for only valid when nfil > nrec */
                       /* not tested for the other cases */
                       nref1 = nfil;
                       nref2 = (layerrefl == 1) ? zlay +1 : nfil + 1;
		       if(nfil < nrec){nref2 = (layerrefl == 1) ? zlay - 1 : nfil - 1;}
                       for(nref = nref1; nref > nref2; nref--){

                        nenn=0;

                        /* keep the ray paths to n=3, 4 */
                       for(m=0; m<nen[ray]-1; m++){
                           na0[k]   =na[m+myray]+incre;    
                           nray0[k] =nray[m+myray];    
                           k++;
                           nenn++;
                        } 

			for(i=0; i<nbounce; i++){
			    type[i]  = rtype[i][j];
//			    printf("%d ",type[i]);     
                         }
//			 printf("\n"); 
		    
		        nen0[n] = nenn;
		        k  = MultiBounces(nrec,nref,zlay,k,nbounce,type,
				   na0,nray0,&ncoun0[n],&nen0[n],layspecial);


                        n++;
                        }
                    }
                }
            }            
        }
	free_imatrix(rtype,0,bounce,1,power(2,bounce+2));
	free_ivector(type,0,bounce);
        myray += nen[ray];                  /* Ling: 02/11/04 */
    }

    *lfinal0=n;
}            

/*
    int BounceRays(int nbounce, int ** rtype, P_SV SV_P

    rtype[i][j]  i -> ray segments j -> number of rays

*/

int BounceRays(nbounce,rtype,P_SV,SV_P)
int nbounce, P_SV, SV_P;
int **rtype;
{
    int i, j, m, num;

    /* first fundamental ray */ 
    for(i=0; i<nbounce; i++)       
         rtype[i][1] = P_SV;

    num = 1;
 
    return num;
}

int MultiBounces(nrec,nfil,zlay,k,bounce,type,na,nray,ncoun,nen,layspecial)
int nrec, nfil, zlay; int k, bounce, *type;
int *na, *nray, *ncoun, *nen;
int layspecial;
{
    int j, i; int start, end1, end2, temp; int going;
    int nen0, k0;

    start = ((nfil > nrec) ? nrec + 1 : nrec-1);
    end1  = ((nfil > nrec) ? nrec + 1 : nrec-1);
    end2  = ((nfil > nrec) ? nfil - 1 : nfil + 1);
    going = ((nfil > nrec) ? DOWN : UP);
    
    if((zlay == nrec) && (nfil > nrec))
	end1 = nrec + 2;


    nen0 = *nen;  k0 = k;
    i = 0;

    while ( i < bounce){

	if(i == bounce-1)
	    /* last segement of the ray */
	    end2  = ((going == UP) ? zlay+1 : zlay);

        if(going == DOWN){
	    for(j=start; j<=end2; j++){
	        na[k0]   = j;
	        nray[k0] = type[i];
	        k0++; nen0++;
            }
	    going = UP;
        } else { 
	    for(j=start; j>=end2; j--){
	        na[k0]   = j;
	        nray[k0] = type[i];
	        k0++; nen0++;
            }
	    going = DOWN;
        }
	if(i==0){
	    start = end2; end2 = end1; 
        } else {
	    temp = start; start = end2; end2 = temp; 
	}

	i++;
    }
    na[k0]   = 1;
    nray[k0] = 4; 
    k0++; nen0++;
    *ncoun  = -going;
    *nen    =  nen0;
    return k0;
}


/* modify the model parameters for each point at m=0*/
/* returns the layer number of the receiver */ 

int modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,
                 jo0,cc0,ss0,dd0,tth0)
float dep; int nrec, nfil; 
int jo; float *cc, *ss, *dd, *tth;
int *jo0; float *cc0, *ss0, *dd0, *tth0;

{   
    int layer, zlay;
    float thickness=0.0;

    if(nrec<nfil){
        for(layer=nrec; layer<nfil+500; layer++){
	    thickness += tth[layer];
	    if(thickness >= dep) break;
        }
        zlay = ( (dep<0.00001) ? nrec : layer+1);
    } else {
        for(layer=nrec-2; layer>nfil-1; layer--){
	    thickness += tth[layer];
	    if(thickness >= dep) break;
        }
        zlay = ( (dep<0.00001) ? nrec : layer+1);
    }

    if(zlay == nrec) thickness = 0.0;

    if(dep > thickness)zlay = nfil;

    if(fabs(thickness-dep) <0.0001 || zlay==nfil){
        for(layer=0; layer<jo; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }            
	if(zlay==nfil)tth0[nfil-1]  = dep-thickness;
        *jo0 =jo; 
    }    
    else {
        for(layer=0; layer<zlay; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }
        
    
         cc0[zlay]  =  cc[zlay-1]+0.00001; 
         ss0[zlay]  =  ss[zlay-1]+0.00001; 
         dd0[zlay]  =  dd[zlay-1]+0.00001; 
         tth0[zlay]   = (nrec<nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep); 
         tth0[zlay-1] = (nrec>nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep); 

         for(layer=zlay+1; layer<jo+1; layer++){
            cc0[layer]  =  cc[layer-1];
            ss0[layer]  =  ss[layer-1];
            dd0[layer]  =  dd[layer-1];
            tth0[layer] = tth[layer-1];
        }

        *jo0 =jo+1;
    }     

    return zlay;
}

int power(int a, int n)
{
    int j, num;

    num = a;
    for(j=2; j<=n; j++)
	num *= a;

    return num;
}

/*subroutine calculating grid point coordinates
Only for the left and right boundaries
Requires:
   gllp: GLL points
   hpn: grid point number per element
   nex, nez: number of horizontal and vertical elements
   eh: element spacing
   bxp: bottom x point
   lzp: left z point 
*/
void  Grid_Point(gllp,bxp,lzp,hpn,eh,nex,nez)
float * gllp;
float * bxp, *lzp;
int hpn;
float eh;
int nex, nez;
{
    int ex,ez,m,index;
    float rx,rz;
    /*bottom point*/
    bxp[0] = 0.0;
    for(ex=0;ex<nex;ex++){
       for(m=1;m<= hpn;m++){
          index=ex*hpn + m;
          rx = gllp[m];
          bxp[index] = ex*eh + (rx+1)/2.*eh;
          if(ex == 10) printf("index=%d rx=%f bxp=%-8.3f\n", index, rx,bxp[index]);
       }
    }
    /*left point*/
    lzp[0] = 0.0;
    for(ez=0;ez<nez;ez++){
       for(m=1;m<= hpn;m++){
          index=ez*hpn + m;
          rz = gllp[m];
          lzp[index] = ez*eh + (rz+1)/2.*eh;
          if(ez == 10) printf("index=%d rz=%f lzp=%-8.3f\n", index, rz,lzp[index]);
       }
    }
}

/* subroutine to change the output to be displacement instead of velocity
nt: time step number
value storing time sequence of velocity at the beginning
-> storing time sequence of displacement 
dt: time step
*/
void Velocity_2_Displace(value,nt,dt)
float * value;
int nt;
{
    int i;
    float tempf = 0.0;
    for(i=0;i<nt;i++){
        tempf += value[i]*dt;
        value[i] = tempf;
    }
}

/* subroutine to change the output to be displacement instead of green function
num: time step number
value storing time sequence of velocity at the beginning
-> storing time sequence of displacement 
ndt: time step
xr,zr: receiver (x,z)
xs,zs: source (x,z)
key: do line -> point source correction when =1
amom: moment of source
dt1,dt2,dt3: source parameters
*/
void Green_2_Displace(value,num,ndt,xs,zs,xr,zr,amom,key,dt1,dt2,dt3,idev)
float * value,ndt;
float xr,zr;
float xs,zs;
float amom,dt1,dt2,dt3;
int key;
int num;
{
     int i,k,tempn;
     float sdt;
     sdt = ndt;
     int nfn, nfa, ifound, nnn,nfad,nnum;
     float cc1, cc2,cc, roo;
     void sstep();
     int convt3_();

     for(k=0;k<num;k++) pp[k] =0;
     for(k=0;k<num;k++) b[k] = value[k];

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

     if(key == 1)  sdt = ndt *2;
     else amom = 1.;

     stime_(&dt1,&dt2,&dt3,&sdt,pp,&nfa);
                  
     nnum = nnn;
     if(key==1 && idev == 1) diff_(&nnum,b,&sdt);

/* modify with 1/sqrt(r) */
     cc1 = xr - xs;
     cc2 = zr - zs;
     roo = sqrt(cc1*cc1 + cc2 *cc2);
     for(k=0; k< nnum; k++) {
          cc = 1.4142/sqrt(roo);
          b[k] *= cc;
     }
  
     for(k=0;k<nnum;k++) b[k] = b[k]*amom;

     convt3_(b,&nnum,pp,&nfa,b,&nnum,&sdt);

     for(k=0;k<num;k++) { 
          tempn= (int) (k*ndt/sdt);
          value[k] = b[tempn];
     }
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

