ccc ******************************
c   Origined by Ling Chen, IGG, CAS
c     set rays for FD
c   Modified by Liang Zhao, IGG, CAS
c     for SH waves
c     SV=3  SH=4 P=5
CCC*******************************
CCC*******************************
        PARAMETER(ICC=500,idown=1,MODP=5,MODSV=3,MODSH=4)
        COMMON/RAYS/NA(2*ICC),NRAY(2*ICC)
	CHARACTER*90  IBJ,rayfile,modelfile
	integer isd,ipd,itd,nlay,nlaySK
	integer ilayd
	
	OPEN(5,FILE='ray_sh_zl.in',STATUS='OLD')
**<1>
*	write(*,*)'input model file name: modelfile'
*	read(5,'(a)') ibj
*	READ(5,'(a)') modelfile
*	WRITE(*,'(a)') modelfile
*<2>
	write(*,*)'output ray file name: rayfile'
	read(5,'(a)') ibj
	READ(5,'(a)') rayfile
	WRITE(*,'(a)') rayfile
*<3>
	write(*,*)'input source layer isd = ? '
	read(5,'(a)') ibj
	READ(5,*) isd
	WRITE(*,*) 'isd =',isd
*<4>
	write(*,*)'input deepest layer ipd = ?'
	read(5,'(a)') ibj
	read(5,*) ipd
	write(*,*) 'ipd =',ipd
*<5>
	write(*,*)'input last layer itd = ?'
	read(5,'(a)') ibj
	read(5,*) itd
	write(*,*) 'itd =',itd
*<6>
	write(*,*)'input number of rays to be considered nlay = ?'
	read(5,'(a)') ibj
	read(5,*) nlay
	write(*,*) 'nlay =',nlay
*<7>
	write(*,*)'input layer number for S->K conversion in PREM'
	read(5,'(a)') ibj
	read(5,*) nlaySK
	write(*,*) 'nlaySK =',nlaySK

	nlay0=nlay/2
	nlayb=ipd - nlay0
        nlaye = ipd + nlay0 - 1
        if(nlaye.ge.nlaySK) then
           write(*,*) 'nlay= ',nlay,'nlayb= ',nlayb,'nlaySK=', nlaySK
           nlaye = nlaySK-1
        endif
        write(*,*) 'nlay= ',nlay,'nlayb= ',nlayb,'nlaye=',nlaye
	close(5)
CCC---- end of parameter input

	OPEN(20,FILE=rayfile,STATUS='UNKNOWN')
	nlayt=nlaye-nlayb + 1
c	write(20,*)nlayt
	write(20,*)nlayt
	write(*,*)'nlayt=',nlayt,'nlayb= ',nlayb,'nlaye=',nlaye

	do 1000 ip=nlayb,nlaye
	  ilayd=ip
	  iseg1=ilayd-isd
	  iseg2=ilayd-itd+1
	  itotal=iseg1+iseg2+1
	  write(*,*)ip,' ilayd =',ilayd,'  itotal =',itotal
	  do j=1,iseg1
	    na(j)=isd+j
            nray(j)=MODSH
	  enddo
	  do j=1,iseg2
	    na(j+iseg1)=ilayd-j+1
            nray(j+iseg1)=MODSH
	  enddo
	  na(itotal)=1
	  nray(itotal)=MODSH
	  write(20,109)itotal,(na(i),i=1,itotal)
	  write(20,110)idown,(nray(i),i=1,itotal)
1000	CONTINUE
109	format(I3,500I4)
110	format(500I2)
	close(20)

999	STOP
	END
