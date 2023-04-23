*$ debug
c     set rays for FD
CCC*******************************
CCC*******************************
        PARAMETER(ICC=500,idown=1,MODP=5,MODS=3)
        COMMON/RAYS/NA(2*ICC),NRAY(2*ICC)
	CHARACTER*80  IBJ,rayfile,modelfile
	integer isd,ipd,itd,nlay,nlaym,nlmi(100),nlmp(100)
	integer ilayd
	
	OPEN(5,FILE='ray_70.in',STATUS='OLD')
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
	nlay0=nlay/2
	nlayb=nlay0-nlay+1
	do i=1,nlay
	  if(ipd+nlayb.lt.itd) then
	    nlayb=nlayb+1
	  else
	    goto 10
	  endif
	enddo
10	nlay=nlay0-nlayb+1
*<7>
	write(*,*)'input layer number for multiple consideration = ?'
	read(5,'(a)') ibj
	read(5,*) nlaym
	write(*,*) 'nlaym =',nlaym
*<8>
	write(*,*)'input layer indexes and phase num for multiples='
	read(5,'(a)') ibj
	nlayt=1
	do i=1,nlaym
	  read(5,*) nlmi(i),nlmp(i)
	  if(nlmp(i).ge.3) nlmp(i)=4
	  nlayt=nlayt+nlmp(i)
	  write(*,*) i,'  nlmi =',nlmi(i),'  nlmp =',nlmp(i)
	enddo

	close(5)
CCC---- end of parameter input

	OPEN(20,FILE=rayfile,STATUS='UNKNOWN')
	nlayt=nlay*nlayt
	write(20,*)nlayt
	write(*,*)'nlayt =',nlayt,' nlayb =',nlayb,' nlay0 =',nlay0

	do 1000 ip=nlayb,nlay0
	  if(itd.eq.1) itd=2
	  ilayd=ipd+ip
	  iseg1=ilayd-isd
	  iseg2=ilayd-itd+1
	  itotal=iseg1+iseg2+1
	  write(*,*)ip,' ilayd =',ilayd,'  itotal =',itotal
	  do j=1,iseg1
	    na(j)=isd+j
	    nray(j)=MODP
	  enddo
	  do j=1,iseg2
	    na(j+iseg1)=ilayd-j+1
	    nray(j+iseg1)=MODP
	  enddo
	  na(itotal)=1
	  nray(itotal)=MODP
	  write(20,109)itotal,(na(i),i=1,itotal)
	  write(20,110)idown,(nray(i),i=1,itotal)
	  
	  itotal0=itotal
	  do 900 ipm=1,nlaym
*Ps
	  if(nlmp(ipm).gt.0) then
	    if(itd.eq.1) itd=2
	    i=nlmi(ipm)-itd+1
	    do j=itotal0-i,itotal0
	      nray(j)=MODS
	    enddo
	    write(20,109)itotal0,(na(j),j=1,itotal0)
	    write(20,110)idown,(nray(j),j=1,itotal0)
	  endif
*PpPs
	  if(nlmp(ipm).gt.1) then
	    itd=1
	    i=nlmi(ipm)-itd+1
	    do j=itotal0-i+1,itotal0
	      nray(j)=MODP
	    enddo
	    do j=itotal0+1,itotal0+i
	      na(j)=na(2*itotal0-j+1)
	      nray(j)=MODP
	    enddo
	    do j=itotal0+i+1,itotal0+i*2
	      na(j)=na(j-2*i)
	      nray(j)=MODS
	    enddo
	    itotal=itotal0+2*i
	    write(20,109)itotal,(na(j),j=1,itotal)
	    write(20,110)idown,(nray(j),j=1,itotal)
	  endif
*PsPs
	  if(nlmp(ipm).gt.2) then
	    do j=itotal0-i+1,itotal0
	      nray(j)=MODS
	    enddo
	    do j=itotal0+1,itotal0+i
	      na(j)=na(2*itotal0-j+1)
	      nray(j)=MODP
	    enddo
	    do j=itotal0+i+1,itotal0+i*2
	      na(j)=na(j-2*i)
	      nray(j)=MODS
	    enddo
	    itotal=itotal0+2*i
	    write(20,109)itotal,(na(j),j=1,itotal)
	    write(20,110)idown,(nray(j),j=1,itotal)
	  endif
*PPss
	  if(nlmp(ipm).gt.2) then
	    do j=itotal0-i+1,itotal0
	      nray(j)=MODP
	    enddo
	    do j=itotal0+1,itotal0+i
	      na(j)=na(2*itotal0-j+1)
	      nray(j)=MODS
	    enddo
	    do j=itotal0+i+1,itotal0+i*2
	      na(j)=na(j-2*i)
	      nray(j)=MODS
	    enddo
	    itotal=itotal0+2*i
	    write(20,109)itotal,(na(j),j=1,itotal)
	    write(20,110)idown,(nray(j),j=1,itotal)
	  endif
900	continue	  

1000	CONTINUE
109	format(I3,500I4)
110	format(500I2)
	close(20)

999	STOP
	END
