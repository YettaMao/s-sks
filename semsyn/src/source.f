
C********************************************************************C
C*                                                                  *C
C*                source.f                                          *C 
C*                                                                  *C
C*       Convert Green's function from hybrid method to SAC files   *C
C*                                                                  *C
C*   Written by: Lianxing Wen                                       *C
C*               Seismological Lab. Caltech                         *C
C*                                                                  *C
C*   Last modified: Oct. 22, 1996                                   *C
C*   Modified by Liang Zhao for anisotropic media                   *C
C*            Oct 17, 2006                                          *C
C*                                                                  *C
C********************************************************************C

      subroutine colb2(nn,datai,signi)
      dimension datai(*)
      n=2**(nn+1)
      j=1
      do 5 i=1,n,2
      if(i-j)1,2,2
    1 tempr=datai(j)
      tempi=datai(j+1)
      datai(j)=datai(i)
      datai(j+1)=datai(i+1)
      datai(i)=tempr
      datai(i+1)=tempi
    2 m=n/2
    3 if(j-m)5,5,4
    4 j=j-m
      m=m/2
      if(m-2)5,3,3
    5 j=j+m
      mmax=2
    6 if(mmax-n)7,10,10
    7 istep=2*mmax
      theta=signi*6.28318531/float(mmax)
      sinth=sin(theta/2.)
      wstpr=-2.0  *sinth*sinth
      wstpi= sin(theta)
      wr=1.
      wi=0.
      do 9 m=1,mmax,2
      do 8 i=m,n,istep
      j=i+mmax
      tempr=wr*datai(j)-wi*datai(j+1)
      tempi=wr*datai(j+1)+wi*datai(j)
      datai(j)=datai(i)-tempr
      datai(j+1)=datai(i+1)-tempi
      datai(i)=datai(i)+tempr
    8 datai(i+1)=datai(i+1)+tempi
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
    9 wi=wi*wstpr+tempr*wstpi+wi
      mmax=istep
      go to 6
   10 return
      end

      subroutine convt3(x,nx,y,ny,z,nz,dt)
      dimension x(*),y(*),z(*),c(80200),d(80200)
      nz=max0(nx,ny)
      call log2fd(nz,n,l2n)
      l2n=l2n+1
      n=2*n
      nx1=nx+1
      ny1=ny+1
      do 1 i=1,nx
      j2=2*i
      j1=j2-1
      c(j1)=x(i)
1     c(j2)=0.
      do 2 i=nx1,n
      j2=2*i
      j1=j2-1
      c(j1)=0.
2     c(j2)=0.
      do 3 i=1,ny
      j2=2*i
      j1=j2-1
      d(j1)=y(i)
3     d(j2)=0.
      do 4 i=ny1,n
      j2=2*i
      j1=j2-1
      d(j1)=0.
4     d(j2)=0.
      call colb(l2n,c,-1.)
      call colb(l2n,d,-1.)
      nhalf=n/2
      ncent=nhalf+1
      do 6 i=1,ncent
       j2=2*i
      j1=j2-1
      e1=c(j1)*d(j1)-c(j2)*d(j2)
      e2=c(j1)*d(j2)+c(j2)*d(j1)
      c(j1)=e1
6     c(j2)=e2
      call conj(c,ncent)
      call colb(l2n,c,1.)
      fscl=dt/float(n)
      do 7 i=1,nz
      j1=2*i-1
7     z(i)=c(j1)*fscl
      return
      end

      subroutine conj2(c,ncent)
      dimension c(*)
      nhm1=ncent-2
      ncent2=ncent*2
      do 1 i=1,nhm1
      l2=2*i
      l1=l2-1
      l3=l2+1
      k1=ncent2+l1
      k2=ncent2+l2
      j1=ncent2-l3
      j2=ncent2-l2
      c(k1)=c(j1)
1     c(k2)=-c(j2)
      return
      end
      subroutine log2fd2(np,n,l2n)
      n1=0
      n=1
      l2n=0
1     continue
      if((np.gt.n1).and.(np.le.n)) go to 2
      n1=n1*2
      n=n*2
      l2n=l2n+1
      go to 1
2     return
      end

	subroutine stime(dt1,dt2,dt3,dt,f,nf)
	dimension f(20000)

c            write(*,*) dt1,dt2,dt3,dt
	    n1  = int(dt1/dt+0.1)
	    n2  = int(dt2/dt+0.1)
	    n3  = int(dt3/dt+0.1)
	    n12 = n1+n2+1
	    n   = n1+n2+n3+1

	    nf =n

	    sum=0.
	    f(1)=0.0
	    do i=2,n
	         k=i
	        f(k)=1.
	        if(i.lt.(n1+1))f(k)=(i-1.0)*(1./n1)
	        if(i.gt.n12)f(k)=1.0-(i-n12)*(1./n3)
 	        sum = sum+(f(k)+f(k-1))*0.5
            enddo

	    sum=sum*dt
	    do i=1,nf
 	        f(i)=f(i)/sum
            enddo
	return
	end

      subroutine fa(del,j2,f)
	dimension f(20000)
      do 5 j=2,j2   
      f(j)  = 1./sqrt((j-1)*del)    
 5      continue  
      f(1)  = (11.-4.*2.**.5)/sqrt(2.*del)  
      return
      end

      subroutine diff(n,pp,dt)
      real pp(*)
      y1=pp(1)
      pp(1)=0.0
      do jk=2,n
          y2=pp(jk)
          pp(jk)=(pp(jk)-y1)/dt
          y1=y2
      enddo
      return
      end    
