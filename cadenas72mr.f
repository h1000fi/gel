      subroutine cadenas72mr(chains,ncha, rpos, zpos, long1)

      use mchains
      use mparameters
      use mncells
      use msegme
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mporesystem

      implicit none

      integer ncha
      integer chains(ncha_max,maxlong)
      integer long1
      
      real*8 y(maxlong),z(maxlong)
      
      real*8 rvect, zvect
      
      integer i,state,ii,j,ive,jve
      real*8 rn,state1,sitheta,cotheta,dista

      real*8 siphip,cophip
      character*1 test
      real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
      
      real*8 x(3),xend(3,maxlong),xendr(3,maxlong), xendr2(3,maxlong)
      real*8 theta
      real*8 zpos, rpos
      
      real*8 tempp(3)

      logical   outside
      external  outside

      integer iR, iZ

      real*8 xo(3)
      
      sitheta=sin(68.0*pi/180.0)
      cotheta=cos(68.0*pi/180.0)
      siphip=sin(120.0*pi/180.0)
      cophip=cos(120.0*pi/180.0)

      theta = pi
            
 223  x(1)=lseg
      x(2)=0.0
      x(3)=0.0
      
      xend(1,1)=lseg
      xend(2,1)=0.0
      xend(3,1)=0.0
      
      tt(1,1)=cotheta
      tt(1,2)=sitheta
      tt(1,3)=0.0
      tt(2,1)=sitheta
      tt(2,2)=-cotheta
      tt(2,3)=0.0
      tt(3,1)=0.0
      tt(3,2)=0.0
      tt(3,3)=-1.0
      
      tp(1,1)=cotheta
      tp(1,2)=sitheta
      tp(1,3)=0.0
      tp(2,1)=sitheta*cophip
      tp(2,2)=-cotheta*cophip
      tp(2,3)=siphip
      tp(3,1)=sitheta*siphip
      tp(3,2)=-cotheta*siphip
      tp(3,3)=-cophip
      
      tm(1,1)=cotheta
      tm(1,2)=sitheta
      tm(1,3)=0.0
      tm(2,1)=sitheta*cophip
      tm(2,2)=-cotheta*cophip
      tm(2,3)=-siphip
      tm(3,1)=-sitheta*siphip
      tm(3,2)=cotheta*siphip
      tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
      state1=0.0
      
      m(1,1)=cotheta
      m(1,2)=sitheta
      m(1,3)=0.0
      
      m(2,1)=cos(state1)*sitheta
      m(2,2)=-cos(state1)*cotheta
      m(2,3)=sin(state1)
      m(3,1)=sin(state1)*sitheta
      m(3,2)=-sin(state1)*cotheta
      m(3,3)=-cos(state1)
      
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
      
      xend(1,2)=lseg+x(1)
      xend(2,2)=x(2)
      xend(3,2)=x(3)

      do 10 i=3,long1
         
 123     rn=rands(seed)
         state=int(rn*3)
c     print*,'state',state

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
            
            call mrrrr(m,tt,mm)
            do 30 ii=1,3
               do 40 j=1,3
                  m(ii,j)=mm(ii,j)
 40            continue
 30         continue
            
            
         elseif (state.eq.1) then
            
            call mrrrr(m,tp,mm)
            do 31 ii=1,3
               do 41 j=1,3
                  m(ii,j)=mm(ii,j)
 41            continue
 31         continue

         elseif (state.eq.2) then
            
            call mrrrr(m,tm,mm)
            do 32 ii=1,3
               do 42 j=1,3
                  m(ii,j)=mm(ii,j)
 42            continue
 32         continue
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
         
c     if (xend(1,i).lt.0.0) goto 222
         
 10   continue
      
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chequea cadenas
! Selfavoiding entre segmentos
 
      dista=0.0
      do 300 ive=4,long1
         do 310 jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
            
               goto 222
            endif
 310     continue
 300  continue

      ncha=0
      do ii=1,300

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random rotation
!

      call rota36(xend,xendr,long1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Displace chains to (rpos, zpos)
! 
        
      tempp(1) = xendr(1, 1)
      tempp(2) = xendr(2, 1)
      tempp(3) = xendr(3, 1)
 
      do i=1,long1 
 
         xendr(1, i) =  xendr(1, i) + rpos - tempp(1) - 1e-4
         xendr(2, i) = xendr(2, i) - tempp(2)
         xendr(3, i) =  xendr(3, i) + zpos - tempp(3)
       
         
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
! Put segment in lattice / check coordintes      
!     
         do i = 1, long1
         do j = 1, 3
         in1(i, j) = xendr(j, i)
         xo(j) = xendr(j, i)
         enddo
         if(outside(xo)) goto 400
         enddo

         ncha=ncha+1

         do j = 1, long1
         iR = int(sqrt(in1(j,1)**2+in1(j,2)**2)/delta)+1
         iZ = int(in1(j,3)/delta)+1
         if ((iR.gt.dimR).or.(iZ.lt.1).or.(iZ.gt.dimZ)) then
         print*, 'Error in creador'
         print*, 'Increase system size'
         endif
         chains(ncha,j) = matriz(iR,iZ)
         enddo
     
         if (ncha.ge.25) goto 402
         
 400  enddo
 402  if (ncha.eq.0) goto 223

      return
      end

