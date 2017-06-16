      subroutine rota36(xend,xendr,n)
      
      use  mparameters
      use mparameters_chain
      use mparameters_monomer
      use mrands

      implicit none
 
      real*8 xend(3,maxlong),xendr(3,maxlong)
      character*1 test
      integer n
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      
      real*8 radio, vect
      real*8 alfa, gama, cga, a, b, c
      integer i
 

      
      fac=rands(seed)
      fac1=rands(seed)
      fac2=rands(seed)
      alfa=fac*2*pi
      cbe=2.0d0*fac1-1.0d0
      gama=fac2*2*pi

      sbe=(1-cbe**2)**0.5
      cal=cos(alfa)
      sal=sin(alfa)
      cga=cos(gama)
      sga=sin(gama)

      do 1 i=1,n

         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)
     &        -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+
     &        b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
 
 1    continue

      return
      end
