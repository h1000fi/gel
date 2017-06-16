      SUBROUTINE splines(x,y,n,yp1,ypn,y2)
      
      
      
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)

      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then !The lower boundary condition is set either to be natural
      y2(1)=0.
      u(1)=0.
      else !or else to have a specified first derivative.
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1 !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.  
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then !The upper boundary condition is set either to be ~Snatural~T
      qn=0.
      un=0.
      else !or else to have a specified first derivative.
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1 ! This is the backsubstitution loop of the tridiagonal algorithm
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END

        double precision function Factorcurv (i)

        use mparameters
        use mgraft
        use mparameters_chain
        use mparameters_monomer

        implicit none
        integer i
        real*8 radio

        radio = posgraft(1)

        factorcurv = 1/(2*pi*(delta**3)*(dfloat(i) - 0.5d0))
        return
        end
