      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)


c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
c xai~Rs in order), and given the array y2a(1:n), which is the output from spline above,
cand given a value of x, this routine returns a cubic-spline interpolated value y.
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
c We will find the right place in the table by means of bisection.
cThis is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely
c spaced, one would do better to store previous values of
c klo and khi and test if they remain appropriate on the
c next call.
      khi=n
    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      ENDIF ! klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
      if (h.eq.0.0)stop 
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     & ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END


