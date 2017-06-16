      subroutine fitspl (yspl, xspl, mulist, points)

      implicit none
           
      integer points, nspl, i, n
      parameter (nspl=1000)
      REAL*8 yspl(nspl), xspl(nspl), yp1, ypn
      REAL*8 ireal, nsplreal
      REAL*8 mulist(2,nspl), x(nspl), y(nspl), y2(nspl)  ! Dimensiona por exceso

      n = points - 1

      do i = 1, n
      x(i) = mulist(1, i)
      y(i) = mulist(2, i)
      end do

      yp1 = 1E31
      ypn = 1E31
      call splines(x, y, n, yp1, ypn, y2)

      nsplreal = nspl
      do  i = 1, nspl
      ireal = i
      xspl(i) = ((ireal-1)/(nsplreal-1))*(x(n)-x(1))+x(1)
      call splint(x, y, y2, n, xspl(i), yspl(i))
      end do
      end

