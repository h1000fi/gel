      subroutine mrrrr(a,b,c)

      implicit none
      
      real*8 a(3,3),b(3,3),c(3,3)

      integer i, j, k

      do 1 i=1,3
         do 1 j=1,3
            c(i,j)=0
 1    continue

      do 2 i=1,3
         do 2 j=1,3
            do 2 k=1,3
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2    continue

      return
      end
