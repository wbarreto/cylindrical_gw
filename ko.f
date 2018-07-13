      subroutine kodis(f,n)
      implicit none
      integer i,n
      double precision f(n),fdis(n)
      double precision dis,fact,a,b,c,shape
      dis=2.0d0 !aqui estaba en 0.1
      fact = dis/16.d0
      do i = 3,n-2
         fdis(i) = f(i-2)-4.d0*f(i-1)+6.d0*f(i)-4.d0*f(i+1)+f(i+2)
      end do
      a = 4
      b = n - 3
      c = 256.d0/(a-b)**8
      do i = 3, n-2
         shape = c * (i-a)**4 * (i-b)**4
         f(i) = f(i) - fact * fdis(i) * shape
      end do
      end subroutine
