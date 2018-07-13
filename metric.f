      subroutine metric_A(n,ny,c,A,dA)
      implicit none
      include 'base.inc'
      integer n,ny,i,k
      double precision c(2*(n+1)), sum1, sum2, A(2*ny), dA(2*ny)
      do k=1,2*ny
        sum1=0.d0
        sum2=0.d0
        do i=1,2*(n+1)
          if(ny.eq.n) then
            sum1=sum1+SBA(k,i)*c(i)
            sum2=sum2+SBDA(k,i)*c(i)
          else
            sum1=sum1+SBA(k+2*n+2,i)*c(i)
            sum2=sum2+SBDA(k+2*n,i)*c(i)
          end if
        end do
        A(k)=1.d0+sum1
        dA(k)=sum2
      end do
      end subroutine metric_A

      subroutine metric_delta(n,ny,b,delta,ddelta)
      implicit none
      include 'base.inc'
      integer n,ny,i,k
      double precision b(2*(n+1)), sum1, sum2, delta(2*ny), ddelta(2*ny)
      do k=1,2*ny
        sum1=0.d0
        sum2=0.d0
        do i=1,2*(n+1)
          if(ny.eq.n) then
            sum1=sum1+SBD(k,i)*b(i)
            sum2=sum2+SBDD(k,i)*b(i)
          else
            sum1=sum1+SBD(k+2*n+2,i)*b(i)
            sum2=sum2+SBDD(k+2*n,i)*b(i)
          end if
        end do
        delta(k)=sum1
        ddelta(k)=sum2
      end do
      end subroutine metric_delta
