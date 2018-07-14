      subroutine metric_PSI(n,ny,a,PSI,dPSI,ddPSI)
      implicit none
      include 'base.inc'
      integer n,ny,i,k
      double precision a(2*(n+1)), sum1, sum2, PSI(2*ny), dPSI(2*ny),
     .                 ddPSI(2*ny),sum3

!  Domain 1
      do k=1,n
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        do i=1,n+1
            sum1=sum1+SBPSI(k,i)*a(i)
            sum2=sum2+SBDPSI(k,i)*a(i)
            sum3=sum3+SBDDPSI(k,i)*a(i)
        end do
        PSI(k)=sum1
        dPSI(k)=sum2
        ddPSI(k)=sum3
      end do
! Domain 2
      do k=n+1,2*n
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        do i=n+2,2*(n+1)
            sum1=sum1+SBPSI(k,i)*a(i)
            sum2=sum2+SBDPSI(k,i)*a(i)
            sum3=sum3+SBDDPSI(k,i)*a(i)
        end do
        PSI(k)=sum1
        dPSI(k)=sum2
        ddPSI(k)=sum3
      end do
      end subroutine metric_PSI

      subroutine metric_OMEGA(n,ny,b,OMEGA,dOMEGA,ddOMEGA)
      implicit none
      include 'base.inc'
      integer n,ny,i,k
      double precision b(2*(n+1)),sum1,sum2,OMEGA(2*ny),dOMEGA(2*ny),
     .                 ddOMEGA(2*ny),sum3

! I tested that the split of two domain it is not necessary here...
! Because the way in which the super base is constructed (using the
! mask trick).

! Domain 1
      do k=1,n
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        do i=1,n+1
            sum1=sum1+SBOMEGA(k,i)*b(i)
            sum2=sum2+SBDOMEGA(k,i)*b(i)
            sum3=sum3+SBDDOMEGA(k,i)*b(i)
        end do
        OMEGA(k)=sum1
        dOMEGA(k)=sum2
        ddOMEGA(k)=sum3
      end do
! Domain 2
      do k=n+1,2*n
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        do i=n+2,2*(n+1)
            sum1=sum1+SBOMEGA(k,i)*b(i)
            sum2=sum2+SBDOMEGA(k,i)*b(i)
            sum3=sum3+SBDDOMEGA(k,i)*b(i)
        end do
        OMEGA(k)=sum1
        dOMEGA(k)=sum2
        ddOMEGA(k)=sum3
      end do
      end subroutine metric_OMEGA
