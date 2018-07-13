      subroutine field_Pi(n,ny,a,Pi,dPi)
      implicit none
      include 'base.inc'
      integer n,ny,i,k,j
      double precision a(2*(n+1)), sum1, sum2
      double precision Pi(2*ny), dPi(2*ny)
      do k=1,2*ny
        sum1=0.d0
        sum2=0.d0
        do i=1,2*(n+1)
          if(ny.eq.n) then 
            sum1=sum1+SBPI(k,i)*a(i)
            sum2=sum2+SBDPI(k,i)*a(i)
          else
            sum1=sum1+SBPI(k+2*n+2,i)*a(i)
            sum2=sum2+SBDPI(k+2*n,i)*a(i)
          end if
        end do
        Pi(k)=sum1
        dPi(k)=sum2
      end do
      end subroutine field_Pi

      subroutine field_Phi(n,ny,f,Phi,dPhi)
      implicit none
      include 'base.inc'
      integer n,ny,i,k
      double precision f(2*n+1), sum1, sum2, Phi(2*ny), dPhi(2*ny)
      do k=1,2*ny
        sum1=0.d0
        sum2=0.d0
        do i=1,2*n+1
          if(ny.eq.n) then
            sum1=sum1+SBPHI(k,i)*f(i)
            sum2=sum2+SBDPHI(k,i)*f(i)
          else
            sum1=sum1+SBPHI(k+2*n+1,i)*f(i)
            sum2=sum2+SBDPHI(k+2*n,i)*f(i)
          end if
        end do
        Phi(k)=sum1
        dPhi(k)=sum2
      end do
      end subroutine field_Phi

!      subroutine field_iPhi(n,ny,f,iPhi)
!      implicit none
!      integer n,ny,i,k
!      double precision f(n+1), sum1, iPhi(ny+1)
!      double precision ichi,pi,y
!      pi=4.d0*datan(1.d0)
!      do k=1,ny+1
!        sum1=0.d0
!        y=dcos(dble(k)*pi/(ny+2))
!        do i=1,n+1
!          sum1=sum1+f(i)*ichi(i-1,y)
!        end do
!        iPhi(k)=sum1
!      end do
!      end subroutine field_iPhi
