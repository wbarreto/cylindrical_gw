      subroutine initial(n,m,A0,B0,a,b)

      implicit none

      include 'base.inc'

      double precision A0,B0
    

      double precision PSIEXT(2*(n+1)),dPSIEXT(2*(n+1))
      double precision OMEGAEXT(2*(n+1))

      double precision M_BASEXT(2*(n+1),2*(n+1))
      double precision copy_PSIEXT(2*(n+1)),copy_OMEGAEXT(2*(n+1))

      double precision a(2*(n+1)),sum,PSIa(2*(n+1))
      double precision b(2*(m+1)),OMEGAa(2*(n+1))
      double precision sum_error

      double precision pi,error,y

      integer n,m,k,i,j,nn

      pi=4.d0*datan(1.d0)


!..........................
! First the vector a_n(u) 
!..........................
! We get a_n(u) from the initial datum PSI0,
! projecting with the basis.
!

      PSIEXT=0.d0
      do i=1,2*n

        y=SGRID(i)

!        x=pi*(y+1.d0)/4.d0 ! Is this valid in the cylindrical case?

        PSIEXT(i)=(1.d0/8.d0)*A0*(y**3.d0)*dexp(-4.d0*(y-2.d0)**2.d0) 
!        write(*,*) y,PSIEXT(i)
      end do

      copy_PSIEXT=PSIEXT

      do i=1,2*(n+1)
        do j=1,2*(n+1)
          M_BASEXT(j,i)=SBPSI(j,i)
        end do
      end do

      call gaussj(M_BASEXT,2*(n+1),2*(n+1),PSIEXT,1,2*(1+n))

! .............Here the  vector a_n(t).............

      a=PSIEXT

!      do i=1, 2*(n+1)
!        write(*,*) i,a(i)
!      end do

!..................................................

! .............<Begin>................
! Only for checking  (just delete c as commented)

! Reconstruction of the initial PSI

c      sum_error=0.d0
c      do i=1,2*n 
c        sum=0.d0
c        do j=1,2*(n+1) 
c          sum=sum+SBPSI(i,j)*a(j)
c        end do
c        PSIa(i)=sum
c        error=(PSIa(i)-copy_PSIEXT(i))
c        sum_error=sum_error+error**2
c        write(*,*) SGRID(i),PSIa(i),copy_PSIEXT(i),error
c     end do

! Derivative of PSI with respect y...  (left to construct if necessary)

c      do i=1,2*n
c        y=SGRID(i)

c        dPIEXT(i)=...
c        sum=0.d0
c        do j=1,2*(n+1)
c          sum=sum+SBDPSI(i,j)*a(j)
c        end do
c        dPSIa(i)=sum
!        write(*,*) SGRID(i),dPSIEXT(i),dPSIa(i)
c      end do
c...................<end>........................

!............................
! Second, the vector b_n(u)
!............................

      OMEGAEXT=0.d0

      do i=1,2*n
        y=SGRID(i)
        OMEGAEXT(i)=B0*(y**3.d0)*dexp(-4.d0*(y-2.d0)**2.d0)/(1.d0+y*y)
      end do

      copy_OMEGAEXT=OMEGAEXT

      do i=1,2*(n+1)
        do j=1,2*(n+1)
          M_BASEXT(j,i)=SBOMEGA(j,i)
        end do
      end do
      call gaussj(M_BASEXT,2*(n+1),2*(n+1),OMEGAEXT,1,2*(1+n))
      
      b=OMEGAEXT

! Only for test <Begin>
! Reconstruction of OMEGA 

c      sum_error=0.d0
c     do i=1,2*n 
c       sum=0.d0
c       do j=1,2*(n+1) 
c         sum=sum+SBOMEGA(i,j)*b(j)
c       end do
c       OMEGAa(i)=sum
c       error=(OMEGAa(i)-copy_OMEGAEXT(i))
c       sum_error=sum_error+error**2
c       write(*,*) SGRID(i),OMEGAa(i),copy_OMEGAEXT(i),error
c     end do

!                <End>

      end subroutine initial
