      subroutine initial(n,m,A0,B0,a,b)

      implicit none

      include 'base.inc'

      double precision A0,B0

!    <begin> <Playground zone defintions>
      double precision T,dT,ddT,x
!    <end>
    

      double precision PSIEXT(2*(n+1)),dPSIEXT(2*(n+1))
      double precision ddPSIEXT(2*(n+1))
      double precision OMEGAEXT(2*(n+1)),dOMEGAEXT(2*(n+1))
      double precision ddOMEGAEXT(2*(n+1))

      double precision M_BASEXT(2*(n+1),2*(n+1))
      double precision copy_PSIEXT(2*(n+1)),copy_OMEGAEXT(2*(n+1))

      double precision a(2*(n+1)),sum,PSIa(2*(n+1)),dPSIa(2*(n+1))
      double precision ddPSIa(2*(n+1))
      double precision b(2*(n+1)),OMEGAa(2*(n+1)), dOMEGAa(2*(n+1))
      double precision ddOMEGAa(2*(n+1))
      double precision sum_error

      double precision pi,error,y,f,df,ddf

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

        PSIEXT(i)=(1.d0/8.d0)*A0*(y**3.d0)*dexp(-4.d0*(y-2.d0)**2.d0) 
!        write(*,*) y,PSIEXT(i)
      end do

      copy_PSIEXT=PSIEXT

      PSIEXT(2*n+1)=0.d0
      PSIEXT(2*n+2)=0.d0

      do i=1,2*(n+1)
        do j=1,2*n
          M_BASEXT(j,i)=SBPSI(j,i)
        end do
      end do

      do i=1,2*(n+1)
         M_BASEXT(2*n+1,i)=SBPSI(2*n+1,i)
         M_BASEXT(2*n+2,i)=SBPSI(2*n+2,i)
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

      sum_error=0.d0
      do i=1,2*n 
        sum=0.d0
        do j=1,2*(n+1) 
          sum=sum+SBPSI(i,j)*a(j)
        end do
        PSIa(i)=sum
        error=(PSIa(i)-copy_PSIEXT(i))
        sum_error=sum_error+error**2
!        write(*,*) SGRID(i),PSIa(i),copy_PSIEXT(i),error
      end do


! Derivative of PSI with respect y... 
      sum_error=0.d0
      do i=1,2*n
        y=SGRID(i)
        dPSIEXT(i)=3.d0*copy_PSIEXT(i)/y-8.d0*(y-2.d0)*copy_PSIEXT(i)
        sum=0.d0
        do j=1,2*(n+1)
          sum=sum+SBDPSI(i,j)*a(j)
        end do
        dPSIa(i)=sum
        error=(dPSIa(i)-dPSIEXT(i))
        sum_error=sum_error+error**2
!        write(*,*) SGRID(i),dPSIEXT(i),dPSIa(i),error
      end do

! Derivative of PSI with respect y twice ... 
      sum_error=0.d0
      do i=1,2*n
        y=SGRID(i)
        ddPSIEXT(i)=(3.d0/y)*(dPSIEXT(i)-copy_PSIEXT(i)/y)
     .             -8.d0*(copy_PSIEXT(i)+(y-2.d0)*dPSIEXT(i))
        sum=0.d0
        do j=1,2*(n+1)
          sum=sum+SBDDPSI(i,j)*a(j)
        end do
        ddPSIa(i)=sum
        error=(ddPSIa(i)-ddPSIEXT(i))
        sum_error=sum_error+error**2
!        write(*,*) SGRID(i),ddPSIEXT(i),ddPSIa(i),error
      end do
c...................<end>........................

!............................
! Second, the vector b_n(u)
!............................

      OMEGAEXT=0.d0

      do i=1,2*n
        y=SGRID(i)
        OMEGAEXT(i)=B0*(y**3.d0)*dexp(-4.d0*(y-2.d0)**2.d0)/(1.d0+y*y)
      end do

      OMEGAEXT(2*n+1)=0.d0
      OMEGAEXT(2*n+1)=0.d0

      copy_OMEGAEXT=OMEGAEXT

      do i=1,2*(n+1)
        do j=1,2*n
          M_BASEXT(j,i)=SBOMEGA(j,i)
        end do
      end do

      do i=1,2*(n+1)
         M_BASEXT(2*n+1,i)=SBOMEGA(2*n+1,i)
         M_BASEXT(2*n+2,i)=SBOMEGA(2*n+2,i)
      end do
      
      call gaussj(M_BASEXT,2*(n+1),2*(n+1),OMEGAEXT,1,2*(1+n))
      
      b=OMEGAEXT

! Only for test <Begin>
! Reconstruction of OMEGA 

      sum_error=0.d0
      do i=1,2*n 
       sum=0.d0
       do j=1,2*(n+1) 
         sum=sum+SBOMEGA(i,j)*b(j)
       end do
       OMEGAa(i)=sum
       error=(OMEGAa(i)-copy_OMEGAEXT(i))
       sum_error=sum_error+error**2
!       write(*,*) SGRID(i),OMEGAa(i),copy_OMEGAEXT(i),error
      end do

! Derivative of OMEGA with respect y... 

      OMEGAEXT=copy_OMEGAEXT

      sum_error=0.d0
      do i=1,2*n
        y=SGRID(i)
        f=y**3.d0/(1.d0+y**2.d0)
        df=y*y*(3.d0+y*y)/(1+y*y)**2.d0
        dOMEGAEXT(i)=(df/f-8.d0*(y-2.d0))*OMEGAEXT(i)
        sum=0.d0
        do j=1,2*(n+1)
          sum=sum+SBDOMEGA(i,j)*b(j)
        end do
        dOMEGAa(i)=sum
        error=(dOMEGAa(i)-dOMEGAEXT(i))
        sum_error=sum_error+error**2
!        write(*,*) SGRID(i),dOMEGAEXT(i),dOMEGAa(i),error
      end do

! Derivative of OMEGA with respect y twice ... 
      sum_error=0.d0
      do i=1,2*n
        y=SGRID(i)
        f=y**3.d0/(1.d0+y**2.d0)
        df=y*y*(3.d0+y*y)/(1+y*y)**2.d0
        ddf=6.d0*y/(1.d0+y*y)-14.d0*y*y*y/(1.d0+y*y)**2.d0
     .      +8.d0*y**5.d0/(1.d0+y*y)**3.d0
!        ddOMEGAEXT(i)=(ddf/f-8.d0-(df/f)**2.d0)*OMEGAEXT(i)
!     .               +dOMEGAEXT(i)*(df/f-8.d0*(y-2.d0))
        ddOMEGAEXT(i)=OMEGAEXT(i)*( (ddf+16.d0*df*(2.d0-y))/f
     .               +64.d0*y*(y-4.d0)+248.d0  )
        sum=0.d0
        do j=1,2*(n+1)
          sum=sum+SBDDOMEGA(i,j)*b(j)
        end do
        ddOMEGAa(i)=sum
        error=(ddOMEGAa(i)-ddOMEGAEXT(i))
        sum_error=sum_error+error**2
!        write(*,*) SGRID(i),ddOMEGAEXT(i),ddOMEGAa(i),error
      end do

!                <End>

!     <begin> <Playground zone>

!         do nn=29,29
!           do i=-90,90
!             x=dble(i)/1.d2
!             write(*,*) x,T(nn,x),dT(nn,x),ddT(nn,x)
!           end do
!           write(*,*)
!         end do

!     <end>

      end subroutine initial
