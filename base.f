      subroutine frame(n,m,y0,L0,DIM,
     .                 SBPSI,SBDPSI,
     .                 SBOMEGA,SBDOMEGA,
     .                 SGRID)

      implicit none

      double precision y0,L0

      double precision y1n(n),y2n(n)

      double precision ji,y

      double precision BPSI1N(n,n+1),BPSI2N(n,n+1)

      double precision xi21,xi22

      double precision MPSI11(0:n),MPSI12(0:n)
      double precision MPSI21(0:n),MPSI22(0:n)

      double precision BDPSI1N(n,n+1),BDPSI2N(n,n+1)

      double precision dxi21,dxi22

      double precision BOMEGA1N(n,n+1),BOMEGA2N(n,n+1)

      double precision chi31,chi32

      double precision MOMEGA11(0:n),MOMEGA12(0:n)
      double precision MOMEGA21(0:n),MOMEGA22(0:n)

      double precision BDOMEGA1N(n,n+1),BDOMEGA2N(n,n+1)

      double precision dchi31, dchi32

      double precision SBPSI(DIM,DIM), SBDPSI(DIM,DIM)
      double precision SBOMEGA(DIM,DIM), SBDOMEGA(DIM,DIM)

      double precision SGRID(DIM)
 
      double precision pi

      integer n,DIM,m
      integer i,k,nn

      pi=4.d0*datan(1.d0)

!------------------------------------------
! Grids for each domain and each truncation
!------------------------------------------

! Grid D1 N

      do k=1,n
        ji=dcos(pi*dble(n+1-k)/dble(n+1))
        y1n(k)=0.5d0*y0*(1.d0+ji)
      end do

      do k=1,n
        ji=dcos(pi*dble(n+1-k)/dble(n+1))
        y2n(k)=y0+L0*(1.d0+ji)/(1.d0-ji)
      end do

! Super Grid

       SGRID=0.d0

       do i=1,n
         SGRID(i)=y1n(i)
         SGRID(i+n)=y2n(i)
       end do

!------------------------------------------
! Base for PSI
!------------------------------------------
! Base for PSI D1 N (BPSI1N)

      do nn=0,n
        do k=1,n
          y=y1n(k)
          BPSI1N(k,nn+1)=xi21(nn,y0,y)
        end do
      end do 

! Base for PSI D2 N (BPSI2N)

      do nn=0,n
        do k=1,n
          y=y2n(k)
          BPSI2N(k,nn+1)=xi22(nn,y0,L0,y)
        end do
      end do

!------------------------------------------
! Base for \partial_y PSI
!------------------------------------------

! Base for DPSI D1 N (BDPSI1N)

      do nn=0,n
        do k=1,n
          y=y1n(k)
          BDPSI1N(k,nn+1)=dxi21(nn,y0,y)
        end do
      end do

! Base for DPSI D2 N

      do nn=0,n
        do k=1,n
          y=y2n(k)
          BDPSI2N(k,nn+1)=dxi22(nn,y0,L0,y)
        end do
      end do

!----------------------------------------------
! Matching for PSI and DPSI
!----------------------------------------------
!
      do nn=0,n

        MPSI11(nn)=xi21(nn,y0,y0)
        MPSI12(nn)=-xi22(nn,y0,L0,y0)
        MPSI21(nn)=dxi21(nn,y0,y0)
        MPSI22(nn)=-dxi22(nn,y0,L0,y0)

      end do

!----------------------------------------------
! Super Base for PSI
!----------------------------------------------

! set to zero

      SBPSI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBPSI(k,nn+1)=BPSI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBPSI(k+n,nn+n+2)=BPSI2N(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,n
          SBPSI(2*n+1,nn+1)  = MPSI11(nn)
          SBPSI(2*n+1,nn+n+2)= MPSI12(nn)
          SBPSI(2*n+2,nn+1)  = MPSI21(nn)
          SBPSI(2*n+2,nn+n+2)= MPSI22(nn)
      end do

!----------------------------------------------
! Super Base for DPSI
!----------------------------------------------

! Set to zero

      SBDPSI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBDPSI(k,nn+1)=BDPSI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBDPSI(k+n,nn+n+2)=BDPSI2N(k,nn+1)
        end do
      end do


!------------------------------------------
! Base for OMEGA
!------------------------------------------
! Base for OMEGA D1 N (BOMEGA1N) 
 
      do nn=0,n
         do k=1,n
           y=y1n(k)
           BOMEGA1N(k,nn+1)=chi31(nn,y0,y)
         end do
      end do
  
! Base for OMEGA D2 N (BOMEGA2N)
  
      do nn=0,n
         do k=1,n
           y=y2n(k)
           BOMEGA2N(k,nn+1)=chi32(nn,y0,L0,y)
         end do
      end do

!------------------------------------------
! Base for \partial_y OMEGA
!------------------------------------------

! Base for DOMEGA D1 N (BDOMEGA1N)

      do nn=0,n
        do k=1,n
          y=y1n(k)
          BDOMEGA1N(k,nn+1)=dchi31(nn,y0,y)
        end do
      end do

! Base for DOMEGA D2 N

      do nn=0,n
        do k=1,n
          y=y2n(k)
          BDOMEGA2N(k,nn+1)=dchi32(nn,y0,L0,y)
        end do
      end do

!----------------------------------------------
! Matching for OMEGA 
!----------------------------------------------
!
      do nn=0,n

        MOMEGA11(nn) = chi31(nn,y0,y0)
        MOMEGA12(nn) =-chi32(nn,y0,L0,y0)
        MOMEGA21(nn) = dchi31(nn,y0,y0)
        MOMEGA22(nn) =-dchi32(nn,y0,L0,y0)

      end do

!----------------------------------------------
! Super Base for OMEGA
!----------------------------------------------

! set to zero

      SBOMEGA=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBOMEGA(k,nn+1)=BOMEGA1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBOMEGA(k+n,nn+n+2)=BOMEGA2N(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,n
          SBOMEGA(2*n+1,nn+1)  = MOMEGA11(nn)
          SBOMEGA(2*n+1,nn+n+2)= MOMEGA12(nn)
          SBOMEGA(2*n+2,nn+1)  = MOMEGA21(nn)
          SBOMEGA(2*n+2,nn+n+2)= MOMEGA22(nn)
      end do

!----------------------------------------------
! Super Base for DOMEGA
!----------------------------------------------

! Set to zero

      SBDOMEGA=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBDOMEGA(k,nn+1)=BDOMEGA1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBDOMEGA(k+n,nn+n+2)=BDOMEGA2N(k,nn+1)
        end do
      end do

      end subroutine frame

      double precision function T(n,y)
        implicit none
        integer n
        double precision y
        T=dcos(dble(n)*dacos(y))
      end function T

      double precision function dT(n,y)
        implicit none
        integer n
        double precision x,y
        x=dacos(y)
        if((y.eq.-1.d0).or.(y.eq.1.d0)) then
          dT=dble(n*n)*dcos(dble(n*x))/dcos(x)
        else
          dT=dble(n)*dsin(dble(n)*x)/dsin(x)
        end if
      end function dT

      double precision function TL1(n,y0,y)
        implicit none
        double precision T,y,ji,y0
        integer n
        ji=2.d0*y/y0-1.d0
        TL1=T(n,ji)
      end function TL1
  
      double precision function dTL1(n,y0,y)
        implicit none
        double precision dT,y0,y,ji
        integer n
        ji=2.d0*y/y0-1.d0
        dTL1=(2.d0/y0)*dT(n,ji)
      end function dTL1

      double precision function TL2(n,y0,L0,y)
        implicit none
        double precision T,y,ji,y0,L0
        integer n
        ji=(y-y0-L0)/(y-y0+L0)
        TL2=T(n,ji)
      end function TL2

      double precision function dTL2(n,y0,L0,y)
        implicit none
        double precision dT,y0,y,ji,L0
        integer n
        ji=(y-y0-L0)/(y-y0+L0)
        dTL2=(2.d0*L0/(y-y0+L0)**2.d0)*dT(n,ji)
      end function dTL2

      double precision function xi21(n,y0,y)
        implicit none
        integer n
        double precision TL1,y,y0
        xi21=0.5d0*(TL1(n+1,y0,y)+TL1(n,y0,y))
      end function xi21

      double precision function dxi21(n,y0,y)
        implicit none
        integer n
        double precision dTL1,y,y0
        dxi21=0.5d0*(dTL1(n+1,y0,y)+dTL1(n,y0,y))
      end function dxi21

      double precision function chi2(n,y0,y)
        implicit none
        integer n
        double precision chi0,y,y0,a,b
        a=1.d0/2.d0
        b=a*(3.d0+2.d0*n)/(1.d0+2.d0*n)
        chi2=a*chi0(n+1,y0,y)+b*chi0(n,y0,y)
      end function chi2

      double precision function dchi2(n,y0,y)
        implicit none
        integer n
        double precision dchi0,y,y0,a,b
        a=1.d0/2.d0
        b=a*(3.d0+2.d0*n)/(1.d0+2.d0*n)
        dchi2=a*dchi0(n+1,y0,y)+b*dchi0(n,y0,y)
      end function dchi2

      double precision function chi0(n,y0,y)
        implicit none
        integer n
        double precision y,y0,TL1
        chi0=0.5d0*(TL1(n+1,y0,y)+TL1(n,y0,y)) !Equal to xi21
      end function chi0

      double precision function dchi0(n,y0,y)
        implicit none
        integer n
        double precision y,y0,dTL1,a,b
        dchi0=0.5d0*(dTL1(n+1,y0,y)+dTL1(n,y0,y))
      end function dchi0

      double precision function chi31(n,y0,y)
        implicit none
        integer n
        double precision chi2,y0,y,a,b
        a=(2.d0*n*n+5.d0*n+3.d0)/(2.d0*(n+1.d0)**2+5.d0*(n+1)+3.d0)
        b=-0.25d0
        chi31=a*b*chi2(n+1,y0,y)+b*chi2(n,y0,y)
      end function chi31

      double precision function dchi31(n,y0,y)
        implicit none
        integer n
        double precision dchi2,y0,y,a,b
        a=(2.d0*n*n+5.d0*n+3.d0)/(2.d0*(n+1.d0)**2+5.d0*(n+1)+3.d0)
        b=-0.25d0
        dchi31=a*b*dchi2(n+1,y0,y)+b*dchi2(n,y0,y)
      end function dchi31

      double precision function xi22(n,y0,L0,y)
        implicit none
        integer n
        double precision TL2,y0,y,L0
        xi22=TL2(n,y0,L0,y)
      end function xi22
 
      double precision function dxi22(n,y0,L0,y)
        implicit none
        integer n
        double precision dTL2,y0,y,L0
        dxi22=dTL2(n,y0,L0,y)
      end function dxi22

      double precision function chi32(n,y0,L0,y)
        implicit none
        integer n
        double precision TL2,y0,y,L0
        chi32=TL2(n,y0,L0,y)
      end function chi32

      double precision function dchi32(n,y0,L0,y)
        implicit none
        integer n
        double precision dTL2,y0,y,L0
        dchi32=dTL2(n,y0,L0,y)
      end function dchi32
