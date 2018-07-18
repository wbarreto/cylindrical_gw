      subroutine write_files(fp_data,fp_seq,
     .                       n,m,y0,L0,a,b,time,itime,iref,done)
      implicit none
      include 'base.inc'
      integer fp_data, fp_seq, fp_asymp
      logical done,detectnan
      integer i,n,itime,m,iref
      double precision a(2*(n+1)),b(2*(n+1))
      double precision PSI(2*n*iref),dPSI(2*n*iref),ddPSI(2*n)
      double precision time,OMEGA(2*n),dOMEGA(2*n),ddOMEGA(2*n)
      double precision pi,y(2*n*iref)
      double precision GAMMA(2*n),bm,y0,L0
      double precision da(2*(n+1)),db(2*(n+1))
      double precision news1,news2


      character filename*11
      
      pi=4.d0*datan(1.d0)

      write (filename(1:7),'(i7)') itime + 1000000
      filename(1:1) = 'p'
      filename(8:11) = '.dat'

      open (unit = fp_data, file = filename,  status = 'unknown')

      write (fp_seq, '(a15,e21.13,a10,a11,a5)')
     &      'set title "t = ', time, '" ; plot "',
     &       filename, '" w l'
      call flush (fp_seq)

      call metric_PSI(n,n*iref,a,PSI,dPSI,ddPSI)

      call metric_OMEGA(n,n,b,OMEGA,dOMEGA,ddOMEGA)

      do i=1,2*n*iref
        y(i)=SGRID(i)
      end do

      write (fp_data, '(2e21.13)')
     &    (y(i),PSI(i), i = 1,2*n*iref )

      close (unit = fp_data)

      end subroutine write_files

      subroutine write_asymp (n,y0,L0,a,b,fp_asymp,time)
        implicit none
        include 'base.inc'
        integer i,j,n,fp_asymp
        double precision y(2*n)
        double precision a(2*(n+1)),PSI(2*n),dPSI(2*n),ddPSI(2*n)
        double precision b(2*(n+1)),OMEGA(2*n),dOMEGA(2*n),
     .                   ddOMEGA(2*n)
        double precision x(2*n),y0,L0,dydx(2*n),dx
        double precision GAMMAy(2*n),GAMMAx(2*n),GAMMAxh
        double precision GAMMA(2*n),GAMMAa,bm

        double precision da(2*(n+1)),db(2*(n+1))
        double precision SUM1,SUM2,SUM3,SUM4
        double precision PSIa,OMEGAa,PSIua,OMEGAua,ya
        double precision news1,news2,time

        double precision dda(2*(n+1)),ddb(2*(n+1))
        double precision GAMMAua,bmloss
        double precision dPSIa,dOMEGAa

        double precision PSIuu(2*n),OMEGAuu(2*n),PSIuua,OMEGAuua
        double precision wsre,wsim

        call metric_PSI(n,n,a,PSI,dPSI,ddPSI)
        call metric_OMEGA(n,n,b,OMEGA,dOMEGA,ddOMEGA)

        do i=1,2*n
          y(i)=SGRID(i)
          if (y(i).le.y0) then
            x(i)=2.d0*y(i)/y0-1.d0
            dydx(i)=0.5d0*y0
          else
            x(i)=(y(i)-y0-L0)/(y(i)-y0+L0)
            dydx(i)=2.d0*L0/(1.d0-x(i))**2.d0
          end if
        end do

        GAMMAy=0.5d0*y*(dPSI/y-PSI/y**2.d0)**2.d0
     .        +dexp(4.d0*PSI/y)*(y*dOMEGA+OMEGA)**2.d0/(8.d0*y**3.d0)

        GAMMAx=dydx*GAMMAy

! Integration of GAMMA using finite differences----
        GAMMA(1)=0.d0
        do i=2,2*n
           if (i.eq.n+1) then
             dx=x(2)-x(1)
           else
             dx=x(i)-x(i-1)
           end if
           GAMMAxh=0.5d0*(GAMMAx(i)+GAMMAx(i-1))
           GAMMA(i)=GAMMA(i-1) + GAMMAxh*dx
        end do

        GAMMAa=GAMMA(2*n)

! Bondi mass---------------
        bm=0.5d0*GAMMAa
!--------------------------

        call dynsys(n,n,L0,a,b,da,db)

        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        sum4=0.d0
        do i=1,2*(n+1)
           sum1=sum1 + SBPSI(2*n,i)*da(i)
           sum2=sum2 + SBOMEGA(2*n,i)*db(i)
           sum3=sum3 + SBPSI(2*n,i)*a(i)
           sum4=sum4 + SBOMEGA(2*n,i)*b(i)
        end do
        ya=SGRID(2*n)
        PSIua=sum1
        OMEGAua=sum2
        PSIa=sum3
        OMEGAa=sum4
        dPSIa=dPSI(2*n)
        dOMEGAa=dOMEGA(2*n)
! Here we can use an alternative for the Bondi mass loss---
        GAMMAua=PSIua*(dPSIa/ya-PSIa/ya**2.d0)-2.d0*PSIua*PSIua
     .         +dexp(4.d0*PSIa/ya)*(OMEGAua*(dOMEGAa*ya+OMEGAa)
     .         -2.d0*ya*ya*OMEGAua*OMEGAua)/(4.d0*ya*ya) 
        bmloss=0.5d0*GAMMAua
! News functions -------------------------------------
        news1=PSIua
        news2=0.5d0*dexp(2.d0*PSIa/ya)*
     .        (2.d0*PSIua*OMEGAa/ya+OMEGAua)
!-------------------------------------------------

        call dynsys(n,n,L0,da,db,dda,ddb)
        do i=1,2*n
          sum1=0.d0
          sum2=0.d0
          do j=1,2*(n+1)
             sum1=sum1+SBPSI(i,j)*dda(j)
             sum2=sum2+SBOMEGA(i,j)*ddb(j)
          end do
          PSIuu(i)=sum1
          OMEGAuu(i)=sum2
        end do
        PSIuua=PSIuu(2*n)
        OMEGAuua=OMEGAuu(2*n)

! Weyl scalar---------------------------------------------------------
        wsre=dexp(-2.d0*GAMMAa)*(2.d0*PSIua*GAMMAua-PSIuua)
        wsim=0.5d0*dexp(-2.d0*GAMMAa)*(-2.d0*OMEGAua*GAMMAua+OMEGAuua)
!---------------------------------------------------------------------
        write(fp_asymp,'(7e21.13)') time,bm,news1,news2,bmloss,wsre,wsim

      end subroutine write_asymp

