      subroutine write_files(fp_data,fp_seq,fp_asymp,
     .                       n,m,a,f,c,b,time,itime,iref,done)
      implicit none
      include 'base.inc'
      integer fp_data, fp_seq, fp_asymp
      logical done,detectnan
      integer i,n,itime,m,iref
      double precision a(2*(n+1)),f(2*n+1),c(2*(m+1)),b(2*(m+1))
      double precision AA(2*n*iref),dAA(2*n*iref)

      double precision time,Pi(2*n),dPi(2*n)
      double precision Phi(2*n),dPhi(2*n),iPhi(2*n)
      double precision mass(2*n)
      double precision delta(2*n),ddelta(2*n),dotiPhi(2*n)
      double precision pi2,y(2*n*iref),x(2*n*iref),ADMass
      character filename*11
      
      pi2=4.d0*datan(1.d0)

      write (filename(1:7),'(i7)') itime + 1000000
      filename(1:1) = 'p'
      filename(8:11) = '.dat'

      open (unit = fp_data, file = filename,  status = 'unknown')

      write (fp_seq, '(a15,e21.13,a10,a11,a5)')
     &      'set title "t = ', time, '" ; plot "',
     &       filename, '" w l'
!      call flush (fp_seq)

      call field_Pi(n,n,a,Pi,dPi)
!      call kodis(Pi,2*n)

      call field_Phi(n,n,f,Phi,dPhi)
!      call field_iPhi(n-1,n,f,iPhi)

!      write(*,*) "before metric in io.f..."
      call metric_A(m,n*iref,c,AA,dAA)
!      write(*,*) "after metric in io.f...."

      call metric_delta(m,n,b,delta,ddelta)

!      iPhi=4.d0*iPhi/pi2
!      write(*,*) "before using sgrid..."
      do i=1,2*n*iref
        y(i)=SGRID(i)
        x(i)=pi2*(y(i)+1.d0)/4.d0
      end do
!      write(*,*) "after using sgrid..."

!      iPhi(1)=0.d0
!      do i=1,2*n-1
!        iPhi(i+1)=iPhi(i)+0.5d0*(Phi(i)+Phi(i+1))/dble(n)
!      end do

!      do i=n,2*n-1
!        iPhi(i+1)=iPhi(i)+0.5d0*(Phi(i)+Phi(i+1))/dble(n)
!      end do

      iPhi(2*n)=0.d0
      do i=2*n,2,-1
        iPhi(i-1)=iPhi(i)+0.5d0*(Phi(i)+Phi(i-1))*(y(i-1)-y(i))
      end do
!      do i=n,1,-1
!        iPhi(i+1)=iPhi(i)+0.5d0*(Phi(i)+Phi(i+1))/dble(n)
!      end do

      iPhi=pi2*iPhi/4.d0

      mass=0.5d0*dtan(x)*(1.d0-AA)/dcos(x)**2.d0

      ADMass=mass(2*n)

! calculates

      dotiPhi=AA*dexp(-delta)*Pi

      write(*,*) time, minval(AA), 0.25d0*pi2*(1.d0+y(minloc(AA)))
      detectnan=isnan(minval(AA))
      if (detectnan.eqv..true.) then
         done=.true.
      end if

      write (fp_data, '(2e21.13)')
     &    (y(i),AA(i), i = 1,2*n*iref )
!      write (fp_data, '(2e21.13)')
!     &    (y(i),iPhi(i), i = 1,2*n*iref )


!      write (fp_asymp, '(5e21.13)') 
!     &    time,iPhi(2),dotiPhi(2),ADMass,y(2)

      write (fp_asymp, '(5e21.13)') 
     &    time,Phi(1),Pi(1),iPhi(1),dotiPhi(1)

      close (unit = fp_data)

      end subroutine write_files
