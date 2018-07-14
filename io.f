      subroutine write_files(fp_data,fp_seq,fp_asymp,
     .                       n,m,a,b,time,itime,iref,done)
      implicit none
      include 'base.inc'
      integer fp_data, fp_seq, fp_asymp
      logical done,detectnan
      integer i,n,itime,m,iref
      double precision a(2*(n+1)),b(2*(n+1))
      double precision PSI(2*n*iref),dPSI(2*n*iref),ddPSI(2*n)
      double precision time,OMEGA(2*n),dOMEGA(2*n),ddOMEGA(2*n)
      double precision mass(2*n)
      double precision pi,y(2*n*iref),ADMass
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
     &    (y(i),OMEGA(i), i = 1,2*n*iref )

      close (unit = fp_data)

      end subroutine write_files
