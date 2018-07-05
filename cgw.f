      program cgw
!---------------------------.-------------------------------------------------
! Rio axial code to simulate cylindrical gravitational waves
! Referential code: ads with dd
! Referential maple script provided by HdO by 09.06.18
! Initial date: 29.06.18
! Final date (First run with certified output from the Maple Script): 
! Updates: 
!-----------------------------------------------------------------------------

      implicit none

      integer itimemax,itime,idump,iref

      integer n,m,i,j

      include 'base.inc'

      double precision timemax,h,time,A0,B0,S,r0
      double precision a(NMAX),b(NMAX)
      double precision da(NMAX),db(NMAX)
      double precision amr,y0,L0

      double precision btemp(NMAX)


      integer fp_seq, fp_data, fp_asymp
 
      logical done
     
      namelist/input/n,m,A0,B0,S,r0,timemax,
     .               h,itime,idump,iref,amr,y0,L0

      read (*,input) 

!      fp_data=21
!      fp_seq=28
!      fp_asymp=35

!      open(unit=fp_seq,file='seq.gnu',status='unknown')
!      open(unit=fp_asymp,file='asymp.gr',status='unknown')

      a=0.d0
      b=0.d0
      da=0.d0
      db=0.d0

      SBPSI=0.d0
      SBDPSI=0.d0
      SBDDPSI=0.d0
      SBOMEGA=0.d0
      SBDOMEGA=0.d0
      SBDDOMEGA=0.d0
      SGRID=0.d0

      call frame(n,m,y0,L0,NMAX,
     .           SBPSI,SBDPSI,SBDDPSI,
     .           SBOMEGA,SBDOMEGA,SBDDOMEGA,
     .           SGRID)


!      do i=1,2*(m+1)
!        write(*,*) i, SBD(2*m+1,i),SBD(2*m+2,i)
!      end do

      call initial(n,m,A0,B0,a,b)

!      btemp=b

!      call write_files(fp_data,fp_seq,fp_asymp,
!     .                 n,m,a,f,c,b,0.d0,itime,iref,done)


      itimemax=timemax/h

      done=.false.

      do while (.not. done)

       time=dble(itime+1)*h

        
       call dynsys(n,m,a,b,da,db)

       call rk4(n,m,h,a,b,da,db,a,b)

       if(time.ge.timemax) done=.true.

       itime=itime+1

!        if(mod(itime,idump).eq.0) then
!          call write_files(fp_data,fp_seq,fp_asymp,
!     .                     n,m,a,f,c,b,time,itime,iref,done)
!        end if

      end do
 
!      close(fp_seq)
!      close(fp_asymp)

      end program cgw    

