      program ads

      implicit none

      integer itimemax,itime,idump,iref

      integer n,m,i,j

      include 'base.inc'

      double precision timemax,h,time,A0,sigma,S,r0
      double precision a(NMAX),b(NMAX),c(NMAX),f(NMAX)
      double precision da(NMAX),dc(NMAX),df(NMAX)
      double precision amr,y0

      double precision btemp(NMAX)


      integer fp_seq, fp_data, fp_asymp
 
      logical done
     
      namelist/input/n,m,A0,sigma,S,r0,timemax,
     .               h,itime,idump,iref,amr,y0

      read (*,input) 

      m=n*13/10

      fp_data=21
      fp_seq=28
      fp_asymp=35

      open(unit=fp_seq,file='seq.gnu',status='unknown')
      open(unit=fp_asymp,file='asymp.gr',status='unknown')

      a=0.d0
      b=0.d0
      c=0.d0
      f=0.d0
      da=0.d0
      dc=0.d0
      df=0.d0

      SBPI=0.d0
      SBDPI=0.d0
      SBPHI=0.d0
      SBDPHI=0.d0
      SBA=0.d0
      SBDA=0.d0
      SBD=0.d0
      SBDD=0.d0
      SGRID=0.d0

      call frame(n,m,y0,NMAX,
     .           SBPI,SBDPI,
     .           SBPHI,SBDPHI,
     .           SBA,SBDA,
     .           SBD,SBDD,
     .           SGRID)

!      write(*,*) " after frame.."

!      do i=1,2*(m+1)
!        write(*,*) i, SBD(2*m+1,i),SBD(2*m+2,i)
!      end do

      call initial(n,m,A0,sigma,S,r0,a,b,c,f,y0)

      btemp=b

!      write(*,*) " after initial.."

      call write_files(fp_data,fp_seq,fp_asymp,
     .                 n,m,a,f,c,b,0.d0,itime,iref,done)

!      write(*,*) "after write..."

      itimemax=timemax/h


      done=.false.

! La siguiente llamada es solo para verificar que
! esta subrutina hace lo esperado... se debe eliminar
! una vez verificado esto.

!      call constraint(n,m,f,a,c,b)

!      do i =1,2*(m+1)
!      write(*,*) i,dabs(btemp(i)-b(i))
!      end do 
! Verificado, está dentro del error. 


      do while (.not. done)

!       write(*,*) time,a(1),f(1),b(1),c(1) 

       time=dble(itime+1)*h

        
       call dynsys(n,m,f,b,a,c,df,da,dc)
!       write(*,*) "after dynsys..."

       call rk4(n,m,h,b,f,a,c,df,da,dc,f,a,c)
!       write(*,*) "after rk..."

       call constraint(n,m,f,a,c,b)
!       write(*,*) "after constraint..."

       if(time.ge.timemax) done=.true.

       itime=itime+1

        if(mod(itime,idump).eq.0) then
          call write_files(fp_data,fp_seq,fp_asymp,
     .                     n,m,a,f,c,b,time,itime,iref,done)
        end if

      end do
 
!      call write_files(fp_data,fp_seq,fp_asymp,
!     .                 n,m,a,f,c,b,time,itime,iref,done)

      close(fp_seq)
      close(fp_asymp)

      end     

