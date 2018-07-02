      subroutine initial(n,m,A0,sigma,S,r0,a,b,c,f,y0)

      implicit none

      include 'base.inc'

      double precision A0,sigma,S,r0,y0
    
      double precision M_TL1(n+1,n+1),M_psi2(n+1,n+1)
      double precision dM_TL1(n+1,n+1),dM_psi2(n+1,n+1)
      double precision copy1(n+1,n+2),copy2(n+1,n+1)
      double precision T,PI01(n+1),PI02(n+1)
      double precision TL1,TL2,dTL1,dT,dTL2
      double precision copyPI01(n+1),copypI02(n+1)
      double precision sum1(n+1),sum2(m+1)
      double precision M_psi(n+1,n+1),M_psinv(n+1,n+1)

      double precision M_dpsid1(m+1,m+1),M_dpsid2(m+1,m+1)

      double precision M_dpsid(m+1,m+1),M_dpsidinv(m+1,m+1)
      double precision M_psia(m+1,m+1),M_dpsiainv(m+1,m+1)
    
      double precision M_psia1(m+1,m+1),M_dpsia1(m+1,m+1)
      double precision M_psia2(m+1,m+1),M_dpsia2(m+1,m+1)

      double precision copy(n+1,n+1),M_dpsia(m+1,m+1)
      double precision copy2delta(m+1,m+1)
      double precision copyPI0(n+1)
      double precision psi2,y,x,identity(n+1,n+1)
      double precision dpsi2
      double precision constr11(n+1),constr12(n+1)
      double precision constr21(n+1),constr22(n+1)

      double precision constr11delta(m+1),constr12delta(m+1)
      double precision constr21delta(m+1),constr22delta(m+1)


      double precision constr11A(m+1),constr12A(m+1)
      double precision constr21A(m+1),constr22A(m+1)

      double precision PSIEXT(2*(n+1)),dPIEXT(2*(n+1))

      double precision M_BASEXT(2*(n+1),2*(n+1)),PI0(2*(n+1))
!      double precision dM_BASEXT(2*(n+1),2*(n+1))
      double precision copy_M_BASEXT(2*(n+1),2*(n+1))
!      double precision SBPI(DIM,DIM)
!      double precision SBDPI(DIM,DIM)
      double precision copy_PSIEXT(2*(n+1))

      double precision M_EXTENDELTA(2*(m+1),2*(m+1))
      double precision M_EXTENDA(2*(m+1),2*(m+1))

      double precision TERMEXT(2*(m+1))
      double precision temp1ext(2*(m+1)),temp2ext(2*(m+1))

      double precision psi_a,dpsi_a
      double precision psia1,dpsia1
      double precision a(2*(n+1)),sum,PSIa(2*(n+1))
      double precision b(2*(m+1)),c(2*(m+1))
      double precision dPIa(2*(n+1))
      double precision f(2*n+1)
      double precision cs(m+1),term1(m+1),copyterm(m+1),sum_error
      double precision term2(m+1)
      double precision PI01delta(m+1),PI02delta(m+1)
      double precision PIEXTD(2*m),csext(2*m),s2ext(2*m)
      double precision PI01a(m+1),PI02a(m+1)
      double precision s2(m+1)
      double precision temp1,temp2

      double precision recons((n+1)*10,n+1), PIa_rec((n+1)*10) 
      double precision PI0_rec((n+1)*10)

      double precision pi,ji,error

      integer n,m,k,i,j,nn

      pi=4.d0*datan(1.d0)


!..........................
! First the vector a_n(t) 
!..........................
! We get a_n(t) from the initial datum PSI0,
! projecting with the basis.
!
      PSIEXT=0.d0

      do i=1,2*n
        y=SGRID(i)
        x=pi*(y+1.d0)/4.d0
        PSIEXT(i)=(1.d0/8.d0)*A0*(x**3.d0)*dexp(-4.d0*(x-2.d0)**2.d0) 
!        write(*,*) x,PSIEXT(i)
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

!      do i=1, 2*(n+1)  ! same than Maple script for A0=0.3, L0=3.5=y0?
!        write(*,*) i,a(i)
!      end do

!..................................................

! .............<Begin>................
! For checking  (just delete c as commented)

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
        write(*,*) SGRID(i),PSIa(i),copy_PSIEXT(i),error
      end do
!      write(*,*) n,sum_error
!      write(*,*) '........ '

! Derivative of PI with respect y.

c      do i=1,2*n
c        y=SGRID(i)
c        x=pi*(y+1.d0)/4.d0

c        dPIEXT(i)=-8.d0*pi*dtan(x)*
c     .            (1.d0+dtan(x)**2.d0)*A0*
c     .            dexp(-16.d0*dtan(x)**2.d0)
c        sum=0.d0
c        do j=1,2*(n+1)
c          sum=sum+SBDPI(i,j)*a(j)
c        end do
c        dPIa(i)=sum
!        write(*,*) SGRID(i),dPIEXT(i),dPIa(i)
c      end do
c...................<end>........................

!............................
! Second, the vector f_n(t)
!............................
!      f=0.d0    ! f has dimension 2n+1 (one collocation point less in D2)
!..........................

!..........................
! Third, the vector b_n(t) 
!..........................
! For this we use the constrain equation,
! taking advantage of the initial datum PI0
! and the fact that in our case Phi0=0.
!
! Domain 1

!      do nn=0,m
!        do k=1,m
!          ji=dcos(pi*dble(k)/dble(m+1))
!          y=0.5d0*(ji*(1.d0+y0)-1.d0+y0)
!          M_dpsid1(k,nn+1)=dpsia1(nn,y0,y)
!        end do
!      end do

!      PI01delta=0.d0

!      do i=1,m
!        ji=dcos(pi*dble(i)/dble(m+1))
!        y=0.5d0*(ji*(1.d0+y0)+y0-1.d0)
!        x=pi*(y+1.d0)/4.d0
!        PI01delta(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
!        cs(i)=dcos(x)*dsin(x)
!      end do

!      term1=-pi*0.25d0*PI01delta*PI01delta*cs

! Domain 2

!      do nn=0,m
!        if(nn.eq.0) then
!         sum2(nn)=.0d0
!        else
!         sum2(nn)=-dTL2(nn,y0,1.d0)/dTL2(m+1,y0,1.d0)
!        end if
!        do k=1,m
!          ji=dcos(pi*dble(k)/dble(m+1))
!          y=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
!          if (nn.eq.0) then
!            M_dpsid2(k,nn+1)=0.d0
!          else
!            M_dpsid2(k,nn+1)=dTL2(nn,y0,y)+sum2(nn)*dTL2(m+1,y0,y) 
!          end if
!        end do
!      end do

!      PI02delta=0.d0

!      do i=1,m
!        ji=dcos(pi*dble(i)/dble(m+1))
!        y=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
!        x=pi*(y+1.d0)/4.d0
!        PI02delta(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
!        cs(i)=dcos(x)*dsin(x)
!      end do

!      term2=-pi*0.25d0*PI02delta*PI02delta*cs
      
! Matching conditions between D1 and D2

!      do nn=0,m

!        constr11delta(nn)=psia1(nn,y0,y0)
!        constr12delta(nn)=-(TL2(nn,y0,y0)+sum2(nn)*TL2(m+1,y0,y0))

!        constr21delta(nn)=dpsia1(nn,y0,y0)

!        if (nn.eq.0) then
!          constr22delta(nn)=0.d0
!        else
!          constr22delta(nn)=-(dTL2(nn,y0,y0)+sum2(nn)*dTL2(m+1,y0,y0))
!        end if
!        write(*,*) constr21delta(nn),constr22delta(nn)
!      end do

! Extended issues...

!      M_EXTENDELTA=0.d0

!      do nn=0,m
!        do k=1,m
!           M_EXTENDELTA(k,nn+1)= M_dpsid1(k,nn+1)
!        end do
!      end do

!      do nn=0,m
!        do k=1,m
!           M_EXTENDELTA(k+m,nn+m+2)=M_dpsid2(k,nn+1)
!        end do
!      end do


!      do nn=0,m
!          M_EXTENDELTA(2*m+1,nn+1)  = constr11delta(nn)
!          M_EXTENDELTA(2*m+1,nn+m+2)= constr12delta(nn)
!          M_EXTENDELTA(2*m+2,nn+1)  = constr21delta(nn)
!          M_EXTENDELTA(2*m+2,nn+m+2)= constr22delta(nn)
!      end do

!      TERMEXT=0.d0

!      do k=1,m   
!        TERMEXT(k)=term1(k)
!        TERMEXT(k+m)=term2(k)
!      end do

!      do i=1,2*m
!        write(*,*) i, TERMEXT(i)
!      end do



! All the business in what follows....................

!      TERMEXT=0.d0

!      do i=1,2*m
!        y=SGRID(2*n+i)
!        x=pi*(y+1.d0)/4.d0
!        PIEXTD(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
!        csext(i)=dcos(x)*dsin(x)
!      end do 

!!      TERMEXT=-pi*0.25d0*PIEXTD*PIEXTD*csext
!      TERMEXT=-PIEXTD*PIEXTD*csext  !!! OJO: asi es como coincide con el script
                                    !!!  de Maple (????)

!      do i=1,2*m
!        write(*,*) i,TERMEXT(i)
!      end do

!      M_EXTENDELTA=0.d0

!      do j=1,2*(m+1)
!        do i=1,2*m
!!          M_EXTENDELTA(i,j)=SBDD(i,j) 
!          M_EXTENDELTA(i,j)=4.d0*SBDD(i,j)/pi  ! This makes sense...
!        end do
!      end do
   
!      write(*,*) (M_EXTENDELTA(2,i),i=1,2*(m+1))


!      do nn=1,2*(m+1)
!         M_EXTENDELTA(2*m+1,nn)=SBD(2*m+1,nn)
!         M_EXTENDELTA(2*m+2,nn)=SBD(2*m+2,nn)
!      end do

!      do i=1,2*(m+1)
!        write(*,*) M_EXTENDELTA(2*m+1,i)
!      end do

! I have to prove that this way of business works......
! (done) the output using the two business models is the same.
! Tudo bem. 

!      call gaussj(M_EXTENDELTA,2*(m+1),2*(m+1),TERMEXT,1,2*(1+m)) 

!.................... 
!      b=TERMEXT
!....................

!      do i=1,2*(m+1)
!       write(*,*) i,b(i)
!      end do

!...........................
! Fourth, the vector c_n(t)
!...........................

!Domain 1

c      do nn=0,m
c        do k=1,m
c          ji=dcos(pi*dble(k)/dble(m+1))
c          y=0.5d0*(ji*(1.d0+y0)-1.d0+y0)
c          M_psia1(k,nn+1)=psia1(nn,y0,y)
c          M_dpsia1(k,nn+1)=dpsia1(nn,y0,y)
c        end do
c      end do
 
c      PI01a=0.d0

c      do i=1,m
c          ji=dcos(pi*dble(i)/(m+1))
c          y=0.5d0*(ji*(1.d0+y0)-1.d0+y0)
c          x=pi*(y+1.d0)/4.d0
c          PI01a(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
c          cs(i)=dcos(x)*dsin(x)
c          s2(i)=1.d0+2.d0*dsin(x)**2.d0
c      end do

!      do i=1,m
c      do i=1,m+1
c        do j=1,m
c          temp1=pi*0.25d0*PI01a(j)*PI01a(j)*cs(j)
c          temp2=pi*0.25d0*s2(j)/cs(j)
c          M_dpsia1(j,i)=M_dpsia1(j,i)
c     .                 +M_psia1(j,i)*(temp1+temp2)
c        end do
c      end do

c      term1=-pi*0.25d0*PI01a*PI01a*cs

! Domain 2

c      do nn=0,m
c        do k=1,m
c          ji=dcos(pi*dble(k)/dble(m+1))
c          y=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
c          M_psia2(k,nn+1)=psi2(nn,y0,y)
c          M_dpsia2(k,nn+1)=dpsi2(nn,y0,y)
c        end do
c      end do


c      PI02a=0.d0

c      do i=1,m
c          ji=dcos(pi*dble(i)/dble(m+1))
c          y=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
c          x=pi*(y+1.d0)/4.d0
c          PI02a(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
c          cs(i)=dcos(x)*dsin(x)
c          s2(i)=1.d0+2.d0*dsin(x)**2.d0
c      end do

!      do i=1,m
c      do i=1,m+1
c        do j=1,m
c          temp1=pi*0.25d0*PI02a(j)*PI02a(j)*cs(j)
c          temp2=pi*0.25d0*s2(j)/cs(j)
c          M_dpsia2(j,i)=M_dpsia2(j,i)
c     .                 +M_psia2(j,i)*(temp1+temp2)
c        end do
c      end do

c      term2=-pi*0.25d0*PI02a*PI02a*cs

! Matching conditions on y_0 between D1 and D2

c      do nn=0,m

c        constr11A(nn)=psia1(nn,y0,y0)
c        constr12A(nn)=-psi2(nn,y0,y0)
c        constr21A(nn)=dpsia1(nn,y0,y0)
c        constr22A(nn)=-dpsi2(nn,y0,y0)
c      end do

! Extension issues...

c      M_EXTENDA=0.d0

c      do nn=0,m
c        do k=1,m
c           M_EXTENDA(k,nn+1)= M_dpsia1(k,nn+1)
c        end do
c      end do

c      do nn=0,m
c        do k=1,m
c           M_EXTENDA(k+m,nn+m+2)=M_dpsia2(k,nn+1)
c        end do
c      end do

c      do nn=0,m
c          M_EXTENDA(2*m+1,nn+1)  = constr11A(nn)
c          M_EXTENDA(2*m+1,nn+m+2)= constr12A(nn)
c          M_EXTENDA(2*m+2,nn+1)  = constr21A(nn)
c          M_EXTENDA(2*m+2,nn+m+2)= constr22A(nn)
c      end do

c      write(*,*) (M_EXTENDA(2*m+1,i),i=1,2*(m+1))
c      write(*,*) "..."

c      TERMEXT=0.d0

c      do k=1,m
c        TERMEXT(k)=term1(k)
c        TERMEXT(k+m)=term2(k)
c      end do
!----------------------------------------
! Other way for business...

!      TERMEXT=0.d0

!      do i=1,2*m
!        y=SGRID(2*n+i)
!        x=pi*(y+1.d0)/4.d0
!        PIEXTD(i)=A0*dexp(-dtan(x)**2.d0/sigma**2.d0)
!        csext(i)=dcos(x)*dsin(x)
!        s2ext(i)=1.d0+2.d0*dsin(x)**2.d0
!      end do

!      TERMEXT=-pi*0.25d0*PIEXTD*PIEXTD*csext 

!      temp1ext=-TERMEXT
!      temp2ext=pi*0.25d0*s2ext/csext

!      M_EXTENDA=0.d0

!      do j=1,2*(m+1)
!        do i=1,2*m
!          M_EXTENDA(i,j)=SBDA(i,j)+SBA(i,j)*(temp1ext(i)+temp2ext(i))
!        end do
!      end do
  
!      do nn=1,2*(m+1)
!         M_EXTENDA(2*m+1,nn)=SBA(2*m+1,nn)
!         M_EXTENDA(2*m+2,nn)=SBA(2*m+2,nn)
!      end do

!      write(*,*) (SBA(2*m+1,i),i=1,2*(m+1))

!      call gaussj(M_EXTENDA,2*(m+1),2*(m+1),TERMEXT,1,2*(1+m)) 

!.................... 
!      c=TERMEXT
!....................

!      do i=1,2*(m+1)  ! igual que el script de Maple para A0=11
!        write(*,*) i,c(i)
!      end do

      end subroutine initial
