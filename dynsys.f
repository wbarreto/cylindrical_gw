      subroutine dynsys(n,m,f,b,a,c,df,da,dc)
        implicit none
        include 'base.inc'
        integer n,m,k,l,i,j
        double precision a(2*(n+1)),f(2*n+1),c(2*(m+1))
        double precision da(2*(n+1)),df(2*n+1),dc(2*(m+1))
        double precision b(2*(m+1))

        double precision Pi_I(2*n),AA_I(2*n),delta_I(2*n)
        double precision dPi_I(2*n),dAA_I(2*n),ddelta_I(2*n)

        double precision Phi_II(2*n), dPhi_II(2*n)
        double precision AA_II(2*n), dAA_II(2*n)
        double precision delta_II(2*n),ddelta_II(2*n)

        double precision AA_III(2*m), dA_III(2*m)
        double precision delta_III(2*m), ddelta_III(2*m)
        double precision Phi_III(2*m), dPhi_III(2*m)
        double precision Pi_III(2*m), dPi_III(2*m)

        double precision RHS1(2*n+1),M_chi(2*n+1,2*n+1)
        double precision tanterm(2*n),RHS2(2*(n+1)) 
        double precision M_psi(2*(n+1),2*(n+1))
        double precision RHS3(2*(m+1)),sc(2*m)
        double precision M_psia(2*(m+1),2*(m+1))
 
        double precision pi2

        double precision y, x
        double precision chi,psi
        double precision psi_A

!        double precision ji,y0

!        y0=-0.9d0  ! ojo con este parámetro aquí
!                    debe pasar como un argumento
!                    en la subrutina.

        pi2=4.d0*datan(1.d0)
!
!.......
! First dynamical equation 
!
! We have to be carefull with dimensions here; for that
! reason we choose label this set of equations with I.
!

        call metric_A(m,n,c,AA_I,dAA_I)
        call metric_delta(m,n,b,delta_I,ddelta_I)
        call field_Pi(n,n,a,Pi_I,dPi_I)
!
!..........................................
! Este bloque en principio no hace falta,
! es para la escritura del test.
!
!        do i=1,2*n
!        if (i.le.n) then
!          ji=dcos(pi2*dble(i)/dble(n+1)) !ojo ya esta secuencia no se usa
!          y=0.5d0*(ji*(1.d0+y0)+y0-1.d0)
!        else
!          ji=dcos(pi2*dble(1+i-(n+1))/dble(n+1))
!          y=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
!        end if
!          write(*,*) y, Pi_I(i)
!        end do
!..........................................

        RHS1=4.d0*dexp(-delta_I)*(dAA_I*Pi_I
     .            -ddelta_I*AA_I*Pi_I+AA_I*dPi_I)/pi2

!        write(*,*) "RHS1(2*n+1)=", RHS1(2*n+1)

        RHS1(2*n+1)=0.d0

        do i=1,2*n+1
          do j=1,2*n+1
             M_chi(j,i)=SBPHI(j,i)
          end do
        end do

! Ojo, confirmar la siguiente sintaxis....
! que sustituye lo anterior. Es decir, un arreglo
! de dimension menor puede igualarse a uno de dimension
! mayor.  Y para el caso de un arreglo mayor igualado
! a uno de dimensión menor los elementos fuera del rango
! se hacen cero....

!        M_chi=SBPHI

!........................................
  
c        do l=0,n-1
c          do k=1,n
c            y=dcos(pi2*dble(k)/dble(n+1))
c            M_chi(k,l+1)=chi(l,y) ! Which does not change in time
c          end do                  ! then can be placed out of here
c        end do

!        write(*,*) (M_chi(i,2),i=1,n)
       
        call gaussj(M_chi,2*n+1,2*n+1,RHS1,1,2*n+1)

        df=RHS1

!        do i=1,2*n+1
!          write(*,*) i, df(i),SGRID(i) 
!        end do
!
!......
! Second dynamical equation (for dot a)
!
        call metric_A(m,n,c,AA_II,dAA_II)
        call metric_delta(m,n,b,delta_II,ddelta_II)
        call field_Phi(n,n,f,Phi_II,dPhi_II)

! The next one does not depend on time
!
         do i=1,2*n
            y=SGRID(i)
            x=(1.d0+y)*pi2/4.d0
            tanterm(i)=2.d0*(1.d0+dtan(x)**2.d0)/dtan(x)
         end do

!         write(*,*) (tanterm(i),i=1,2*n)
!

        RHS2=(4.d0/pi2)*dexp(-delta_II)*(dAA_II*Phi_II
     .       -ddelta_II*AA_II*Phi_II
     .       +AA_II*dPhi_II) + 
     .       dexp(-delta_II)*AA_II*Phi_II*tanterm

        RHS2(2*n+1)=0.d0
        RHS2(2*n+2)=0.d0

c        do l=0,n
c          do k=1,n+1
c            y=dcos(pi2*dble(k)/dble(n+2))
c            M_psi(k,l+1)=psi(l,y) ! Which does not change in time
c          end do                  ! then can be placed out of here
c        end do

!        write(*,*) (M_psi(i,n+1),i=1,n+1)

         do i=1,2*(n+1)
           do j=1,2*(n+1)
             M_psi(j,i)=SBPI(j,i)
           end do
         end do

        call gaussj(M_psi,2*(n+1),2*(n+1),RHS2,1,2*(1+n))

        da=RHS2

!        do i=1,n+1
!          write(*,*) i,da(i)
!        end do

!
!......
! Third dynamical equation

        call metric_A(m,m,c,AA_III,dA_III)
        call metric_delta(m,m,b,delta_III,ddelta_III)
        call field_Pi(n,m,a,Pi_III,dPi_III)
        call field_Phi(n,m,f,Phi_III,dPhi_III)

        do i=1,2*m
           y=SGRID(2*n+i)
           x=(1.d0+y)*pi2/4.d0
           sc(i)=dsin(x)*dcos(x)
        end do

        RHS3 = -2.d0 * AA_III**2.d0 * dexp(-delta_III) *
     .         Pi_III * Phi_III * sc

        RHS3(2*m+1)=0.d0
        RHS3(2*m+2)=0.d0

c        do l=0,m
c          do k=1,m+1
c            y=dcos(pi2*dble(k)/dble(m+2))
c            M_psia(k,l+1)=psi_A(l,y) ! Which does not change in time
c          end do                     ! then can be placed out of here
c        end do

         do i=1,2*(m+1)
           do j=1,2*(m+1)
             M_psia(j,i)=SBA(j,i)
           end do
         end do

        call gaussj(M_psia,2*(m+1),2*(m+1),RHS3,1,2*(1+m))

        dc=RHS3

!        do i=1,2*(m+1)
!          write(*,*) i, dc(i),SGRID(i+2*n) 
!        end do

      end subroutine dynsys
