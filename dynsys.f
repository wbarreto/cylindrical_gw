      subroutine dynsys(n,m,L0,a,b,da,db)
        implicit none
        include 'base.inc'
        integer n,m,k,l,i,j
        double precision a(2*(n+1)),b(2*(n+1))
        double precision da(2*(n+1)),db(2*(n+1))

        double precision PSI(2*n),OMEGA(2*n),y(2*n)
        double precision dPSI(2*n),dOMEGA(2*n)
        double precision ddPSI(2*n),ddOMEGA(2*n)

        double precision RHS1(2*n), RHS2(2*n) 
        double precision M_EXT(4*(n+1),4*(n+1))
        double precision RHS_EXT(4*(n+1))
 
        double precision pi,L0

        pi=4.d0*datan(1.d0)
!
!   Dynamical system as one...
!
        call metric_PSI(n,n,a,PSI,dPSI,ddPSI)
        call metric_OMEGA(n,n,b,OMEGA,dOMEGA,ddOMEGA)

!   First evolution equation
!       y=SGRID ! does not work because SGRID has dimension NMAX and y has dim 2*n

!   What follows (with exception of y(i) and RHS_EXT(i)) was written dummy, i.e.,
!      RHS1=0.25d0*(ddPSI-dPSI/y+PSI/y**2.d0), because all the 
!   objects have the same dimension.

        do i=1,2*n 

          y(i)=SGRID(i)

          RHS1(i)=0.25d0*(ddPSI(i)-dPSI(i)/y(i)+PSI(i)/y(i)**2.d0)
      
          RHS2(i)=-0.125d0*dexp(4.d0*PSI(i)/y(i))*((OMEGA(i)
     .         +y(i)*dOMEGA(i))**2.d0)/(y(i)**3.d0)

          RHS_EXT(i)=RHS1(i)+RHS2(i)

        end do

!        write(*,'(i3,e16.8)') (i,RHS_EXT(i),i=1,2*n)

        M_EXT=0.d0


        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(j,i)=y(j)*SBDPSI(j,i)
          end do
        end do

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(j,i+2*(n+1))=(-dexp(4.d0*PSI(j)/y(j))
     .                  *(OMEGA(j)+y(j)*dOMEGA(j))
     .                  /(2.d0*y(j)))*SBOMEGA(j,i)
          end do
        end do

!        do j=1,2*n
!           write(*,*) (M_EXT(j,i),i=1,4*(n+1))
!        write(*,*) "----------"
!        end do

! In similar way follows the second evolution equation

        do i=1,2*n
          RHS1(i)=0.25d0*(ddOMEGA(i)-dOMEGA(i)/y(i)
     .            -3.d0*OMEGA(i)/y(i)**2.d0)
          RHS2(i)=(OMEGA(i)+y(i)*dOMEGA(i))*(dPSI(i)/y(i)
     .            -PSI(i)/y(i)**2.d0)/y(i)
          RHS_EXT(i+2*n)=RHS1(i)+RHS2(i)
        end do

!        write(*,'(i3,e16.8)') (i+2*n,RHS_EXT(i+2*n),i=1,2*n)

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(2*n+j,i)=2.d0*(y(j)*dOMEGA(j)+OMEGA(j))
     .                      *SBPSI(j,i)/y(j)
          end do
        end do

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(2*n+j,i+2*(n+1))=y(j)*SBDOMEGA(j,i)
     .                            +2.d0*y(j)*(dPSI(j)/y(j)-PSI(j)
     .                            /y(j)**2.d0)*SBOMEGA(j,i)
          end do
        end do

!        do j=2*n+1,4*n
!           write(*,*) (M_EXT(j,i),i=1,4*(n+1))
!           write(*,*) "----- "
!        end do

        RHS_EXT(4*n+1)=0.d0
        RHS_EXT(4*n+2)=0.d0
        RHS_EXT(4*n+3)=0.d0
        RHS_EXT(4*n+4)=0.d0

        do i=1,2*(n+1)
           M_EXT(4*n+1,i)=SBPSI(2*n+1,i) 
           M_EXT(4*n+2,i)=L0*SBPSI(2*n+2,i) 
           M_EXT(4*n+3,2*n+2+i)=SBOMEGA(2*n+1,i)
           M_EXT(4*n+4,2*n+2+i)=L0*SBOMEGA(2*n+2,i)
        end do
!        write(*,*) (i,M_EXT(4*n+4,2*n+2+i),i=1,2*(n+1))

!        write(*,'(5e10.2)') (M_EXT(4*n+1,i),i=1,4*(n+1))
!        print*," "
!        write(*,'(5e10.2)') (M_EXT(4*n+2,i),i=1,4*(n+1))
!        print*," "
!        write(*,'(5e10.2)') (M_EXT(4*n+3,i),i=1,4*(n+1))
!        print*," "
!        write(*,'(5e10.2)') (M_EXT(4*n+4,i),i=1,4*(n+1))

        call gaussj(M_EXT,4*(n+1),4*(n+1),RHS_EXT,1,4*(n+1))

        do i=1,2*(n+1)
          da(i)=RHS_EXT(i)
        end do
        do i=1,2*(n+1)
          db(i)=RHS_EXT(i+2*(n+1))
        end do

      end subroutine dynsys
