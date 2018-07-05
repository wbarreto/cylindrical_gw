      subroutine dynsys(n,m,a,b,da,db)
        implicit none
        include 'base.inc'
        integer n,m,k,l,i,j
        double precision a(2*(n+1)),b(2*(n+1))
        double precision da(2*(n+1)),db(2*n+1)

        double precision PSI(2*n),OMEGA(2*n),y(2*n)
        double precision dPSI(2*n),dOMEGA(2*n)
        double precision ddPSI(2*n),ddOMEGA(2*n)

        double precision RHS1(2*n), RHS2(2*n) 
        double precision M_EXT(4*(n+1),4*(n+1))
        double precision RHS_EXT(4*(n+1))
 
        double precision pi

        pi=4.d0*datan(1.d0)
!
!   Dynamical system as one...
!
        call metric_PSI(n,n,a,PSI,dPSI,ddPSI)
        call metric_OMEGA(n,n,b,OMEGA,dOMEGA,ddOMEGA)

!   First evolution equation
        
        y=SGRID
        RHS1=0.25d0*(ddPSI-dPSI/y+PSI/y**2.d0)
        RHS2=-0.125d0*dexp(4.d0*PSI/y)*(OMEGA+y*dOMEGA)**2.d0/y**3.d0
       
        do i=1,2*n 
          RHS_EXT(i)=RHS1(i)+RHS2(i)
        end do

        M_EXT=0.d0

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(j,i)=y(j)*SBDPSI(j,i)
          end do
        end do

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(j,i+2*n+3)=(-dexp(4.d0*PSI(j)/y(j))
     .                  *(OMEGA(j)+y(j)*dOMEGA(j))
     .                  /(2.d0*y(j)))*SBOMEGA(j,i)
          end do
        end do

! In similar way follows the second evolution equation

        RHS1=0.25d0*y**2.d0*(-3.d0*(OMEGA+y*dOMEGA)/y**4.d0+
     .       (2.d0*dOMEGA+y*ddOMEGA)/y**3.d0)
        RHS2=(OMEGA+y*dOMEGA)*(dPSI/y-PSI/y**2.d0)/y

        do i=1,2*n
          RHS_EXT(i+2*n)=RHS1(i)+RHS2(i)
        end do

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(2*n+j,i)=2.d0*(y(j)*dOMEGA(j)+OMEGA(j))*SBPSI(j,i)
          end do
        end do

        do i=1,2*(n+1)
          do j=1,2*n
             M_EXT(2*n+j,i+2*n+3)=y(j)*SBDOMEGA(j,i)
     .                            +2.d0*y(j)*(dPSI(j)/y(j)-PSI(j)
     .                            /y(j)**2.d0)*SBOMEGA(j,i)
          end do
        end do

!...............................

        RHS_EXT(4*n+1)=0.d0
        RHS_EXT(4*n+2)=0.d0
        RHS_EXT(4*n+3)=0.d0
        RHS_EXT(4*n+4)=0.d0
 
        call gaussj(M_EXT,4*(n+1),4*(n+1),RHS_EXT,1,4*(n+1))

        do i=1,2*(n+1)
          da(i)=RHS_EXT(i)
        end do
        do i=1,2*(n+1)
          db(i)=RHS_EXT(i+2*(n+1))
        end do

      end subroutine dynsys
