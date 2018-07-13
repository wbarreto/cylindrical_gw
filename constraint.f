      subroutine constraint(n,m,f,a,c,b)
      implicit none
      include 'base.inc'
      integer n,m,nn,k,i
      double precision f(2*n+1),a(2*(n+1)),c(2*(m+1)),b(2*(m+1))
      double precision y,x,pi2,M_dpsid(2*(m+1),2*(m+1)),cs(2*m)
      double precision Pi(2*m),dPi(2*m),Phi(2*m),dPhi(2*m)

      pi2=4.d0*datan(1.d0)

      do nn=1,2*(m+1)
        do k=1,2*m
!          M_dpsid(k,nn)=SBDD(k,nn)
          M_dpsid(k,nn)=4.d0*SBDD(k,nn)/pi2
        end do
        M_dpsid(2*m+1,nn)=SBD(2*m+1,nn)
        M_dpsid(2*m+2,nn)=SBD(2*m+2,nn)
      end do

      call field_Pi(n,m,a,Pi,dPi)
      call field_Phi(n,m,f,Phi,dPhi)

      do i=1,2*m
        y=SGRID(2*n+i)
        x=pi2*(y+1.d0)/4.d0
        cs(i)=dcos(x)*dsin(x)
      end do

!      b=-0.25d0*pi2*(Pi*Pi+Phi*Phi)*cs
      b=-(Pi*Pi+Phi*Phi)*cs
      b(2*m+1)=0.d0
      b(2*m+2)=0.d0

      call gaussj(M_dpsid,2*(m+1),2*(m+1),b,1,2*(1+m))

!      do i=1,2*(m+1)
!        write(*,*) i,b(i)
!      end do

      end subroutine constraint
