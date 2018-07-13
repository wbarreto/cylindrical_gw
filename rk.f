      subroutine rk4(n,m,delta,b,x,y,z,dx,dy,dz,xout,yout,zout)

      implicit none

      double precision x(2*n+1),y(2*(n+1)),z(2*(m+1)),b(2*(m+1))
      double precision dx(2*n+1),dy(2*(n+1)),dz(2*(m+1))
      double precision xout(2*n+1),yout(2*(n+1)),zout(2*(m+1))
      double precision xt(2*n+1),yt(2*(n+1)),zt(2*(m+1))
      double precision dxt(2*n+1),dyt(2*(n+1)),dzt(2*(m+1))
      double precision dxm(2*n+1),dym(2*(n+1)),dzm(2*(m+1))
      double precision hh,delta,h6
      integer i,n,m

      hh=delta*.5d0
      h6=delta/6.d0

      do i=1,2*n+1
        xt(i)=x(i)+hh*dx(i)
      end do

      do i=1,2*(n+1)
        yt(i)=y(i)+hh*dy(i)
      end do

      do i=1,2*(m+1)
        zt(i)=z(i)+hh*dz(i)
      end do

      call dynsys(n,m,xt,b,yt,zt,dxt,dyt,dzt)
!      write(*,*) "after dynsys in rk... call 1"

      do i=1,2*n+1
        xt(i)=x(i)+hh*dxt(i)
      end do

      do i=1,2*(n+1)
        yt(i)=y(i)+hh*dyt(i)
      end do

      do i=1,2*(m+1)
        zt(i)=z(i)+hh*dzt(i)
      end do

      call dynsys(n,m,xt,b,yt,zt,dxm,dym,dzm)
!      write(*,*) "after dynsys in rk... call 2"

      do i=1,2*n+1
        xt(i)=x(i)+delta*dxm(i)
        dxm(i)=dxt(i)+dxm(i)
      end do

      do i=1,2*(n+1)
        yt(i)=y(i)+delta*dym(i)
        dym(i)=dyt(i)+dym(i)
      end do

      do i=1,2*(m+1)
        zt(i)=z(i)+delta*dzm(i)
        dzm(i)=dzt(i)+dzm(i)
      end do

      call dynsys(n,m,xt,b,yt,zt,dxt,dyt,dzt)
!      write(*,*) "after dynsys in rk... call 3"

      do i=1,2*n+1
        xout(i)=x(i)+h6*(dx(i)+dxt(i)+2.d0*dxm(i))
      end do

      do i=1,2*(n+1)
        yout(i)=y(i)+h6*(dy(i)+dyt(i)+2.d0*dym(i))
      end do

      do i=1,2*(m+1)
        zout(i)=z(i)+h6*(dz(i)+dzt(i)+2.d0*dzm(i))
      end do

      end subroutine rk4 
