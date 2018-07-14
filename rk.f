      subroutine rk4(n,m,L0,delta,x,y,dx,dy,xout,yout)

      implicit none

      double precision x(2*(n+1)),y(2*(n+1))
      double precision dx(2*(n+1)),dy(2*(n+1))
      double precision xout(2*(n+1)),yout(2*(n+1))
      double precision xt(2*(n+1)),yt(2*(n+1))
      double precision dxt(2*(n+1)),dyt(2*(n+1))
      double precision dxm(2*(n+1)),dym(2*(n+1))
      double precision hh,delta,h6,L0
      integer i,n,m

      hh=delta*.5d0
      h6=delta/6.d0

      do i=1,2*(n+1)
        xt(i)=x(i)+hh*dx(i)
        yt(i)=y(i)+hh*dy(i)
      end do

      call dynsys(n,m,L0,xt,yt,dxt,dyt)

      do i=1,2*(n+1)
        xt(i)=x(i)+hh*dxt(i)
        yt(i)=y(i)+hh*dyt(i)
      end do

      call dynsys(n,m,L0,xt,yt,dxm,dym)

      do i=1,2*(n+1)
        xt(i)=x(i)+delta*dxm(i)
        dxm(i)=dxt(i)+dxm(i)
        yt(i)=y(i)+delta*dym(i)
        dym(i)=dyt(i)+dym(i)
      end do

      call dynsys(n,m,L0,xt,yt,dxt,dyt)

      do i=1,2*(n+1)
        xout(i)=x(i)+h6*(dx(i)+dxt(i)+2.d0*dxm(i))
        yout(i)=y(i)+h6*(dy(i)+dyt(i)+2.d0*dym(i))
      end do

      end subroutine rk4 
