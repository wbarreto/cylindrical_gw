      subroutine frame(n,m,y0,DIM,
     .                 SBPI,SBDPI,
     .                 SBPHI,SBDPHI,
     .                 SBA,SBDA,
     .                 SBD,SBDD,
     .                 SGRID)

      implicit none

      double precision y0
      double precision y1n(n),y2n(n),y1m(m),y2m(m)

      double precision mistery1m(m)

      double precision ji,y

      double precision TL1,psi2,dTL1,dpsi2
      double precision sum1(0:n),sum2(0:m)
      double precision chi2,dchi2
      double precision TL2,dTL2,psia1,dpsia1

      double precision BPI1N(n,n+1),BPI2N(n,n+1)
      double precision BPI1M(m,n+1),BPI2M(m,n+1)

      double precision MPI11(0:n),MPI12(0:n)
      double precision MPI21(0:n),MPI22(0:n)

      double precision BDPI1N(n,n+1),BDPI2N(n,n+1)
      double precision BDPI1M(m,n+1),BDPI2M(m,n+1)

      double precision BPHI1N(n,n+1),BPHI2N(n,n)
      double precision BPHI1M(m,n+1),BPHI2M(m,n)

      double precision MPHI11(0:n),MPHI12(0:n-1)

      double precision BDPHI1N(n,n+1),BDPHI2N(n,n)
      double precision BDPHI1M(m,n+1),BDPHI2M(m,n)

      double precision BA1N(n,m+1),BA2N(n,m+1)
      double precision BA1M(m,m+1),BA2M(m,m+1)

      double precision MA11(0:m),MA12(0:m)
      double precision MA21(0:m),MA22(0:m)

      double precision BDA1N(n,m+1),BDA2N(n,m+1)
      double precision BDA1M(m,m+1),BDA2M(m,m+1)

      double precision BD1N(n,m+1),BD2N(n,m+1)
      double precision BD1M(m,m+1),BD2M(m,m+1)      

      double precision MD11(0:m),MD12(0:m)
      double precision MD21(0:m),MD22(0:m)

      double precision BDD1N(n,m+1),BDD2N(n,m+1)
      double precision BDD1M(m,m+1),BDD2M(m,m+1)      

      double precision SBPI(DIM,DIM), SBDPI(DIM,DIM)
      double precision SBPHI(DIM,DIM), SBDPHI(DIM,DIM)
      double precision SBA(DIM,DIM), SBDA(DIM,DIM)
      double precision SBD(DIM,DIM), SBDD(DIM,DIM)

      double precision copySBD(DIM,DIM)

      double precision SGRID(DIM)
 
      double precision pi

      integer n,m,DIM
      integer i,k,nn

      pi=4.d0*datan(1.d0)

!------------------------------------------
! Grids for each domain and each truncation
!------------------------------------------

! Grid D1 N

      do k=1,n
!        ji=dcos(pi*dble(k)/dble(n+1))
        ji=dcos(pi*dble(n+1-k)/dble(n+1))
        y1n(k)=0.5d0*(ji*(1.d0+y0)-1.d0+y0) 
      end do

! Grid D1 M

      do k=1,m
!        ji=dcos(pi*dble(k)/dble(m+1))
        ji=dcos(pi*dble(m+1-k)/dble(m+1))
        y1m(k)=0.5d0*(ji*(1.d0+y0)-1.d0+y0) 
!        mistery1m(k)=0.5d0*(ji*(1.d0+y0)-1.d0+y0)
      end do

! Grid D2 N

      do k=1,n
!        ji=dcos(pi*dble(k)/dble(n+1))
        ji=dcos(pi*dble(n+1-k)/dble(n+1))
        y2n(k)=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
      end do

! Grid D2 M

      do k=1,m
!        ji=dcos(pi*dble(k)/dble(m+1))
        ji=dcos(pi*dble(m+1-k)/dble(m+1))
        y2m(k)=0.5d0*(ji*(1.d0-y0)+1.d0+y0)
      end do

! Super Grid

       SGRID=0.d0

       do i=1,n
         SGRID(i)=y1n(i)
         SGRID(i+n)=y2n(i)
       end do
       do i=1,m
         SGRID(2*n+i)=y1m(i)
         SGRID(2*n+m+i)=y2m(i)
       end do

!------------------------------------------
! Base for PI
!------------------------------------------
! Base for PI D1 N (BPI1N)

      do nn=0,n
        if (nn.eq.0) then
         sum1(nn)=0.d0
        else
         sum1(nn)=-dTL1(nn,y0,-1.d0)/dTL1(n+1,y0,-1.d0)
        end if
        do k=1,n
          y=y1n(k)
          BPI1N(k,nn+1)=TL1(nn,y0,y)+sum1(nn)*TL1(n+1,y0,y)
        end do
      end do 

! Base for PI D2 N (BPI2N)

      do nn=0,n
        do k=1,n
          y=y2n(k)
          BPI2N(k,nn+1)=psi2(nn,y0,y)
        end do
      end do

! Base for PI D1 M

      do nn=0,n
        if (nn.eq.0) then
         sum1(nn)=0.d0
        else
         sum1(nn)=-dTL1(nn,y0,-1.d0)/dTL1(n+1,y0,-1.d0)
        end if
        do k=1,m
          y=y1m(k)
          BPI1M(k,nn+1)=TL1(nn,y0,y)+sum1(nn)*TL1(n+1,y0,y)
        end do
      end do

! Base for PI D2 M

      do nn=0,n
        do k=1,M
          y=y2m(k)
          BPI2M(k,nn+1)=psi2(nn,y0,y)
        end do
      end do

!------------------------------------------
! Base for \partial_y PI
!------------------------------------------

! Base for DPI D1 N (BDPI1N)

      do nn=0,n
        do k=1,n
          y=y1n(k)
          if(nn.eq.0) then
            BDPI1N(k,nn+1)=0.d0
          else
            BDPI1N(k,nn+1)=dTL1(nn,y0,y)+sum1(nn)*dTL1(n+1,y0,y)
          end if
        end do
      end do

! Base for DPI D2 N

      do nn=0,n
        do k=1,n
          y=y2n(k)
          BDPI2N(k,nn+1)=dpsi2(nn,y0,y)
        end do
      end do

! Base for DPI D1 M (BDPI1M)

      do nn=0,n
        do k=1,m
          y=y1m(k)
          if(nn.eq.0) then
            BDPI1M(k,nn+1)=0.d0
          else
            BDPI1M(k,nn+1)=dTL1(nn,y0,y)+sum1(nn)*dTL1(n+1,y0,y)
          end if
        end do
      end do

! Base for DPI D2 M

      do nn=0,n
        do k=1,M
          y=y2m(k)
          BDPI2M(k,nn+1)=dpsi2(nn,y0,y)
        end do
      end do

!----------------------------------------------
! Matching for PI and DPI
!----------------------------------------------
!
      do nn=0,n

        MPI11(nn)=TL1(nn,y0,y0)+sum1(nn)*TL1(n+1,y0,y0)
        MPI12(nn)=-psi2(nn,y0,y0)

        if (nn.eq.0) then
          MPI21(nn)=0.d0
        else
          MPI21(nn)= dTL1(nn,y0,y0)+sum1(nn)*dTL1(n+1,y0,y0)
        end if

        MPI22(nn)=-dpsi2(nn,y0,y0)

      end do

!----------------------------------------------
! Super Base for PI
!----------------------------------------------

! set to zero

      SBPI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBPI(k,nn+1)=BPI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBPI(k+n,nn+n+2)=BPI2N(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,n
          SBPI(2*n+1,nn+1)  = MPI11(nn)
          SBPI(2*n+1,nn+n+2)= MPI12(nn)
          SBPI(2*n+2,nn+1)  = MPI21(nn)
          SBPI(2*n+2,nn+n+2)= MPI22(nn)
      end do

! Sector D1 M

      do nn=0,n
        do k=1,m
           SBPI(k+2*n+2,nn+1)=BPI1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,n
        do k=1,m
           SBPI(2*n+2+m+k,nn+n+2)=BPI2M(k,nn+1)
        end do
      end do

!----------------------------------------------
! Super Base for DPI
!----------------------------------------------

! Set to zero

      SBDPI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBDPI(k,nn+1)=BDPI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n
        do k=1,n
           SBDPI(k+n,nn+n+2)=BDPI2N(k,nn+1)
        end do
      end do

! Sector D1 M

      do nn=0,n
        do k=1,m
           SBDPI(k+2*n,nn+1)=BDPI1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,n
        do k=1,m
           SBDPI(2*n+m+k,nn+n+2)=BDPI2M(k,nn+1)
        end do
      end do

!------------------------------------------

!------------------------------------------
! Base for PHI
!------------------------------------------
! Base for PHI D1 N (BPHI1N) 
 
      do nn=0,n
         if (nn.eq.0) then
           sum1(nn)=0.d0
         else
           sum1(nn)=-dTL1(nn,y0,-1.d0)/dTL1(n+1,y0,-1.d0)
         end if
         do k=1,n
           y=y1n(k)
           BPHI1N(k,nn+1)=TL1(nn,y0,y)+sum1(nn)*TL1(n+1,y0,y)
         end do
      end do
  
! Base for PHI D2 N (BPHI2N)
  
      do nn=0,n-1
         do k=1,n
           y=y2n(k)
           BPHI2N(k,nn+1)=chi2(nn,y0,y)
         end do
      end do

! Base for PHI D1 M

      do nn=0,n
        if (nn.eq.0) then
         sum1(nn)=0.d0
        else
         sum1(nn)=-dTL1(nn,y0,-1.d0)/dTL1(n+1,y0,-1.d0)
        end if
        do k=1,m
          y=y1m(k)
          BPHI1M(k,nn+1)=TL1(nn,y0,y)+sum1(nn)*TL1(n+1,y0,y)
        end do
      end do

! Base for PHI D2 M

      do nn=0,n-1
        do k=1,M
          y=y2m(k)
          BPHI2M(k,nn+1)=chi2(nn,y0,y)
        end do
      end do

!------------------------------------------
! Base for \partial_y PHI
!------------------------------------------

! Base for DPHI D1 N (BDPHI1N)

      do nn=0,n
        do k=1,n
          y=y1n(k)
          if(nn.eq.0) then
            BDPHI1N(k,nn+1)=0.d0
          else
            BDPHI1N(k,nn+1)=dTL1(nn,y0,y)+sum1(nn)*dTL1(n+1,y0,y)
          end if
        end do
      end do

! Base for DPHI D2 N

      do nn=0,n-1
        do k=1,n
          y=y2n(k)
          BDPHI2N(k,nn+1)=dchi2(nn,y0,y)
        end do
      end do

! Base for DPHI D1 M (BDPHI1M)

      do nn=0,n
        do k=1,m
          y=y1m(k)
          if(nn.eq.0) then
            BDPHI1M(k,nn+1)=0.d0
          else
            BDPHI1M(k,nn+1)=dTL1(nn,y0,y)+sum1(nn)*dTL1(n+1,y0,y)
          end if
        end do
      end do

! Base for DPHI D2 M

      do nn=0,n-1
        do k=1,m
          y=y2m(k)
          BDPHI2M(k,nn+1)=dchi2(nn,y0,y)
        end do
      end do

!----------------------------------------------
! Matching for PHI 
!----------------------------------------------
!
      do nn=0,n

        MPHI11(nn)=TL1(nn,y0,y0)+sum1(nn)*TL1(n+1,y0,y0)

      end do

      do nn=0,n-1
 
        MPHI12(nn)=-chi2(nn,y0,y0)

      end do

!----------------------------------------------
! Super Base for PHI
!----------------------------------------------

! set to zero

      SBPHI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBPHI(k,nn+1)=BPHI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n-1
        do k=1,n
           SBPHI(k+n,nn+n+2)=BPHI2N(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,n
          SBPHI(2*n+1,nn+1)  = MPHI11(nn)
      end do

      do nn=0,n-1
          SBPHI(2*n+1,nn+n+2)= MPHI12(nn)
      end do

! Sector D1 M

      do nn=0,n
        do k=1,m
!           SBPHI(k+2*n+2,nn+1)=BPHI1M(k,nn+1)
           SBPHI(k+2*n+1,nn+1)=BPHI1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,n-1
        do k=1,m
!           SBPHI(2*n+2+m+k,nn+n+2)=BPHI2M(k,nn+1)
           SBPHI(2*n+1+m+k,nn+n+2)=BPHI2M(k,nn+1)
        end do
      end do

!----------------------------------------------
! Super Base for DPHI
!----------------------------------------------

! Set to zero

      SBDPHI=0.d0

! Sector D1 N

      do nn=0,n
        do k=1,n
           SBDPHI(k,nn+1)=BDPHI1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,n-1
        do k=1,n
           SBDPHI(k+n,nn+n+2)=BDPHI2N(k,nn+1)
        end do
      end do

! Sector D1 M

      do nn=0,n
        do k=1,m
           SBDPHI(k+2*n,nn+1)=BDPHI1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,n-1
        do k=1,m
           SBDPHI(2*n+m+k,nn+n+2)=BDPHI2M(k,nn+1)
        end do
      end do

!------------------------------------------
! Base for A
!------------------------------------------

! Base for A D1 M (BA1M) 

      do nn=0,m
        do k=1,m
          y=y1m(k)
          BA1M(k,nn+1)=psia1(nn,y0,y)
        end do
      end do

! Base for A D2 M (BA2M)

      do nn=0,m
        do k=1,m
          y=y2m(k)
          BA2M(k,nn+1)=psi2(nn,y0,y)
        end do
      end do

! Base for A D1 N (BA1N)

      do nn=0,m
        do k=1,n
          y=y1n(k)
          BA1N(k,nn+1)=psia1(nn,y0,y)
        end do
      end do

! Base for A D2 N (BA2N)

      do nn=0,m
        do k=1,n
          y=y2n(k)
          BA2N(k,nn+1)=psi2(nn,y0,y)
        end do
      end do

!------------------------------------------
! Base for \partial_y A
!------------------------------------------

! Base for DA D1 M 

      do nn=0,m
        do k=1,m
          y=y1m(k)
          BDA1M(k,nn+1)=dpsia1(nn,y0,y)
        end do
      end do

! Base for DA D2 M

      do nn=0,m
        do k=1,m
          y=y2m(k)
          BDA2M(k,nn+1)=dpsi2(nn,y0,y)
        end do
      end do

! Base for DA D1 N

      do nn=0,m
        do k=1,n
          y=y1n(k)
          BDA1N(k,nn+1)=dpsia1(nn,y0,y)
        end do
      end do

! Base for DA D2 N

      do nn=0,m
        do k=1,n
          y=y2n(k)
          BDA2N(k,nn+1)=dpsi2(nn,y0,y)
        end do
      end do

!------------------------------------------
! Matching for A and DA
!------------------------------------------
      do nn=0,m

        MA11(nn)=psia1(nn,y0,y0)
        MA12(nn)=-psi2(nn,y0,y0)
        MA21(nn)=dpsia1(nn,y0,y0)
        MA22(nn)=-dpsi2(nn,y0,y0)

      end do


!------------------------------------------
! Super Base for A
!------------------------------------------

! set to zero

      SBA=0.d0

! Sector D1 M

      do nn=0,m
        do k=1,m
           SBA(k,nn+1)=BA1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,m
        do k=1,m
           SBA(k+m,nn+m+2)=BA2M(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,m
          SBA(2*m+1,nn+1)  = MA11(nn)
          SBA(2*m+1,nn+m+2)= MA12(nn)
          SBA(2*m+2,nn+1)  = MA21(nn)
          SBA(2*m+2,nn+m+2)= MA22(nn)
      end do
! Sector D1 N

      do nn=0,m
        do k=1,n
           SBA(k+2*m+2,nn+1)=BA1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,m
        do k=1,n
           SBA(2*m+2+n+k,nn+m+2)=BA2N(k,nn+1)
        end do
      end do

!------------------------------------------
! Super Base for DA
!------------------------------------------

! Set to zero

      SBDA=0.d0

! Sector D1 M

      do nn=0,m
        do k=1,m
           SBDA(k,nn+1)=BDA1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,m
        do k=1,m
           SBDA(k+m,nn+m+2)=BDA2M(k,nn+1)
        end do
      end do

! Sector D1 N

      do nn=0,m
        do k=1,n
           SBDA(k+2*m,nn+1)=BDA1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,m
        do k=1,n
           SBDA(2*m+n+k,nn+m+2)=BDA2N(k,nn+1)
        end do
      end do

!------------------------------------------
! Base for Delta
!------------------------------------------

! Base for Delta  D1 M (BD1M) 

      do nn=0,m
        do k=1,m
          y=y1m(k)
          BD1M(k,nn+1)=psia1(nn,y0,y)
        end do
      end do

! Base for Delta D2 M (BD2M)

      do nn=0,m
        if(nn.eq.0) then
         sum2(nn)=.0d0
        else
         sum2(nn)=-dTL2(nn,y0,1.d0)/dTL2(m+1,y0,1.d0)
        end if
        do k=1,m
          y=y2m(k)
          BD2M(k,nn+1)=TL2(nn,y0,y)+sum2(nn)*TL2(m+1,y0,y)
        end do
      end do

! Base for Delta D1 N (BD1N)

      do nn=0,m
        do k=1,n
          y=y1n(k)
          BD1N(k,nn+1)=psia1(nn,y0,y)
        end do
      end do

! Base for Delta D2 N (BD2N)

      do nn=0,m
        if(nn.eq.0) then
         sum2(nn)=.0d0
        else
         sum2(nn)=-dTL2(nn,y0,1.d0)/dTL2(m+1,y0,1.d0)
        end if
        do k=1,n
          y=y2n(k)
          BD2N(k,nn+1)=TL2(nn,y0,y)+sum2(nn)*TL2(m+1,y0,y)
        end do
      end do
!------------------------------------------
! Base for \partial_y Delta (DDelta)
!------------------------------------------

! Base for DDelta D1 M (BDD1M)

      do nn=0,m
        do k=1,m
!          y=mistery1m(k)
           y=y1m(k)
          BDD1M(k,nn+1)=dpsia1(nn,y0,y)
        end do
      end do

! Base for DDelta D2 M (BDD2M)

      do nn=0,m
        if(nn.eq.0) then
         sum2(nn)=.0d0
        else
         sum2(nn)=-dTL2(nn,y0,1.d0)/dTL2(m+1,y0,1.d0)
        end if
        do k=1,m
          y=y2m(k)
          if (nn.eq.0) then
            BDD2M(k,nn+1)=0.d0
          else
            BDD2M(k,nn+1)=dTL2(nn,y0,y)+sum2(nn)*dTL2(m+1,y0,y)
          end if
        end do
      end do


! Base for DDelta D1 N (BDD1N)

      do nn=0,m
        do k=1,n
          y=y1n(k)
          BDD1N(k,nn+1)=dpsia1(nn,y0,y)
        end do
      end do

! Base for DDelta D2 N (BDD2N)

      do nn=0,m
        if(nn.eq.0) then
         sum2(nn)=.0d0
        else
         sum2(nn)=-dTL2(nn,y0,1.d0)/dTL2(m+1,y0,1.d0)
        end if
        do k=1,N
          y=y2n(k)
          if (nn.eq.0) then
            BDD2N(k,nn+1)=0.d0
          else
            BDD2N(k,nn+1)=dTL2(nn,y0,y)+sum2(nn)*dTL2(m+1,y0,y)
          end if
        end do
      end do

!----------------------------------------------
! Matching for Delta and DDelta
!----------------------------------------------
!
      do nn=0,m

        MD11(nn)=psia1(nn,y0,y0)
!        write(*,*) MD11(nn), "In frame 0"
        MD12(nn)=-(TL2(nn,y0,y0)+sum2(nn)*TL2(m+1,y0,y0))

        MD21(nn)=dpsia1(nn,y0,y0)

        if (nn.eq.0) then
          MD22(nn)=0.d0
        else
          MD22(nn)=-(dTL2(nn,y0,y0)+sum2(nn)*dTL2(m+1,y0,y0))
        end if

      end do

!----------------------------------------------
! Super Base for Delta
!----------------------------------------------
! set to zero

      SBD=0.d0

! Sector D1 M

      do nn=0,m
        do k=1,m
           SBD(k,nn+1)=BD1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,m
        do k=1,m
           SBD(k+m,nn+m+2)=BD2M(k,nn+1)
        end do
      end do

! Sector Matching at y0

      do nn=0,m
          SBD(2*m+1,nn+1)  = MD11(nn)
!          write(*,*) SBD(2*m+1,nn+1), "In frame 1"
          SBD(2*m+1,nn+m+2)= MD12(nn)
          SBD(2*m+2,nn+1)  = MD21(nn)
          SBD(2*m+2,nn+m+2)= MD22(nn)
      end do

!      do nn=0,m
!          write(*,*) SBD(2*m+1,nn+1), copySBD(2*m+1,nn+1), "In frame 2"
!      end do

! Sector D1 N

      do nn=0,m
        do k=1,n
           SBD(k+2*m+2,nn+1)=BD1N(k,nn+1)
        end do
      end do

!      do nn=0,m
!          write(*,*) SBD(2*m+1,nn+1), copySBD(2*m+1,nn+1), "In frame 3"
!      end do

! Sector D2 N

      do nn=0,m
        do k=1,n
           SBD(2*n+2+m+k,nn+m+2)=BD2N(k,nn+1)
        end do
      end do

!   This is temporal...

      copySBD=SBD

!   ...while the error is solved.

!----------------------------------------------
! Super Base for DDelta
!----------------------------------------------
! Set to zero

      SBDD=0.d0

! Sector D1 M

      do nn=0,m
        do k=1,m
           SBDD(k,nn+1)=BDD1M(k,nn+1)
        end do
      end do

! Sector D2 M

      do nn=0,m
        do k=1,m
           SBDD(k+m,nn+m+2)=BDD2M(k,nn+1)
        end do
      end do


! Sector D1 N

      do nn=0,m
        do k=1,n
           SBDD(k+2*m,nn+1)=BDD1N(k,nn+1)
        end do
      end do

! Sector D2 N

      do nn=0,m
        do k=1,n
           SBDD(2*m+n+k,nn+m+2)=BDD2N(k,nn+1)
        end do
      end do

   
!----------------------------------------------
      end subroutine frame


      double precision function T(n,y)
        implicit none
        integer n
        double precision y
        T=dcos(dble(n)*dacos(y))
      end function T

      double precision function dT(n,y)
        implicit none
        integer n
        double precision x,y
        x=dacos(y)
        if((y.eq.-1.d0).or.(y.eq.1.d0)) then
          dT=dble(n*n)*dcos(dble(n*x))/dcos(x)
        else
          dT=dble(n)*dsin(dble(n)*x)/dsin(x)
        end if
      end function dT

      double precision function iT(n,y)
        implicit none
        integer n
        double precision T,y

        if (n.eq.0) then
          iT=y
          return
        end if
        if (n.eq.1) then
          iT=0.5d0*y*y
          return
        end if
        if (n.ge.2) then
          iT=0.5d0*(T(n+1,y)/dble(n+1)+T(n-1,y)/dble(n-1))
        end if
      end function iT


      double precision function TL1(n,y0,y)
        implicit none
        double precision T,y,ji,y0
        integer n
        ji=(2.d0*y+1.d0-y0)/(1.d0+y0)
        TL1=T(n,ji)
      end function TL1
  
      double precision function dTL1(n,y0,y)
        implicit none
        double precision dT,y0,y,ji
        integer n
        ji=(2.d0*y+1.d0-y0)/(1.d0+y0)
        dTL1=(2.d0/(1.d0+y0))*dT(n,ji)
      end function dTL1

      double precision function TL2(n,y0,y)
        implicit none
        double precision T,y,ji,y0
        integer n
        ji=(2.d0*y-1.d0-y0)/(1.d0-y0)
        TL2=T(n,ji)
      end function TL2

      double precision function dTL2(n,y0,y)
        implicit none
        double precision dT,y0,y,ji
        integer n
        ji=(2.d0*y-1.d0-y0)/(1.d0-y0)
        dTL2=(2.d0/(1.d0-y0))*dT(n,ji)
      end function dTL2


      double precision function chi1(n,y0,y)
        implicit none
        integer n
        double precision TL1,y,y0
        chi1=0.5d0*(TL1(n+1,y0,y)+TL1(n,y0,y))
      end function chi1


      double precision function dchi1(n,y0,y)
        implicit none
        integer n
        double precision dTL1,y,y0
        dchi1=0.5d0*(dTL1(n+1,y0,y)+dTL1(n,y0,y))
      end function dchi1

      double precision function chi2(n,y0,y)
        implicit none
        integer n
        double precision psiaux2,y,y0,a,b
        b=1.d0
        a=-(2.d0*dble(n)+1.d0)*b/(2.d0*dble(n)+3.d0)
        chi2=a*psiaux2(n+1,y0,y)+b*psiaux2(n,y0,y)
      end function chi2

      double precision function dchi2(n,y0,y)
        implicit none
        integer n
        double precision dpsiaux2,y,y0,a,b
        b=1.d0
        a=-(2.d0*dble(n)+1.d0)*b/(2.d0*dble(n)+3.d0)
        dchi2=a*dpsiaux2(n+1,y0,y)+b*dpsiaux2(n,y0,y)
      end function dchi2


      double precision function psia1(n,y0,y)
        implicit none
        integer n
        double precision y,y0,chi1,a,b
        b=0.5d0
        a=(2.d0*dble(n)+1.d0)*b/(2.d0*dble(n)+3.d0)
        psia1=a*chi1(n+1,y0,y)+b*chi1(n,y0,y)
      end function psia1

      double precision function dpsia1(n,y0,y)
        implicit none
        integer n
        double precision y,y0,dchi1,a,b
        b=0.5d0
        a=(2.d0*dble(n)+1.d0)*b/(2.d0*dble(n)+3.d0)
        dpsia1=a*dchi1(n+1,y0,y)+b*dchi1(n,y0,y)
      end function dpsia1

      double precision function psiaux2(n,y0,y)
        implicit none
        integer n
        double precision TL2,y0,y
        psiaux2=0.5d0*(TL2(n+1,y0,y)-TL2(n,y0,y))
      end function psiaux2

      double precision function dpsiaux2(n,y0,y)
        implicit none
        integer n
        double precision dTL2,y0,y
        dpsiaux2=0.5d0*(dTL2(n+1,y0,y)-dTL2(n,y0,y))
      end function dpsiaux2

      double precision function psi2(n,y0,y)
        implicit none
        integer n
        double precision chi2,y0,a,b,y
        b=0.5d0
        a=-(dble(n)+1.d0)*(2.d0*dble(n)+1.d0)*b/
     .     ( (dble(n)+2.d0)*(2.d0*dble(n)+3.d0) )
        psi2=a*chi2(n+1,y0,y) + b*chi2(n,y0,y)
      end function psi2
 
      double precision function dpsi2(n,y0,y)
        implicit none
        integer n
        double precision dchi2,y0,a,b,y
        b=0.5d0
        a=-(dble(n)+1.d0)*(2.d0*dble(n)+1.d0)*b/
     .     ( (dble(n)+2.d0)*(2.d0*dble(n)+3.d0) )
        dpsi2=a*dchi2(n+1,y0,y) + b*dchi2(n,y0,y)
      end function dpsi2

! La siguiente función está definida en el 
! script de Maple pero no se usa...

      double precision function psidelta2(n,y0,y)
        implicit none
        integer n
        double precision a,b,c,TL2,y0,y
        c=0.25d0
        b=0.5d0
        a=-(c*dble(n)**2+b*(dble(n)+1.d0)**2)/((dble(n)+2.d0)**2)
        psidelta2=a*TL2(n+2,y0,y)+b*TL2(n+1,y0,y)+c*TL2(n,y0,y)
      end function psidelta2
