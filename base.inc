      integer NMAX,ene,eme
      parameter(NMAX=300)
      parameter(ene=2*(40+1))
      parameter(eme=2*(52+1))
      double precision SBPI(NMAX,NMAX),SBDPI(NMAX,NMAX)
      double precision SBPHI(NMAX,NMAX),SBDPHI(NMAX,NMAX)
      double precision SBA(NMAX,NMAX),SBDA(NMAX,NMAX)
      double precision SBD(NMAX,NMAX),SBDD(NMAX,NMAX)
      double precision SGRID(NMAX)
!      double precision SBPI(ene,ene),SBDPI(ene,ene)
!      double precision SBPHI(ene,ene-1),SBDPHI(ene,ene-1)
!      double precision SBA(eme,eme),SBDA(eme,eme)
!      double precision SBD(eme,eme),SBDD(eme,eme)
!      double precision SGRID(ene+eme-2)
      common/fields1/SBPI,SBDPI
      common/fields2/SBPHI,SBDPHI
      common/fields3/SBA,SBDA
      common/fields4/SBD,SBDD
      common/grid/SGRID
