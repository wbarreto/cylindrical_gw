      integer NMAX
      parameter(NMAX=300)
      double precision SBPI(NMAX,NMAX),SBDPI(NMAX,NMAX)
      double precision SBPHI(NMAX,NMAX),SBDPHI(NMAX,NMAX)
      double precision SBA(NMAX,NMAX),SBDA(NMAX,NMAX)
      double precision SBD(NMAX,NMAX),SBDD(NMAX,NMAX)
      double precision SGRID(NMAX)
      common/fields1/SBPI,SBDPI
      common/fields2/SBPHI,SBDPHI
      common/fields3/SBA,SBDA
      common/fields4/SBD,SBDD
      common/grid/SGRID
