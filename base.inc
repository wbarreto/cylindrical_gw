      integer NMAX
      parameter(NMAX=300)
      double precision SBPSI(NMAX,NMAX),SBDPSI(NMAX,NMAX)
      double precision SBOMEGA(NMAX,NMAX),SBDOMEGA(NMAX,NMAX)
      double precision SGRID(NMAX)
      common/fields1/SBPSI,SBDPSI
      common/fields2/SBOMEGA,SBDOMEGA
      common/grid/SGRID
