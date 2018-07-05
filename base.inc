      integer NMAX
      parameter(NMAX=300)
      double precision SBPSI(NMAX,NMAX),SBDPSI(NMAX,NMAX),
     .                 SBDDPSI(NMAX,NMAX)
      double precision SBOMEGA(NMAX,NMAX),SBDOMEGA(NMAX,NMAX),
     .                 SBDDOMEGA(NMAX,NMAX)
      double precision SGRID(NMAX)
      common/fields1/SBPSI,SBDPSI,SBDDPSI
      common/fields2/SBOMEGA,SBDOMEGA,SBDDOMEGA
      common/grid/SGRID
