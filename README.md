# cylindrical_gw

Developer(s): W. Barreto
Initial date (of developing):      29.06.18
Final date (ready for production): 13.07.18 

CHANGELOG:

Evolving (compiling well and living code):

- makefile
- control deck
- driver (cgw.f)
- grid, basis and frame as an unchanged strutcture (base.f)
- maximal dimensions and allocations (base.inc)
- Solver for a linear system of equations (gaussj.f) (from Numerical Recipes)
- Initial data (initial.f)
- Dynamical system (dynsys.f) 
- Ordinary differential equations solver (rk.4)
- Metric functions and its spatial derivatives (metric.f)
- Input (initial spectral modes or integrated ones) Output (metric functions evolution;
  Bondi mass and Weyl scalar). (io.f) (write_files) (here is the calculation post-integration
  at each time-step); (write_asymp).  
NOTES:

Check only one of these stages:

- still developing &  debugging 
- still not used in production (check)
- in production 

TO DO (next tasks):

- Adapt ko.f (standard Kreiss-Oliger noise filter, if necessary), 
- Write the pseudo code;
- Algorithmic;
- Code's history, referential Maple script and pre-codes;
- F90, C and Python versions;
 
