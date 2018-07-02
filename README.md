# cylindrical_gw

Developer(s): W. Barreto
Initial date (of developing):      29.06.18
Final date (ready for production): 

CHANGELOG:

Evolving (compiling well and living code):

- makefile
- control deck
- driver (cgw.f)
- grid, basis and frame as an unchanged strutcture (base.f)
- maximal dimensions and allocations (base.inc)
- Solver for a linear system of equations (gaussj.f) (from Numerical Recipes)
- Initial data (initial.f)

NOTES:

Check only one of these stages:

- still developing (check) (50%)
- still debugging
- still not used in production
- in production 

TO DO (next tasks):

- Adapt metric.f, dynsys.f, rk.f, ko.f (standard Kreiss-Oliger noise filter,
  if necessary), io.f; 
- write the pseudo code;
- Algorithmic;
- Code's history, referential Maple script and pre-codes;
- F90, C and Python versions;
 
