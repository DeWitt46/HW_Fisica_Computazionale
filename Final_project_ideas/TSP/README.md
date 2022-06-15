# SETUP
- Open GA_for_TSP.f90 code to check internal algorithm parameters; <br>
- Modify RND_PLACE_FLAG according to the configuration you want to run: .TRUE. for random generation, .FALSE. for an existing one; <br>
- Compile GA_for_TSP.f90 code with some FORTRAN compiler (gfortran worked for me); <br>
- Run "a.exe".

NOTE: Configurations used in the report are in "CONFIGURATION" folder. If you want to run that, just type "a.exe < configuration_file" <br>
after checked that RND_PLACE_FLAG is equal to .FALSE.
