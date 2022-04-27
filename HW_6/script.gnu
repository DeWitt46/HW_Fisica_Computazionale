set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data of es_1
#-------------------------------------------------


set key font "Verdana,11"



set xtics font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[-0.2:0.2]
#set yrange[0:1]
set pointsize 1.1
#set lw 1.3



#----------------------------------------------------------------------------------------------
#ISOTROPIC ENERGY
#EXP:100000   7.0000000000000000       0.49999999999999989       0.50000000000000000        0.0000000000000000       0.32124000000000003
#PAR:100000   7.0000000000000000        2.0533333460489906       0.59104909734616917       0.46106521369143294       0.36564000000000002     
#----------------------------------------------------------------------------------------------


#set title 'Harmonic oscillator energy with Variational Monte Carlo' font "Verdana,15"
#set ylabel 'Energy' font "Verdana,13"


#set xlabel '{/Symbol b} parameter' font "Verdana,13"
#set yrange[-0.3:2.1]
#set xrange[0.1:0.9]


#a     = 2.0533333460489906
#b     = 0.49999999999999989
#E_par = 0.59104909734616917
#E_exp =  0.50000000000000000 
#E_th = 0.5
#f(x) = E_exp
#g(x) = E_par
#h(x) = E_th


#p 'Metropolis_exp.dat' u 3:4 w p pt 5 lc "red" title 'VMC with exponential eigenstate'
#replot 'Metropolis_exp.dat' u 3:4:5 w yerrorbars pt 13 lc "red" notitle
##set arrow from b,-1 to b,3 nohead lc "red"
#replot f(x) lw 2 lc "red" title 'Minimum'
#replot h(x) lc "black" title 'Referement Value'


#set xlabel 'a parameter' font "Verdana,13"
#set yrange[-0.2:1.85]
#set xrange[1.:3.1]

#p 'Metropolis_par.dat' u 3:4 w p pt 5 lc "blue" title 'VMC with parabolic eigenstate'
#replot 'Metropolis_par.dat' u 3:4:5 w yerrorbars pt 13 lc "blue" notitle
##set arrow from a,-0.5 to a,2. nohead lc "blue"
#replot g(x) lw 2 lc "blue" title 'Minimum'
#replot h(x) lw 2 lc "black" title 'Referement Value'


#POTENTIAL PLOT
#set title 'Harmonic oscillator eigenstates comparison' font "Verdana,15"
#set ylabel '' font "Verdana,13"
#set xlabel 'x' font "Verdana,13"
#set xrange[-a:a]
#set yrange[0.:1.25]

#p 'pot.dat' u 1:2 w l  lw 2 lc "gold" title '0.5x^2'
#replot 'pot.dat' u 1:3 w l  lw 2 lc "red" title 'Exponential eigenstate ({/Symbol b} = 0.5)'
#replot f(x) lc "red" lw 2 title 'E_{exp}'
#replot 'pot.dat' u 1:4 w l  lw 2 lc "blue" title 'Parabolic eigenstate'
#replot g(x) lc "blue" lw 2 title 'E_{par}'
#replot 'pot.dat' u 1:3 w l lc "black" title 'Theoretical eigenstate'
#replot h(x) lc "black" title 'E_{th}'






#----------------------------------
#ANISOTROPIC ENERGY
#EXP:100000   6.0000000000000000       0.62426667372385658       0.57059759044518532        8.8407887996059156E-002  0.33446999999999999
#PAR:100000   6.0000000000000000        1.8166666666666667       0.67317738082389866       0.53596668453098051       0.37498999999999999
#----------------------------------

set title 'Anharmonic oscillator energy with Variational Monte Carlo' font "Verdana,15"
set ylabel 'Energy' font "Verdana,13"
#set xlabel '{/Symbol b} parameter' font "Verdana,13"
#set yrange[-0.4:1.9]
#set xrange[0.25:1.]


a     = 1.8166666666666667
b     = 0.62426667372385658
E_par = 0.67317738082389866
E_exp = 0.57059759044518532
E_pert = 0.5937
f(x) = E_exp
g(x) = E_par
h(x) = E_pert


#p 'Metropolis_aniso_exp.dat' u 3:4 w p pt 5 lc "red" title 'VMC with exponential eigenstate'
#replot 'Metropolis_aniso_exp.dat' u 3:4:5 w yerrorbars pt 13 lc "red" notitle
##set arrow from b,-1 to b,3 nohead lc "red"
#replot f(x) lw 2 lc "red" title 'Minimum'
#replot h(x) lw 2 lc "black" title 'Perturbation theory'


set xlabel 'a parameter' font "Verdana,13"
set yrange[-0.4:2.3]
set xrange[1.017:2.617]

p 'Metropolis_aniso_par.dat' u 3:4 w p pt 5 lc "blue" title 'VMC with parabolic eigenstate'
replot 'Metropolis_aniso_par.dat' u 3:4:5 w yerrorbars pt 13 lc "blue" notitle
##set arrow from a,-0.8 to a,3.5 nohead lc "blue"
replot g(x) lw 2 lc "blue" title 'Minimum'
replot h(x) lw 2 lc "black" title 'Perturbation theory'


#POTENTIAL PLOT
#set title 'Anharmonic oscillator eigenstates comparison' font "Verdana,15"
#set ylabel '' font "Verdana,13"
#set xlabel 'x' font "Verdana,13"
#set xrange[-a:a]
#set yrange[0.:1.25]

#p 'Aniso_pot.dat' u 1:2 w l  lw 2 lc "gold" title '0.5x^2 + 0.125x^4'
#replot 'Aniso_pot.dat' u 1:3 w l  lw 2 lc "red" title 'Exponential eigenstate ({/Symbol b} = 0.62)'
#replot f(x) lc "red" lw 2 title 'E_{exp}'
#replot 'Aniso_pot.dat' u 1:4 w l  lw 2 lc "blue" title 'Parabolic eigenstate'
#replot g(x) lc "blue" lw 2 title 'E_{par}'
#replot 'pot.dat' u 1:3 w l  lw 2 lc "black" title 'Perturbation eigenstate ({/Symbol b} = 0.5)'
#replot h(x) lc "black" lw 2 title 'E_{perth}'






