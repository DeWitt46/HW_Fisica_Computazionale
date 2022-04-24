set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data in "es_2.dat"
#-------------------------------------------------


set key font "Verdana,11"



set xtics font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[-0.2:0.2]
#set yrange[0:1]
set pointsize 1.3
#set lw 1.3

#----------------------------------
#HISTOGRAMS
#----------------------------------
#set title 'Frequencies histogram from Metropolis algorithm' font "Verdana,15"
#set yrange[0:0.45]
#set ylabel 'counts (normalized at 1)' font "Verdana,13"
#set xlabel 'x' font "Verdana,13"


#p 'hist_metropolis.dat' i 0 u 1:2 w boxes lw 2 lc "gold" title 'N = 100'
#replot 'hist_metropolis.dat' i 1 u 1:2 w boxes lw 2 lc "red" title 'N = 1000'
#replot 'hist_metropolis.dat' i 2 u 1:2 w boxes lw 2 lc "green" title 'N = 10000'
#replot 'hist_metropolis.dat' i 3 u 1:2 w boxes lw 2 lc "blue" title 'N = 100000'
#replot 'hist_metropolis.dat' i 3 u 1:3 w boxes lw 2 lc "black" title 'GAUSSIAN REFEREMENT'


#----------------------------------
#RESULTS
#----------------------------------
#set title 'Acceptance Ratio distribution for fixed N = 100000' font "Verdana,15"
#set ylabel 'Acceptance Ratio' font "Verdana,13"
#set xlabel 'log_{10}({/Symbol d}/{/Symbol s})' font "Verdana,13"
#set xrange[0.5:1.2]

#f(x) =  a*x + b
#fit[0.5:1.2] f(x) 'gauss_metropolis.dat' u (log10($3)):4 via a,b

#p 'gauss_metropolis.dat' u (log10($3)):4 lc "gold" pt 13 notitle
#replot f(x) w l lc "black" notitle







#set title 'Numerical esteem and analytical first moments comparison' font "Verdana,15"
#set ylabel 'Relative deviations (in percentage)' font "Verdana,13"
#set xlabel 'N' font "Verdana,13"
#set xrange[0:30500]
#set xtics 0,5000,30000

#f(x) = 5
#p 'gauss_metropolis.dat' u 1:(abs($5-$6)*100) w lp lc "gold" pt 5 title '|<x> - <x>_{th}|'
#replot 'gauss_metropolis.dat' u 1:(abs(($7-$8)/$8)*100) w lp lc "red" pt 13 title '|(<x^2> - <x^2>_{th})/<x^2>_{th}|'
#plot 'gauss_metropolis.dat' u 1:(abs(($9-$10)/$10)*100) w lp lc "gold" pt 15 title '|{/Symbol s}^2 -\
#{/Symbol s}^2_{th}/{/Symbol s}^2_{th}|'
#replot 'gauss_metropolis.dat' u 1:(abs($11-$12)*100) w lp lc "red" pt 15 title '|(<x^3> - <x^3>_{th})|'
#replot 'gauss_metropolis.dat' u 1:(abs(($13-$14)/$14)*100) w lp lc "purple" pt 15 title '|<x^4> - <x^4>_{th}/<x^4>_{th}|'
#replot f(x) lc 'black' title "tollerance"







#set title 'Auto-Correlation function for different {/Symbol d}/{/Symbol s} fractions' font "Verdana,15"
#set ylabel 'C(j)' font "Verdana,13"
#set xlabel 'j' font "Verdana,13"
#set yrange[-0.1:1.1]
#set xtics 0,5000,30000

#f(x) = 0
#p 'corr_metropolis.dat' i 0 u 4:5 lc "gold" pt 5 title '{/Symbol d}/{/Symbol s} = 2'
#replot 'corr_metropolis.dat' i 1 u 4:5 lc "red" pt 13 title '{/Symbol d}/{/Symbol s} = 5'
#replot 'corr_metropolis.dat' i 2 u 4:5 lc "green" pt 15 title '{/Symbol d}/{/Symbol s} = 8'
#replot f(x) lc 'black' notitle

        

set title 'Auto-Correlation function comparison' font "Verdana,15"
set ylabel 'C(j)' font "Verdana,13"
set xlabel 'j' font "Verdana,13"

f(x) = 0
p 'corr_metropolis.dat' i 2 u 4:5 lc "gold" pt 5 title 'Metropolis, {/Symbol d}/{/Symbol s} = 8'
replot 'corr_metropolis.dat' i 2 u 4:6 lc "red" pt 13 title 'Box-Muller'
replot f(x) lc 'black' notitle


#----------------------------------
#Lorentz distribution
#----------------------------------
#set title 'Frequencies histogram for Lorentz distribution' font "Verdana,15"
#set ylabel 'counts (normalized at 1)' font "Verdana,13"
#set xlabel 'x_N' font "Verdana,13"
#set yrange[0:50]
#set xrange[-10:10]

#f(x) = 1/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
#g(x) = y0/(pi*((x0 - x)**2 + y0**2))
#fit f(x) 'hist.dat' i 23 u 1:4 via sigma, mu
#x0 = -0.02
#y0 = 0.08
#fit g(x) 'hist.dat' i 23 u 1:2 via x0, y0

#p 'hist.dat' i 16 u 1:2 w boxes lc "gold" title 'N = 500'
#replot 'hist.dat' i 23 u 1:2 w boxes lc "red" title 'N = 2000'
#replot 'hist.dat' i 23 u 1:3 w boxes lc "black" title 'GAUSSIAN REFEREMENT'
#replot f(x) lc "black" title 'GAUSSIAN FIT'
#replot g(x) lc "black" title 'CAUCHY FIT'

#set xlabel 'r_i' font "Verdana,13"
#set xrange[-10:10]
#p 'hist.dat' i 23 u 1:4 w boxes lc "red" title 'm=2000'
#replot 'hist.dat' i 23 u 1:5 w l lc "black" title 'LORENTZ REFEREMENT'
#replot f(x) lc "blue" title 'GAUSSIAN FIT'


#set title 'CLT results for Lorentz distribution' font "Verdana,15"
#set ylabel 'Relative deviations (percentage)' font "Verdana,13"
#set xlabel 'N' font "Verdana,13"
#set xrange[450:2050]

#f(x) =  0

#p 'es_2.dat' i 2 u 1:(($2-$3)*100) w lp lc "gold" pt 5 title '{/Symbol m}_x - {/Symbol m}'
#replot 'es_2.dat' i 2 u 1:((($4-$5)/$5)*100) w lp lc "red" pt 13 title '({/Symbol s}_x - {/Symbol s})/{/Symbol s}'
#replot 'es_2.dat' i 2 u 1:((($6-$7)/$7)*100) w lp lc "green" pt 15 title '(<z^4> - 3<z^2>^2)/<z^4>'
#replot f(x) w l lc "black" notitle

