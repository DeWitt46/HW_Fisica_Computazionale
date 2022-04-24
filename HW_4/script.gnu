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


#----------------------------------
#Uniform distribution
#----------------------------------
#set title 'Frequencies histogram for uniform distribution' font "Verdana,15"
#set xrange[-0.1:0.1]
#set ylabel 'counts (normalized at 1)' font "Verdana,13"
#set xlabel 'x_N' font "Verdana,13"

#p 'hist.dat' i 2 u 1:2 w boxes lc "gold" title 'N = 850'
#replot 'hist.dat' i 7 u 1:2 w boxes lc "red" title 'N = 2000'
#replot 'hist.dat' i 7 u 1:3 w boxes lc "black" title 'GAUSSIAN REFEREMENT'


#set title 'CLT results for uniform distribution' font "Verdana,15"
#set ylabel 'Relative deviations (percentage)' font "Verdana,13"
#set xlabel 'N' font "Verdana,13"
#set xrange[450:2050]
#f(x) =  0

#p 'es_2.dat' i 0 u 1:(($2-$3)*100) w lp lc "gold" pt 5 title '{/Symbol m} - {/Symbol m}_x'
#replot 'es_2.dat' i 0 u 1:((($4-$5)/$4)*100) w lp lc "red" pt 13 title '({/Symbol s} - {/Symbol s}_x)/{/Symbol s}'
#replot 'es_2.dat' i 0 u 1:((($6-$7)/$7)*100) w lp lc "green" pt 15 title '(<z^4> - 3<z^2>^2)/<z^4>'
#replot f(x) w l lc "black" notitle




#----------------------------------
#Exponential distribution
#----------------------------------
#set title 'Frequencies histogram for exponential distribution' font "Verdana,15"
#set ylabel 'counts (normalized at 1)' font "Verdana,13"
#set xlabel 'x_N' font "Verdana,13"
#set xrange[0.8:1.2]

#p 'hist.dat' i 8 u 1:2 w boxes lc "gold" title 'N = 500'
#replot 'hist.dat' i 14 u 1:2 w boxes lc "red" title 'N = 2000'
#replot 'hist.dat' i 14 u 1:3 w boxes lc "black" title 'GAUSSIAN REFEREMENT'


#set title 'CLT results for exponential distribution' font "Verdana,15"
#set ylabel 'Relative deviations (percentage)' font "Verdana,13"
#set xlabel 'N' font "Verdana,13"
#set xrange[450:2050]

#f(x) =  0

#p 'es_2.dat' i 0 u 1:(($2-$3)*100) w lp lc "gold" pt 5 title '{/Symbol m} - {/Symbol m}_x'
#replot 'es_2.dat' i 0 u 1:((($4-$5)/$5)*100) w lp lc "red" pt 13 title '({/Symbol s} - {/Symbol s}_x)/{/Symbol s}'
#replot 'es_2.dat' i 0 u 1:((($6-$7)/$7)*100) w lp lc "green" pt 15 title '(<z^4> - 3<z^2>^2)/<z^4>'
#replot f(x) w l lc "black" notitle



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


set title 'CLT results for Lorentz distribution' font "Verdana,15"
set ylabel 'Relative deviations (percentage)' font "Verdana,13"
set xlabel 'N' font "Verdana,13"
set xrange[450:2050]

f(x) =  0

p 'es_2.dat' i 2 u 1:(($2-$3)*100) w lp lc "gold" pt 5 title '{/Symbol m}_x - {/Symbol m}'
replot 'es_2.dat' i 2 u 1:((($4-$5)/$5)*100) w lp lc "red" pt 13 title '({/Symbol s}_x - {/Symbol s})/{/Symbol s}'
#replot 'es_2.dat' i 2 u 1:((($6-$7)/$7)*100) w lp lc "green" pt 15 title '(<z^4> - 3<z^2>^2)/<z^4>'
replot f(x) w l lc "black" notitle

