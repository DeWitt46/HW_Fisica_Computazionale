set termoption enhanced

#set title 'z histogram fit for different {/Symbol l} values' font "Verdana,15"
#set title 'z histogram for different {/Symbol l} values' font "Verdana,15"
#set key left font "Verdana,12"
set key right font "Verdana,12"


set autoscale
set pointsize 0.3
set xtics font "Verdana,11"
#set xrange[0:10100]
set ytics font "Verdana,11"




#Plot x_i e x^2_i" di pochi runner
#set xlabel 't (in step)' font "Verdana,11"
#set title 'Walks for different runners' font "Verdana,15"

#set ylabel 'x (in units of l)' font "Verdana,11"
#set yrange[-250:200]
#p 'es_1_medio.dat' i 0 u 2:3 w l lc "gold" title 'runner 1'
#replot 'es_1_medio.dat' i 1 u 2:3 w l lc "green" title 'runner 2'
#replot 'es_1_medio.dat' i 2 u 2:3 w l lc "black" title 'runner 3'
#replot 'es_1_medio.dat' i 4 u 2:5 w l lc "red" title 'avg. 10 runner'
#f(x) = 0
#replot f(x) w l lc "purple" title 'expected'

#set ylabel 'x^2 (in units of l^2)' font "Verdana,11"
#set yrange[0:27500]
#p 'es_1_medio.dat' i 0 u 2:4 w l lc "gold" title 'runner 1'
#replot 'es_1_medio.dat' i 1 u 2:4 w l lc "green" title 'runner 2'
#replot 'es_1_medio.dat' i 2 u 2:4 w l lc "black" title 'runner 3'
#replot 'es_1_medio.dat' i 4 u 2:6 w l lc "red" title 'avg. 10 runner'
#f(x) = x
#replot f(x) w l lc "purple" title 'expected'






#Plot valori medi al variare di N
#set title 'Mean square displacement distribution fit' font "Verdana,15"
#set ylabel 'ln({{/Symbol s}_N}^2)' font "Verdana,13"
#set xlabel 'ln(N)' font "Verdana,13"
#set xrange[0:10500]
#set yrange[-8:4]

#f(x) = a*x + b
#fit f(x) 'es_1_medio.dat' i 3 u (2*log($2)):(log($7)) via a,b
#p 'es_1_medio.dat' i 3 u (2*log($2)):(log($7)) w p lc "black" pt 5 title 'numerical data'
#replot f(x) w l lw 2 lc "purple" title 'fit'




#p 'es_1_medio.dat' i 0 u 2:5 w l lc "gold" pt 5 title 'x_i medio 250 run'
#replot 'es_1_medio.dat' i 1 u 2:5 w p lc pt 7 "green" title 'x_i medio 500 run'
#replot 'es_1_medio.dat' i 2 u 2:5 w p lc pt 13 "black" title 'x_i medio 750 run'
#replot 'es_1_medio.dat' i 3 u 2:5 w p lc pt 15 "red" title 'x_i medio 1000 run'

#p 'es_1_medio.dat' i 0 u 2:6 w l lc "gold" pt 5 title 'x^2_i medio 250 run'
#replot 'es_1_medio.dat' i 1 u 2:6 w p lc "green" pt 7 title 'x^2_i medio 500 run'
#replot 'es_1_medio.dat' i 2 u 2:6 w p lc "black" pt 13 title 'x^2_i medio 750 run'
#replot 'es_1_medio.dat' i 3 u 2:6 w p lc "red" pt 15 title 'x^2_i medio 1000 run'

#p 'es_1_medio.dat' i 0 u 2:7 w p lc "gold" pt 5 title '{\Symbol s}^2 250 run'
#replot 'es_1_medio.dat' i 1 u 2:7 w p lc "green" pt 7 title '{\Symbol s}^2 500 run'
#replot 'es_1_medio.dat' i 2 u 2:7 w p lc "black" pt 13 title '{\Symbol s}^2 750 run'
#replot 'es_1_medio.dat' i 3 u 2:7 w p lc "red" pt 15 title '{\Symbol s}^2 1000 run'












#Plot valori medi al variare numero runner
#set title 'Error convergence' font "Verdana,15"
#set xlabel 'Number of runners' font "Verdana,11"
#set ylabel '{/Symbol D} (in percentage)' font "Verdana,11"
#set xrange[20:1000]

#p 'es_1.dat' u 1:2 w p lc "green" pt 5 title 'x_N medio'
#replot 'es_1.dat' u 1:3 w p lc "black" pt 7 title 'x^2_N medio'
#replot 'es_1.dat' u 1:4 w p lc "black" pt 13 title '{/Symbol s}^2'
#p 'es_1.dat' u 1:(100*$5) w p lc "gold" pt 15 title '{/Symbol D} from simulation'
#f(x) = 5
#replot f(x) w l lc "black" title 'Tollerance'




#Plot istogramma
set title 'Frequencies histogram and gaussian distribution comparison' font "Verdana,15"
set ylabel 'counts' font "Verdana,13"
set xlabel 'x_N (in units of l)' font "Verdana,13"
set xrange[-30:30]


p 'hist_1.dat' u 1:2 w boxes lw 2 lc "gold" title 'Frequencies N = 8'
#p 'hist_1.dat' u 1:3 w boxes  lw 2 lc "green" title 'Frequencies N = 16'
replot 'hist_1.dat' u 1:4 w boxes  lw 2 lc "black" title 'Frequencies N = 32'
#replot 'hist_1.dat' u 1:5 w boxes  lw 2 lc "red" title 'Frequencies N = 64'
replot 'hist_1.dat' u 1:6 w l lw 2 lc "purple" title "Gaussian esteem N = 8"
#replot 'hist_1.dat' u 1:7 w l lw 2 lc "purple" title "Gaussian esteem N = 16"
replot 'hist_1.dat' u 1:8 w l lw 2 lc "blue" title "Gaussian esteem N = 32"
#replot 'hist_1.dat' u 1:9 w l lw 2 lc "blue" title "Gaussian esteem N = 64"

