set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data in "es_2.dat"
#-------------------------------------------------


set title 'Error distribution for MC methods' font "Verdana,15"
#set title 'Esteem distribution for MC methods' font "Verdana,15"
set key font "Verdana,11"


set xlabel 'ln(n)' font "Verdana,13"
set xtics font "Verdana,13"
set ylabel 'ln(Error)' font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[5000:100100]
#set yrange[0:1]
set pointsize 1


#----------------------------------
#Esteem Distribution
#----------------------------------
#f(x) = (sqrt(pi)/2)*erf(1.)
#p 'es_2.dat' i 0 u 1:2 w lp lc "gold" pt 5 title 'SAMPLE MEAN'
#replot 'es_2.dat' i 1 u 1:2 w lp lc "red" pt 13 title 'IMPORTANCE SAMPLING'
#replot f(x) lc "black" pt 13 title 'REFEREMENT VALUE'


#---------------------------------
#Error Distribution
#---------------------------------
f(x) =  a*x + b
g(x) =  c*x + d
f1(x) =  e*x + f
g1(x) =  g*x + h

fit f(x) 'plot.dat' i 0 u (log($1)):(log($4)) via a,b
fit g(x) 'plot.dat' i 1 u (log($1)):(log($4)) via c,d
fit f1(x) 'plot.dat' i 0 u (log($1)):(log($5)) via e,f
fit g1(x) 'plot.dat' i 1 u (log($1)):(log($5)) via g,h

p 'plot.dat' i 0 u (log($1)):(log($4)) w p lc "gold" pt 5 title '{/Symbol s}_n**0.5 SAMPLE MEAN'
replot 'plot.dat' i 0 u (log($1)):(log($5)) w p lc "red" pt 7 title '{/Symbol D}_n SAMPLE MEAN'
replot 'plot.dat' i 1 u (log($1)):(log($4)) w p lc "green" pt 13 title '{/Symbol s}_n**0.5 IMPORTANT SAMPLING'
replot 'plot.dat' i 1 u (log($1)):(log($5)) w p lc "black" pt 15 title '{/Symbol D}_n IMPORTANT SAMPLING'

replot f(x) w l lc "gold" notitle
replot f1(x) w l lc "red" notitle
replot g(x) w l lc "green" notitle
replot g1(x) w l lc "black" notitle