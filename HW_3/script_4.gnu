set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data in "es_2.dat"
#-------------------------------------------------


set title 'Error distribution for different error esteem methods' font "Verdana,15"
#set title 'Esteem distribution for MC methods' font "Verdana,15"
set key font "Verdana,11"


set xlabel 'ln(n)' font "Verdana,13"
set xtics font "Verdana,13"
set ylabel 'ln(Error)' font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[5000:100100]
#set yrange[0:1]
set pointsize 1.3


#---------------------------------
#Error Distribution 1
#---------------------------------
#f(x) =  a*x + b
#g(x) =  c*x + d

#fit f(x) 'plot.dat' u (log($1)):(log($3)) via a,b
#fit g(x) 'plot.dat' u (log($1)):(log($5)) via c,d

#p 'plot.dat' u (log($1)):(log($3)) w p lc "gold" pt 5 title '{/Symbol s}_n'
#replot 'plot.dat'u (log($1)):(log($5)) w p lc "red" pt 7 title '{/Symbol D}_n'

#replot f(x) w l lc "gold" notitle
#replot g(x) w l lc "red" notitle



#---------------------------------
#Error Distribution 2
#---------------------------------
f(x) =  a*x + b
g(x) =  c*x + d
h(x) =  e*x + f
i(x) =  g*x + h
j(x) =  i*x + j
k(x) =  k*x + l


fit f(x) 'final.dat' i 0 u (log($1)):(log($3)) via a,b
fit g(x) 'final.dat' i 0 u (log($1)):(log($4)) via c,d
fit h(x) 'final.dat' i 0 u (log($1)):(log($5)) via e,f
fit i(x) 'final.dat' i 1 u (log($1)):(log($3)) via g,h
fit j(x) 'final.dat' i 1 u (log($1)):(log($4)) via i,j
fit k(x) 'final.dat' i 1 u (log($1)):(log($5)) via k,l



p 'final.dat' i 0 u (log($1)):(log($3)) w p lc "gold" pt 5 title '{/Symbol s}_m'
replot 'final.dat' i 0 u (log($1)):(log($5)) w p lc "blue" pt 3 title '{/Symbol s}_n/√n AVG OF AVGS'
replot 'final.dat' i 1 u (log($1)):(log($3)) w p lc "green" pt 13 title '{/Symbol s}_s/√s'
replot 'final.dat' i 1 u (log($1)):(log($5)) w p lc "purple" pt 9 title '{/Symbol s}_n/√n BLOCK AVG'
replot 'final.dat' i 1 u (log($1)):(log($4)) w p lc "black" pt 15 title '{/Symbol D}_n BLOCK AVG'

replot f(x) w l lc "gold" notitle
#replot g(x) w l lc "red" notitle
replot h(x) w l lc "blue" notitle
replot i(x) w l lc "green" notitle
replot j(x) w l lc "black" notitle
replot k(x) w l lc "purple" notitle
