set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data of es_1
#-------------------------------------------------


set key t r font "Verdana,11"



set xtics font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[-0.2:0.2]
#set yrange[0:1]
set pointsize 1.
#set lw 1.3


#----------------------------------------------------------------------------------------------
#PLOT VARYING MC STEP
#----------------------------------------------------------------------------------------------
set ylabel ''
set xlabel 'MC step' font "Verdana,13"


#PBC, T = 2

#set title 'Istantaneous and cumulative variables of interest at T= 2J/k_{B} for PBC' font "Verdana,15"

#ENERGY


#set yrange[-2:1]
#set xrange[0:1000]


#p 'ist_tot.dat' i 0 u 1:2 w lp pt 5 lc "red" title 'Istantaneous E/N (L=4)'
#replot 'ist_cum_tot.dat' i 0 u 1:2 w l lw 2 lc "red" title 'Cumulative E/N (L=4)' 
#replot 'ist_tot.dat' i 2 u 1:2 w lp pt 13 lc "blue" title 'Istantaneous E/N (L=30)' 
#replot 'ist_cum_tot.dat' i 2 u 1:2 w l lw 2 lc "blue" title 'Cumulative E/N (L=30)' 

#MAGNETIZATION PBC

#set ylabel 'M/N' font "Verdana,13"
#set yrange[-0.2:1.85]
#set xrange[1.:3.1]

#replot 'ist_tot.dat' i 0 u 1:3 w p pt 5 lc "gold" title 'Istantaneous M/N (L=4)'
#replot 'ist_cum_tot.dat' i 0 u 1:3 w l lw 2 lc "gold" title 'Cumulative M/N (L=4)' 
#replot 'ist_tot.dat' i 2 u 1:3 w p pt 13 lc "green" title 'Istantaneous M/N (L=30)'
#replot 'ist_cum_tot.dat' i 2 u 1:3 w l lw 2 lc "green" title 'Cumulative M/N (L=30)' 


#ABS(MAGNETIZATION)

#set ylabel '|M|/N' font "Verdana,13"
#set yrange[-0.2:1.85]
#set xrange[1.:3.1]

#replot 'ist_tot.dat' i 0 u 1:4 w lp pt 5 lc "gold" title 'Istantaneous |M|/N (L=4)'
#replot 'ist_cum_tot.dat' i 0 u 1:4 w l lw 2 lc "gold" title 'Cumulative |M|/N (L=4)' 
#replot 'ist_tot.dat' i 2 u 1:4 w lp pt 13 lc "green" title 'Istantaneous |M|/N (L=30)'
#replot 'ist_cum_tot.dat' i 2 u 1:4 w l lw 2 lc "green" title 'Cumulative |M|/N (L=30)' 



#PBC, T = 4

set title 'Istantaneous and cumulative variables of interest at T=4J/k_B for PBC' font "Verdana,15"
#ENERGY

#set ylabel 'E/N' font "Verdana,13"
#set yrange[-2:0.8]
#set xrange[0:1000]


#p 'ist_tot.dat' i 1 u 1:2 w lp pt 5 lc "red" title 'Istantaneous E/N (L=4)'
#replot 'ist_cum_tot.dat' i 1 u 1:2 w l lw 2 lc "red" title 'Cumulative E/N (L=4)' 
#replot 'ist_tot.dat' i 3 u 1:2 w lp pt 13 lc "blue" title 'Istantaneous E/N (L=30)' 
#replot 'ist_cum_tot.dat' i 3 u 1:2 w l lw 2 lc "blue" title 'Cumulative E/N (L=30)' 

#MAGNETIZATION PBC

#set ylabel 'M/N' font "Verdana,13"
#set yrange[-1:1.5]
#set xrange[0:1000]

#plot 'ist_tot.dat' i 1 u 1:3 w lp pt 5 lc "gold" title 'Istantaneous M/N (L=4)'
#replot 'ist_cum_tot.dat' i 1 u 1:3 w l lw 2 lc "gold" title 'Cumulative M/N (L=4)' 
#replot 'ist_tot.dat' i 3 u 1:3 w lp pt 13 lc "green" title 'Istantaneous M/N (L=30)'
#replot 'ist_cum_tot.dat' i 3 u 1:3 w l lw 2 lc "green" title 'Cumulative M/N (L=30)' 


#ABS(MAGNETIZATION)

#set ylabel '|M|/N' font "Verdana,13"
#set yrange[-0.1:1.3]
#set xrange[0:1000]

#plot 'ist_tot.dat' i 1 u 1:4 w lp pt 5 lc "gold" title 'Istantaneous |M|/N (L=4)'
#replot 'ist_cum_tot.dat' i 1 u 1:4 w l lw 2 lc "gold" title 'Cumulative |M|/N (L=4)' 
#replot 'ist_tot.dat' i 3 u 1:4 w lp pt 13 lc "green" title 'Istantaneous |M|/N (L=30)'
#replot 'ist_cum_tot.dat' i 3 u 1:4 w l lw 2 lc "green" title 'Cumulative |M|/N (L=30)' 










#----------------------------------------------------------------------------------------------
#SPIN PLOT
#----------------------------------------------------------------------------------------------
#set xlabel 'x'
#set ylabel 'y'

#set title 'Equilibrium condition spins at T= 2J/k_{B} for PBC' font "Verdana,15"

#set xrange[0:5]
#set yrange[0:5]

#plot 'ising_up_tot.dat' i 0 u 1:2 w p pt 5 lc "purple" title 'Spin up'
#replot 'ising_down_tot.dat' i 0 u 1:2 w p pt 5 lc "black" title 'Spin down'

#set xrange[-3:34]
#set yrange[-3:34]

#plot 'ising_up_tot.dat' i 2 u 1:2 w p pt 5 lc "purple" title 'Spin up'
#replot 'ising_down_tot.dat' i 2 u 1:2 w p pt 5 lc "black" title 'Spin down'




#set title 'Equilibrium condition spins at T= 4J/k_{B} for PBC' font "Verdana,15"

#set xrange[0:5]
#set yrange[0:5]

#plot 'ising_up_tot.dat' i 1 u 1:2 w p pt 5 lc "purple" title 'Spin up'
#replot 'ising_down_tot.dat' i 1 u 1:2 w p pt 5 lc "black" title 'Spin down'

#set xrange[-3:34]
#set yrange[-3:34]

#plot 'ising_up_tot.dat' i 3 u 1:2 w p pt 5 lc "purple" title 'Spin up'
#replot 'ising_down_tot.dat' i 3 u 1:2 w p pt 5 lc "black" title 'Spin down'










#----------------------------------------------------------------------------------------------
#PLOT VARYING TEMPERATURE
#----------------------------------------------------------------------------------------------
set title 'Variables of interest at equilibrium' font "Verdana,15"

#set terminal gif animate delay 20
#set output 'animation.gif'
#set encoding utf8

#set xrange[1:4]
set xlabel 'T' font "Verdana,13"


set arrow nohead from 2.269,-2 to 2.269,0.15 lw 3 lc "black"

set ylabel '<E>/N' font "Verdana,13"
set yrange[-2:0.15]
#plot 'eq_data_tot.dat' i 0 u 1:2:($7**0.5) w yerrorbars pt 5 lc "red" title '<E>/N for L=4'
plot 'eq_data_tot.dat' i 0 u 1:2 w lp pt 5 lc "red" title 'PBC L=4'
#replot 'eq_data_tot.dat' i 1 u 1:2:($7**0.5) w yerrorbars pt 13 lc "blue" title 'PBC L=30'
replot 'eq_data_tot.dat' i 1 u 1:2 w lp pt 13 lc "blue" title '<E>/N for L=30'

replot 'eq_data_tot.dat' i 2 u 1:2 w lp pt 5 lc "gold" title 'OBC L=4'
replot 'eq_data_tot.dat' i 3 u 1:2 w lp pt 13 lc "green" title 'OBC L=30'


#set ylabel '<|M|>/N' font "Verdana,13"
#set yrange[0:1.1]
#plot 'eq_data_tot.dat' i 0 u 1:6 w lp pt 5 lc "red" title 'L=4'
#replot 'eq_data_tot.dat' i 1 u 1:6 w lp pt 13 lc "blue" title 'L=30'

#replot 'eq_data_tot.dat' i 2 u 1:6 w lp pt 5 lc "gold" title 'L=4'
#replot 'eq_data_tot.dat' i 3 u 1:6 w lp pt 13 lc "green" title 'L=30'



#set yrange[0:2]
#set ylabel 'c' font "Verdana,13"
#plot 'eq_data_tot.dat' i 0 u 1:7 w lp pt 5 lc "red" title 'PBC L=4'
#replot 'eq_data_tot.dat' i 1 u 1:7 w lp pt 13 lc "blue" title 'PBC L=30'
#replot 'c_dev_tot.dat' i 0 u 1:2 w lp pt 5 lc "violet" title 'Derivative L=4'
#replot 'c_dev_tot.dat' i 1 u 1:2 w lp pt 13 lc "grey" title 'Derivative L=30'

#replot 'eq_data_tot.dat' i 2 u 1:7 w lp pt 5 lc "gold" title 'OBC L=4'
#replot 'eq_data_tot.dat' i 3 u 1:7 w lp pt 13 lc "green" title 'OBC L=30'


set yrange[0:18]
set ylabel 'χ' font "Verdana,13"
plot 'eq_data_tot.dat' i 0 u 1:8 w lp pt 5 lc "red" title 'PBC L=4'
replot 'eq_data_tot.dat' i 1 u 1:8 w lp pt 13 lc "blue" title 'PBC L=30'

replot 'eq_data_tot.dat' i 2 u 1:8 w lp pt 5 lc "gold" title 'OBC L=4'
replot 'eq_data_tot.dat' i 3 u 1:8 w lp pt 13 lc "green" title 'OBC L=30'






#----------------------------------------------------------------------------------------------
#PLOT OBC
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#PLOT VARYING MC STEP
#----------------------------------------------------------------------------------------------
#set ylabel ''
#set xlabel 'MC step' font "Verdana,13"
#set xrange[0:1000]


#T = 2

#set title 'Istantaneous and cumulative variables of interest at T= 2J/k_{B} for OBC' font "Verdana,15"

#ENERGY


#set yrange[-2:0.5]

#plot 'ist_tot.dat' i 4 u 1:2 w lp pt 5 lc "red" title 'Istantaneous E/N (L=4)'
#replot 'ist_cum_tot.dat' i 4 u 1:2 w l lw 2 lc "red" title 'Cumulative E/N (L=4)' 
#replot 'ist_tot.dat' i 6 u 1:2 w lp pt 13 lc "blue" title 'Istantaneous E/N (L=30)' 
#replot 'ist_cum_tot.dat' i 6 u 1:2 w l lw 2 lc "blue" title 'Cumulative E/N (L=30)' 


#MAGNETIZATION PBC

#set yrange[-0.2:1.85]

#replot 'ist_tot.dat' i 4 u 1:3 w lp pt 5 lc "gold" title 'Istantaneous M/N (L=4)'
#replot 'ist_cum_tot.dat' i 4 u 1:3 w l lw 2 lc "gold" title 'Cumulative M/N (L=4)' 
#replot 'ist_tot.dat' i 6 u 1:3 w lp pt 13 lc "green" title 'Istantaneous M/N (L=30)'
#replot 'ist_cum_tot.dat' i 6 u 1:3 w l lw 2 lc "green" title 'Cumulative M/N (L=30)' 


#ABS(MAGNETIZATION)

#set ylabel '|M|/N' font "Verdana,13"

#replot 'ist_tot.dat' i 4 u 1:4 w lp pt 5 lc "gold" title 'Istantaneous |M|/N (L=4)'
#replot 'ist_cum_tot.dat' i 4 u 1:4 w l lw 2 lc "gold" title 'Cumulative |M|/N (L=4)' 
#replot 'ist_tot.dat' i 6 u 1:4 w lp pt 13 lc "green" title 'Istantaneous |M|/N (L=30)'
#replot 'ist_cum_tot.dat' i 6 u 1:4 w l lw 2 lc "green" title 'Cumulative |M|/N (L=30)' 



#PBC, T = 4

#set title 'Istantaneous and cumulative variables of interest at T=4J/k_B for PBC' font "Verdana,15"
#ENERGY

#set ylabel 'E/N' font "Verdana,13"
#set yrange[-2:0.8]


#plot 'ist_tot.dat' i 5 u 1:2 w lp pt 5 lc "red" title 'Istantaneous E/N (L=4)'
#replot 'ist_cum_tot.dat' i 5 u 1:2 w l lw 2 lc "red" title 'Cumulative E/N (L=4)' 
#replot 'ist_tot.dat' i 7 u 1:2 w lp pt 13 lc "blue" title 'Istantaneous E/N (L=30)' 
#replot 'ist_cum_tot.dat' i 7 u 1:2 w l lw 2 lc "blue" title 'Cumulative E/N (L=30)' 

#MAGNETIZATION PBC

#set ylabel 'M/N' font "Verdana,13"
#set yrange[-1:1.5]
#set xrange[0:1000]

#plot 'ist_tot.dat' i 5 u 1:3 w lp pt 5 lc "gold" title 'Istantaneous M/N (L=4)'
#replot 'ist_cum_tot.dat' i 5 u 1:3 w l lw 2 lc "gold" title 'Cumulative M/N (L=4)' 
#replot 'ist_tot.dat' i 7 u 1:3 w lp pt 13 lc "green" title 'Istantaneous M/N (L=30)'
#replot 'ist_cum_tot.dat' i 7 u 1:3 w l lw 2 lc "green" title 'Cumulative M/N (L=30)' 


#ABS(MAGNETIZATION)

#set ylabel '|M|/N' font "Verdana,13"
#set yrange[-0.1:1.3]
#set xrange[0:1000]

#plot 'ist_tot.dat' i 5 u 1:4 w lp pt 5 lc "gold" title 'Istantaneous |M|/N (L=4)'
#replot 'ist_cum_tot.dat' i 5 u 1:4 w l lw 2 lc "gold" title 'Cumulative |M|/N (L=4)' 
#replot 'ist_tot.dat' i 7 u 1:4 w lp pt 13 lc "green" title 'Istantaneous |M|/N (L=30)'
#replot 'ist_cum_tot.dat' i 7 u 1:4 w l lw 2 lc "green" title 'Cumulative |M|/N (L=30)' 






#set title 'Average energy per spin deviation L distribution for T= 1 J/k_{B} and OBC' font "Verdana,15"

#set ylabel 'ln[abs(-2 - <E>/N)]' font "Verdana,13"
#set xlabel 'L' font "Verdana,13"
#set xrange[4:30]

#f(x) = a*x + b

#fit f(x) 'obc_try.dat' u 1:(log(2 + $2)) via a,b
#plot 'obc_try.dat' u 1:(log(2 + $2)) w lp pt 5 lc "red" title "Ensemble averages"
#replot f(x) lw 2 lc "black" title "fit"









#set terminal gif animate delay 20
#set output 'animation.gif'
#set encoding utf8

#stats 'eq_data.dat' nooutput
#do for [i=1:20] {
#plot 'eq_data.dat' index (i-1) u 1:2 pt 13 lc "red" title '<E>/N'
#replot 'eq_data.dat' index (i-1) u 1:4 pt 13 lc "blue" title '<M>/N'
#pause 10}
#do for [i=1:20] {plot 'eq_data.dat' index (i-1) u 1:7 pt 13 lc "green" title 'c'}
#do for [i=1:20] {plot 'eq_data.dat' index (i-1) u 1:8 pt 13 lc "purple" title 'c_{numerical derivative}'}
#do for [i=1:20] {plot 'eq_data.dat' index (i-1) u 1:9 pt 13 lc "gold" title 'χ'}

#ind = 0 
#a = 0

#do for [ind=0:19] {
#    pause 0.5
#    plot 'eq_data.dat' i ind u 1:2 pt 13 lc "blue" notitle
#}
