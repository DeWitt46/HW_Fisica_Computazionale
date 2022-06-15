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
set pointsize 1.3
#set lw 1.3
CHOICE = 5  #0 FOR PLACE; 1-2 FOR AVG FITNESS; 3 - FOR TIME; 4 - FOR ACCURACY;
             #5 - FOR OPT; 

#----------------------------------------------------------------------------------------------
#PLACE
#----------------------------------------------------------------------------------------------

set style data boxplot
stats 'data_sim.dat' using 4:5 nooutput
f(x) = STATS_mean_x
#plot 'data_sim_PMX_p01.dat' u 1:4:xtic('Col 4')

plot f(x)
if(CHOICE == 0){
                set title 'Place of current run' font "Verdana,15"
                set ylabel 'y' font "Verdana,13"
                set xlabel 'x' font "Verdana,13"
                set key t r font "Verdana,12"
                set pointsize 3.
                
                set yrange[0:16]
                set xrange[0:16]
                
                plot for [indice=0:7] 'pos.dat' i indice u 1:2 pt 5 title ''.(indice+1)

                }

#----------------------------------------------------------------------------------------------
#AVERAGE FITNESS
#----------------------------------------------------------------------------------------------

if(CHOICE == 1){
                set title 'Best fitness distributions for some simulations - ATT\_48' font "Verdana,15"
                set ylabel 'Fitness' font "Verdana,13"
                set xlabel 'Number of generation' font "Verdana,13"
                set key t l font "Verdana,12"
                                
                #set yrange[0.3:0.65]
                #set xrange[0:5]
                
                plot 'data_gen_PMX_att48.dat' i 5 u 1:2 w lp pt 5 lc "red" title 'SIMULATION 6 CX2'
                #replot 'data_gen_CX2_att48.dat' i 6 u 1:2:3 w yerrorbars pt 5 lc "red" notitle
                #replot 'data_gen_OX_att48.dat' i 4 u 1:2 w lp pt 7 lc "blue" title 'SIMULATION 4 OX'
                #replot 'data_gen_OX_att48.dat' i 4 u 1:2:3 w yerrorbars pt 7 lc "blue" notitle
                replot 'data_gen_PMX_att48.dat' i 5 u 1:4 w lp pt 13 lc "green" title 'SIMULATION 5 PMX'
                #replot 'data_gen_PMX_att48.dat' i 5 u 1:2:3 w yerrorbars pt 13 lc "green" notitle
                
                
                }

if(CHOICE == 2){
                set title 'Average fitness of best individuals per generation distributions - ATT\_48' font "Verdana,15"
                set ylabel 'Fitness' font "Verdana,13"
                set xlabel 'Number of simulation' font "Verdana,13"
                set key t r font "Verdana,12"
                                
                #set yrange[9.4:9.9]
                #set xrange[0:5]
                
                plot 'data_sim_CX2_att48.dat' u 1:2 w lp pt 5 lc "red" title 'CX2'
                #replot 'data_sim_CX2_att48.dat' u 1:2:3 w yerrorbars pt 5 lc "red" notitle
                replot 'data_sim_OX_att48.dat' u 1:2 w lp pt 7 lc "blue" title 'OX'
                #replot 'data_sim_OX_att48.dat' u 1:2:3 w yerrorbars pt 7 lc "blue" notitle
                replot 'data_sim_PMX_att48.dat' u 1:2 w lp pt 13 lc "green" title 'PMX'
                #replot 'data_sim_PMX_att48.dat' u 1:2:3 w yerrorbars pt13 lc "green" notitle
                
                
                }
                
                

#----------------------------------------------------------------------------------------------
# GA - BF COMPARISON
#----------------------------------------------------------------------------------------------
if(CHOICE == 3){
                set title 'Time trend comparison' font "Verdana,15"
                set ylabel 'Time (s)' font "Verdana,13"
                set xlabel 'N' font "Verdana,13"


                #set yrange[0.15:0.35]
                #set xrange[0:1210]


                plot 'analysis.dat' i 0 u 1:2 w lp lw 2 pt 5 lc "red" title 'CX2'
                replot 'analysis.dat' i 0 u 1:2:3 w yerrorbars lw 2 lc "red" notitle

                replot 'analysis.dat' i 1 u 1:2 w lp  lw 2 pt 7 lc "blue" title 'Brute-force'
                replot 'analysis.dat' i 1 u 1:2:3 w yerrorbars lw 2 lc "blue" notitle               
                }

if(CHOICE == 4){
                set title 'Optimum distribution' font "Verdana,15"
                set ylabel 'Optimum' font "Verdana,13"
                set xlabel 'N' font "Verdana,13"


                #set yrange[0.15:0.35]
                #set xrange[0:1210]


                plot 'analysis.dat' i 0 u 1:4 w lp lw 2 pt 5 lc "red" title 'CX2'
                replot 'analysis.dat' i 0 u 1:4:5 w yerrorbars lw 2 lc "red" notitle
                }

if(CHOICE == 5){
                set title 'Average fitness displacement from optimal value' font "Verdana,15"
                set ylabel '{/Symbol D} value/fitness (percentage)' font "Verdana,13"
                set xlabel 'N' font "Verdana,13"


                #set yrange[0.15:0.35]
                set xrange[4:11.1]


                plot 'analysis.dat' i 0 u 1:6 w lp lw 2 pt 5 lc "red" title 'CX2'
                replot 'analysis.dat' i 0 u 1:6:7 w yerrorbars lw 2 lc "red" notitle
                }
