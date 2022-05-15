set termoption enhanced

#-------------------------------------------------
#GNUplot script to plot data of es_1
#-------------------------------------------------


set key t r font "Verdana,13"



set xtics font "Verdana,13"
set ytics font "Verdana,13"


set autoscale
#set xrange[-0.2:0.2]
#set yrange[0:1]
set pointsize 1.3
#set lw 1.3
CHOICE = 11  #0-1 FOR DISPLACEMENT; 2-3 FOR D; 4-5 FOR D-VARIANCE; 6-7 FOR D AVG;
             #8-9 FOR D ERR; 10 FOR D-ERR ON N_P; 11 FOR D ON RHO

#----------------------------------------------------------------------------------------------
#DISPLACEMENT
#----------------------------------------------------------------------------------------------
if(CHOICE == 0){
                set title 'Istantaneous mean square displacement per particles for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel 'ln[<R(t)^2>]' font "Verdana,13"
                set xlabel 'ln[MC step]' font "Verdana,13"
                set key t l font "Verdana,12"
                
                
                set yrange[0:5]
                set xrange[0:5]
                
                f1(x) = a*x + b
                f2(x) = c*x + d
                f3(x) = e*x + f
                f4(x) = g*x + h
                f5(x) = i*x + j
                
                fit f1(x) 'disp_data_tot.dat' i 0 u (log($1)):(log($2)) via a,b
                plot 'disp_data.dat' i 0 u (log($1)):(log($2)) w p pt 5 lc "red" title 'RUNNER 1'
                replot f1(x) lw 3 lc "red" notitle
                
                fit f2(x) 'disp_data_tot.dat' i 1 u (log($1)):(log($2)) via c,d
                replot 'disp_data.dat' i 1 u (log($1)):(log($2)) w p pt 7 lc "blue" title 'RUNNER 2'
                replot f2(x) lw 3 lc "blue" notitle

                fit f3(x) 'disp_data_tot.dat' i 2 u (log($1)):(log($2)) via e,f
                replot 'disp_data.dat' i 2 u (log($1)):(log($2)) w p pt 9 lc "gold" title 'RUNNER 3'
                replot f3(x) lw 3 lc "gold" notitle

                fit f4(x) 'disp_data_tot.dat' i 3 u (log($1)):(log($2)) via g,h
                replot 'disp_data.dat' i 3 u (log($1)):(log($2)) w p pt 13 lc "green" title 'RUNNER 4'
                replot f4(x) lw 3 lc "green" notitle

                fit f5(x) 'disp_data_tot.dat' i 4 u (log($1)):(log($2)) via i,j
                replot 'disp_data.dat' i 4 u (log($1)):(log($2)) w p pt 15 lc "purple" title 'RUNNER 5'
                replot f5(x) lw 3 lc "purple" notitle

                }

if(CHOICE == 1){
                set title 'Istantaneous mean square displacement per particles for {/Symbol r} = 0.2' font "Verdana,15"
                set ylabel 'ln[<R(t)^2>]' font "Verdana,13"
                set xlabel 'ln[MC step]' font "Verdana,13"
                set key t l font "Verdana,12"
                                
                set yrange[0:5]
                set xrange[0:5]
                
                f1(x) = a*x + b
                f2(x) = c*x + d
                f3(x) = e*x + f
                f4(x) = g*x + h
                f5(x) = i*x + j
                
                fit f1(x) 'disp_data_tot.dat' i 0 u (log($1)):(log($2)) via a,b
                plot 'disp_data.dat' i 0 u (log($1)):(log($2)) w p pt 5 lc "red" title 'RUNNER 1'
                replot f1(x) lw 3 lc "red" notitle
                
                fit f2(x) 'disp_data_tot.dat' i 1 u (log($1)):(log($2)) via c,d
                replot 'disp_data.dat' i 1 u (log($1)):(log($2)) w p pt 7 lc "blue" title 'RUNNER 2'
                replot f2(x) lw 3 lc "blue" notitle

                fit f3(x) 'disp_data_tot.dat' i 2 u (log($1)):(log($2)) via e,f
                replot 'disp_data.dat' i 2 u (log($1)):(log($2)) w p pt 9 lc "gold" title 'RUNNER 3'
                replot f3(x) lw 3 lc "gold" notitle

                fit f4(x) 'disp_data_tot.dat' i 3 u (log($1)):(log($2)) via g,h
                replot 'disp_data.dat' i 3 u (log($1)):(log($2)) w p pt 13 lc "green" title 'RUNNER 4'
                replot f4(x) lw 3 lc "green" notitle

                fit f5(x) 'disp_data_tot.dat' i 4 u (log($1)):(log($2)) via i,j
                replot 'disp_data.dat' i 4 u (log($1)):(log($2)) w p pt 15 lc "purple" title 'RUNNER 5'
                replot f5(x) lw 3 lc "purple" notitle

                }

#----------------------------------------------------------------------------------------------
#ISTANTANEOUS SELF-DIFFUSIVE COEFFICIENT
#----------------------------------------------------------------------------------------------
if(CHOICE == 2){
                set title 'Istantaneous self-diffusive coefficient per particles for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel 'D(t)' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                set yrange[0.15:0.35]
                set xrange[0:1210]


                plot 'disp_data.dat' i 0 u 1:3 w p pt 5 lc "red" title 'RUNNER 1'
                replot 'disp_data_tot.dat' i 0 u 1:3 w l lw 2 lc "red" notitle
                replot 'disp_data.dat' i 1 u 1:3 w p pt 7 lc "blue" title 'RUNNER 2'
                replot 'disp_data_tot.dat' i 1 u 1:3 w l lw 2 lc "blue" notitle               
                replot 'disp_data.dat' i 2 u 1:3 w p pt 9 lc "gold" title 'RUNNER 3'
                replot 'disp_data_tot.dat' i 2 u 1:3 w l lw 2 lc "gold" notitle
                replot 'disp_data.dat' i 3 u 1:3 w p pt 13 lc "green" title 'RUNNER 4'
                replot 'disp_data_tot.dat' i 3 u 1:3 w l lw 2 lc "green" notitle                
                replot 'disp_data.dat' i 4 u 1:3 w p pt 15 lc "purple" title 'RUNNER 5'
                replot 'disp_data_tot.dat' i 4 u 1:3 w l lw 2 lc "purple" notitle
}

if(CHOICE == 3){
                set title 'Istantaneous self-diffusive coefficient per particles for {/Symbol r} = 0.2' font "Verdana,15"
                set ylabel 'D(t)' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                set yrange[0.05:0.5]
                #set xrange[0:1000]


                plot 'disp_data.dat' i 0 u 1:3 w p pt 5 lc "red" title 'RUNNER 1'
                replot 'disp_data_tot.dat' i 0 u 1:3 w l lw 2 lc "red" notitle
                replot 'disp_data.dat' i 1 u 1:3 w p pt 7 lc "blue" title 'RUNNER 2'
                replot 'disp_data_tot.dat' i 1 u 1:3 w l lw 2 lc "blue" notitle               
                replot 'disp_data.dat' i 2 u 1:3 w p pt 9 lc "gold" title 'RUNNER 3'
                replot 'disp_data_tot.dat' i 2 u 1:3 w l lw 2 lc "gold" notitle
                replot 'disp_data.dat' i 3 u 1:3 w p pt 13 lc "green" title 'RUNNER 4'
                replot 'disp_data_tot.dat' i 3 u 1:3 w l lw 2 lc "green" notitle                
                replot 'disp_data.dat' i 4 u 1:3 w p pt 15 lc "purple" title 'RUNNER 5'
                replot 'disp_data_tot.dat' i 4 u 1:3 w l lw 2 lc "purple" notitle

}



#----------------------------------------------------------------------------------------------
#ISTANTANEOUS SELF-DIFFUSIVE COEFFICIENT AVERAGE
#----------------------------------------------------------------------------------------------
if(CHOICE == 4){
                set title 'Istantaneous self-diffusive coefficient per particles variance for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel '{/Symbol s}^2' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.05:0.5]
                #set xrange[0:1000]


                plot 'disp_data.dat' i 0 u 1:5 w lp lw 2 pt 5 lc "red" title 'RUNNER 1'
                replot 'disp_data.dat' i 1 u 1:5 w lp lw 2 pt 7 lc "blue" title 'RUNNER 2'
                replot 'disp_data.dat' i 2 u 1:5 w lp lw 2 pt 9 lc "gold" title 'RUNNER 3'
                replot 'disp_data.dat' i 3 u 1:5 w lp lw 2 pt 13 lc "green" title 'RUNNER 4'
                replot 'disp_data.dat' i 4 u 1:5 w lp lw 2 pt 15 lc "purple" title 'RUNNER 5'
}

if(CHOICE == 5){
                set title 'Istantaneous self-diffusive coefficient per particles variance for {/Symbol r} = 0.2' font "Verdana,15"
                set ylabel '{/Symbol s}^2' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.05:0.5]
                #set xrange[0:1000]


                plot 'disp_data.dat' i 0 u 1:5 w lp lw 2 pt 5 lc "red" title 'RUNNER 1'
                replot 'disp_data.dat' i 1 u 1:5 w lp lw 2 pt 7 lc "blue" title 'RUNNER 2'
                replot 'disp_data.dat' i 2 u 1:5 w lp lw 2 pt 9 lc "gold" title 'RUNNER 3'
                replot 'disp_data.dat' i 3 u 1:5 w lp lw 2 pt 13 lc "green" title 'RUNNER 4'
                replot 'disp_data.dat' i 4 u 1:5 w lp lw 2 pt 15 lc "purple" title 'RUNNER 5'
}


#----------------------------------------------------------------------------------------------
#ISTANTANEOUS SELF-DIFFUSIVE COEFFICIENT AVERAGE DISTRIBUTION COMPARISON
#----------------------------------------------------------------------------------------------
if(CHOICE == 6){
                set title 'Self-diffusive coefficient per particles Average comparison for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel '<D(t)>' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.05:0.5]
                set xrange[0:600]


                plot 'disp_data.dat' i 0 u 1:3 w p pt 5 lc "black" title 'D(t)'
                replot 'disp_data_tot.dat' i 0 u 1:3 w l lw 2 lc "black" notitle
                replot 'disp_data.dat' i 0 u 1:4:5 w yerrorbars pt 5 lc "cyane" title '<D(t)>_t'
                replot 'disp_data_tot.dat' i 0 u 1:4 w l lw 2 lc "cyane" notitle
                replot 'disp_data.dat' i 0 u 1:6:7 w yerrorbars pt 5 lc "violet" title '<D(t)>^{(I)}'
                replot 'disp_data_tot.dat' i 0 u 1:6 w l lw 2 lc "violet" notitle

}

if(CHOICE == 7){
                set title 'Self-diffusive coefficient per particles Average comparison for {/Symbol r} = 0.2' font "Verdana,15"
                set ylabel '<D(t)>' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.05:0.5]
                #set xrange[0:600]


                plot 'disp_data.dat' i 0 u 1:3 w p pt 5 lc "black" title 'D(t)'
                replot 'disp_data_tot.dat' i 0 u 1:3 w l lw 2 lc "black" notitle
                replot 'disp_data.dat' i 0 u 1:4:5 w yerrorbars pt 5 lc "cyane" title '<D(t)>_t'
                replot 'disp_data_tot.dat' i 0 u 1:4 w l lw 2 lc "cyane" notitle
                replot 'disp_data.dat' i 0 u 1:6:7 w yerrorbars pt 5 lc "violet" title '<D(t)>^{(I)}'
                replot 'disp_data_tot.dat' i 0 u 1:6 w l lw 2 lc "violet" notitle
}


#----------------------------------------------------------------------------------------------
#ISTANTANEOUS SELF-DIFFUSIVE COEFFICIENT
#----------------------------------------------------------------------------------------------
if(CHOICE == 8){
                set title 'Istantaneous self-diffusive coefficient per particles Averages error comparison for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel '{/Symbol s}' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.15:0.35]
                #set xrange[0:1210]


                plot 'disp_data.dat' i 0 u 1:5 w p pt 5 lc "cyane" title '{/Symbol s}_{<D(t)>_t}'
                replot 'disp_data_tot.dat' i 0 u 1:5 w l lw 2 lc "cyane" notitle
                replot 'disp_data.dat' i 0 u 1:7 w p pt 7 lc "violet" title '{/Symbol s}_{<D(t)>^{(I)}}'
                replot 'disp_data_tot.dat' i 0 u 1:7 w l lw 2 lc "violet" notitle               

}

if(CHOICE == 9){
                set title 'Istantaneous self-diffusive coefficient per particles Averages error comparison for {/Symbol r} = 0.03' font "Verdana,15"
                set ylabel '{/Symbol s}' font "Verdana,13"
                set xlabel 'MC step' font "Verdana,13"


                #set yrange[0.15:0.35]
                #set xrange[0:1210]


                plot 'disp_data.dat' i 4 u 1:5 w p pt 5 lc "cyane" title '{/Symbol s}_{<D(t)>_t}'
                replot 'disp_data_tot.dat' i 4 u 1:5 w l lw 2 lc "cyane" notitle
                replot 'disp_data.dat' i 4 u 1:7 w p pt 7 lc "violet" title '{/Symbol s}_{<D(t)>^{(I)}}'
                replot 'disp_data_tot.dat' i 4 u 1:7 w l lw 2 lc "violet" notitle               


}


#----------------------------------------------------------------------------------------------
#SELF-DIFFUSIVE COEFFICIENT ERROR COMPARISON ON NUMBER OF PARTICLES
#----------------------------------------------------------------------------------------------
if(CHOICE == 10){
                set title 'Self-diffusive coefficient per particles Averages error comparison' font "Verdana,15"
                set ylabel 'ln[{/Symbol s}]' font "Verdana,13"
                set xlabel 'ln[N]' font "Verdana,13"


                #set yrange[0.15:0.35]
                #set xrange[0:1210]

                f(x) = a*x + b
                g(x) = c*x + d
                
                fit f(x) 'N_data.dat' i 0 u (log($1)):(log($2)) via a,b
                plot 'N_data.dat' i 0 u (log($1)):(log($2)) w p pt 5 lc "cyane" title '{/Symbol s}_{<D(t)>_t}'
                replot f(x) lw 2 lc "cyane" notitle

                fit g(x) 'N_data.dat' i 0 u (log($1)):(log($3)) via c,d
                replot 'N_data.dat' i 0 u (log($1)):(log($3)) w p pt 7 lc "violet" title '{/Symbol s}_{<D(t)>^{(I)}}'
                replot g(x) lw 2 lc "violet" notitle
}




#----------------------------------------------------------------------------------------------
#SELF-DIFFUSIVE COEFFICIENT ON DENSITY
#----------------------------------------------------------------------------------------------
if(CHOICE == 11){
                set title 'Self-diffusive coefficient per particles on a 50 x 50 lattice' font "Verdana,15"
                set ylabel 'D' font "Verdana,13"
                set xlabel '{/Symbol r}' font "Verdana,13"


                set yrange[0.1:0.25]
                set xrange[0.09:0.73]


                plot 'N_data.dat' i 1 u 1:2:3 w yerrorbars pt 5 lc "red" notitle
                replot 'D.dat' u 1:2 w l lw 2 lc "black" notitle

}

