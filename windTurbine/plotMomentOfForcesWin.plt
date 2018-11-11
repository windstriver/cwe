# Script for gnuplot

set title "Pressure Moment of Forces in -x direction"
set xlabel "Time [s]"
set ylabel "Moment of Forces [N.m]"
set yrange [0:300000]


plot "postProcessing/turbine-forces-and-moments/0/forcesWithoutBrackets.dat" using 1:(-$11) title "Moment of Forces (-x direction)" with lines,
pause mouse