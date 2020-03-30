# use gnuplot -p kalman.gnuplot
set datafile separator ','
plot 'estimation.csv' using 1:2 with lines, '' using 1:3 with lines
