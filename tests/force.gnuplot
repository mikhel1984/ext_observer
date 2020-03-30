# use gnuplot -p force.gnuplot
set datafile separator ','
plot 'force.csv' using 1:2 with lines, '' using 1:3 with lines
