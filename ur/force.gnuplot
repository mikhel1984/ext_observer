# use gnuplot -p force.gnuplot
set datafile separator ','
plot 'external.csv' using 1:2 with lines, '' using 1:3 with lines, '' using 1:4 with lines, '' using 1:5 with lines, '' using 1:6 with lines, '' using 1:7 with lines
