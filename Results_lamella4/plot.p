reset


set style line 1 lc rgb "#ff0000" lt 1 lw 1.5
set style line 2 lc rgb "#ffff00" lt 1 lw 1.5
set style line 3 lc rgb "#00ff00" lt 1 lw 1.5

  
set xrange [0:8]
set yrange [0:8]
set zrange [0:1.1]
set xlabel "r"
set ylabel "z"
set zlabel "phi"

	splot "phi.dat" using 1:2:3 ls 1 title "phiA", \
"phi.dat" using 1:2:4 ls 2 title "phiB", \
"phi.dat" using 1:2:5 ls 3 title "phiC", \
