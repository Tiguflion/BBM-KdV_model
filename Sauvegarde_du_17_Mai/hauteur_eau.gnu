 set title "Hauteur d'eau à T =   0.25000000000000000      "
 set terminal png
 set output "frame/plot.0500.png"
 plot "data/data0500.dat" using 1:2 title "h app"w l, "data/data0500.dat" using 1:3 title "h theo" w l
