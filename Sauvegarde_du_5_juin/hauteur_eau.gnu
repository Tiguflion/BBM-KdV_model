 set title "Hauteur d'eau à T =   0.99833166499832482      "
 set terminal png
 set output "frame/plot.0400.png"
 plot "data/data0400.dat" using 1:2 title "h app"w l, "data/data0400.dat" using 1:3 title "h theo" w l
