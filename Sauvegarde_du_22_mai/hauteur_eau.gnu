 set title "Hauteur d'eau à T =    1.6016016016016023E-002 "
 set terminal png
 set output "frame/plot.0020.png"
 plot "data/data0020.dat" using 1:2 title "h app"w l, "data/data0020.dat" using 1:3 title "h theo" w l
