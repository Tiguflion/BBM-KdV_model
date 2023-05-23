 set title "Hauteur d'eau à T =   0.99032365699031699      "
 set terminal png
 set output "frame/plot.0100.png"
 plot "data/data0100.dat" using 1:2 title "h app"w l, "data/data0100.dat" using 1:3 title "h theo" w l
