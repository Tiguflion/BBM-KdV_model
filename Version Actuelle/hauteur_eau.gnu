 set title "Hauteur d'eau à T =    10.000000000000000      "
 set terminal png
 set output "frame/plot.0101.png"
 plot "data/data0101.dat" using 1:2 title "u"w l, "data/data0101.dat" using 1:3 title "w" w l, "data/data0101.dat" u 4:5 title "domaine" w l linetype 4 linewidth 3
