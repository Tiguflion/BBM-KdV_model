 set title "Hauteur d'eau à T =    9.9759975997595838      "
 set terminal png
 set output "frame/plot.0400.png"
 plot "data/data0400.dat" using 1:2 title "u"w l, "data/data0400.dat" using 1:3 title "w" w l, "data/data0400.dat" u 4:5 title "domaine" w l linetype 4 linewidth 3
