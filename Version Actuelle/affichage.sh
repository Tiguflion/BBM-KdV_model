#!/bin/bash 

function tout(){
cd bash 
       gnuplot -persist hauteur_eau.gnu
       gnuplot -persist debit.gnu
       #gnuplot -persist energie.gnu
       #gnuplot -persist debit_vert.gnu
       #gnuplot -persist sigma.gnu
       #gnuplot -persist maree.gnu
cd ../
       sortie_video;
       echo "Affichage terminé"
}      

function sortie_video(){
	rm -f out.mp4
	./plot
	ffmpeg -r 20 -i frame/plot.%04d.png -pix_fmt yuv420p out.mp4
};



tout;

