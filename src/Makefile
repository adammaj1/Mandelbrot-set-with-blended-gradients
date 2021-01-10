all: 
	printf "make ppm files"
	gcc d.c -lm -Wall -march=native -fopenmp
	time ./a.out > a.txt
	printf "convert all ppm files to png using Image Magic convert"
	chmod +x g.sh
	./g.sh
	printf "delete all ppm and txt files"
	rm *.ppm
	
