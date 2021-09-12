
build:
	g++ Mandelbrot.cpp -std=c++0x -o Mandelbrot

run:
	./Mandelbrot

clean:
	rm -f Mandelbrot

install:
	./install.sh