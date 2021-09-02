// To compile: g++ Mandelbrot.cc -std=c++0x -o Mandelbrot

#include <iostream>
#include <sys/time.h>
#include <stdexcept>

using namespace std;

class Timer {
	timeval* _start_time;
    timeval* _stop_time;

public:

    Timer() {
		this->_start_time = nullptr;
        this->_stop_time = nullptr;
	}

    ~Timer() {
		if (this->_start_time != nullptr) {delete this->_start_time;}
		if (this->_stop_time != nullptr) {delete this->_stop_time;}
	}

    void start() {
        //Start a new timer
        if (this->_start_time != nullptr) {
            throw std::invalid_argument("Timer is running. Use .stop() to stop it");
		}
        this->_stop_time = nullptr;
		this->_start_time = new struct timeval;
		gettimeofday(this->_start_time,nullptr);
	}

    void stop() {
        //Stop the timer, and report the elapsed time
        if (this->_start_time == nullptr) {
            throw std::invalid_argument("Timer is not running. Use .start() to start it");
		}
		this->_stop_time = new struct timeval;
        gettimeofday(this->_stop_time,nullptr);
	}

	double getElapsed_time() {
        double elapsed_time = (this->_stop_time->tv_sec - this->_start_time->tv_sec)*1000000;
		elapsed_time += (this->_stop_time->tv_usec - this->_start_time->tv_usec);
		elapsed_time /= 1000000;
        delete this->_start_time;
		delete this->_stop_time;
		this->_start_time = nullptr;
        this->_stop_time = nullptr;
        return elapsed_time;
	}
};

class MandelbrotPoint {
	int INFINITY = 2;
    int MAX_ITER = 1000;
	long double x;
	long double y;

public:

	//Default constructor
	MandelbrotPoint(double x, double y) {
		this->x = x;
		this->y = y;
	}

	//Should return triplet of color code when finished
	int getColor() {
		long double a = 0;
        long double b = 0;
        int i = 0;
        int N = 0;
        long double a2 = a*a;
        long double b2 = b*b;
        int inf = this->INFINITY * this->INFINITY;
        while (N < inf and i < this->MAX_ITER) {
            long double aa = a;
            a = a2-b2 + this->x;
            b = 2*aa*b + this->y;
            i+=1;
            a2 = a*a;
            b2 = b*b;
            N = a2+b2;
		}
        return i;
	}

};


class MandelbrotSet {
	long double _xmin;
	long double _xmax;
	long double _ymin;
	long double _ymax;
	int _width;
	int _height;

public:

    //x,y coordinate of center of image, width and height the number of pixels to compute
    MandelbrotSet(long double xmin,long double xmax,long double ymin,long double ymax,int width=640,int height=480) {
        this->_xmin = xmin;
        this->_xmax = xmax;
        this->_ymin = ymin;
        this->_ymax = ymax;
        this->_width = width;
        this->_height = height;
        //this->_data = nullptr;
	}

    void compute() {
        //this->_data = np.zeros((this->_width, this->_height, 3), dtype=np.uint8)
        long double xstep = (this->_xmax - this->_xmin) / this->_width;
        long double ystep = (this->_ymax - this->_ymin) / this->_height;
        int i = 0;
        while (i<this->_width) {
            int j = 0;
            while (j<this->_height) {
                MandelbrotPoint m = MandelbrotPoint(this->_xmin+i*xstep, this->_ymin+j*ystep);
                m.getColor();
				//this->_data[i, j] = m.getColor()
                j+=1;
			}
            i+=1;
		}
	}

    void exportToImage() {
        //return Image.fromarray(this->_data)
	}

};


int main() {
	int NB_TEST = 10;
    double total_time = 0;
	for (size_t n = 0; n < NB_TEST; n++)
	{	
		Timer t = Timer();
		t.start();
		MandelbrotSet mSet = MandelbrotSet(-2,1,-1.5,1.1000,1000);
		mSet.compute();
		t.stop();
		double e = t.getElapsed_time();
		std::cout << "Boucle " << n << " ran in " << e << "seconds" << std::endl;
		total_time += e;
	}
	std::cout << "Have run  " << NB_TEST << " Mandelbrot Set calculations in " << total_time/NB_TEST << " seconds average" << std::endl;
	return 0;
}




/*

int mandelbrot(double real, double imag) {
	int limit = 100;
	double zReal = real;
	double zImag = imag;

	for (int i = 0; i < limit; ++i) {
		double r2 = zReal * zReal;
		double i2 = zImag * zImag;
		
		if (r2 + i2 > 4.0) return i;

		zImag = 2.0 * zReal * zImag + imag;
		zReal = r2 - i2 + real;
	}
	return limit;
}


int main() {
	
	int width = 379; //number of characters fitting horizontally on my screen 
	int heigth = 98; //number of characters fitting vertically on my screen
		
	double x_start = -2.0;
	double x_fin = 1.0;
	double y_start = -1.0;
	double y_fin = 1.0;
	
	//~ double x_start = -0.25;
	//~ double x_fin = 0.05;
	//~ double y_start = -0.95;
	//~ double y_fin = -0.75;
	
	//~ double x_start = -0.13;
	//~ double x_fin = -0.085;
	//~ double y_start = -0.91;
	//~ double y_fin = -0.88;
	
	double dx = (x_fin - x_start)/(width - 1);
	double dy = (y_fin - y_start)/(heigth - 1);

	string char_ = "\u2588";

	string black = "\033[22;30m";
	string red = "\033[22;31m";
	string l_red = "\033[01;31m";
	string green = "\033[22;32m";
	string l_green = "\033[01;32m";
	string orange = "\033[22;33m";
	string yellow = "\033[01;33m";
	string blue = "\033[22;34m";
	string l_blue = "\033[01;34m";
	string magenta = "\033[22;35m";
	string l_magenta = "\033[01;35m";
	string cyan = "\033[22;36m";
	string l_cyan = "\033[01;36m";
	string gray = "\033[22;37m";
	string white = "\033[01;37m";

	for (int i = 0; i < heigth; i++) {
		for (int j = 0; j < width; j++) {
			
			double x = x_start + j*dx; // current real value
			double y = y_fin - i*dy; // current imaginary value
			
			int value = mandelbrot(x,y);
			
			if (value == 100) {cout << " ";}
			else if (value > 90) {cout << red << char_;}
			else if (value > 70) {cout << l_red << char_;}
			else if (value > 50) {cout << orange << char_;}
			else if (value > 30) {cout << yellow << char_;}
			else if (value > 20) {cout << l_green << char_;}
			else if (value > 10) {cout << green << char_;}
			else if (value > 5) {cout << l_cyan << char_;}
			else if (value > 4) {cout << cyan << char_;}
			else if (value > 3) {cout << l_blue << char_;}
			else if (value > 2) {cout << blue << char_;}
			else if (value > 1) {cout << magenta << char_;}
			else {cout << l_magenta << char_;}
			
			cout << "\033[0m";
		}
		cout << endl;
	}

	return 0;
}

*/