// To compile: g++ Mandelbrot.cc -std=c++0x -o Mandelbrot

#include <iostream>
#include <string>
#include <iostream>
#include <sys/time.h>
#include <stdexcept>
#include <stdio.h>

using namespace std;

////////////////////////////////////////////////////////
// Bitmap file creation helper functions
////////////////////////////////////////////////////////

const int BYTES_PER_PIXEL = 3; /// red, green, & blue
const int FILE_HEADER_SIZE = 14;
const int INFO_HEADER_SIZE = 40;

void generateBitmapImage(unsigned char* image, int height, int width, const char* imageFileName);
unsigned char* createBitmapFileHeader(int height, int stride);
unsigned char* createBitmapInfoHeader(int height, int width);

void generateBitmapImage (unsigned char* image, int height, int width, const char* imageFileName)
{
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

	remove(imageFileName);
    FILE* imageFile = fopen(imageFileName, "wb");

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i*widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

unsigned char* createBitmapFileHeader (int height, int stride)
{
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0,0,     /// signature
        0,0,0,0, /// image file size in bytes
        0,0,0,0, /// reserved
        0,0,0,0, /// start of pixel array
    };

    fileHeader[ 0] = (unsigned char)('B');
    fileHeader[ 1] = (unsigned char)('M');
    fileHeader[ 2] = (unsigned char)(fileSize      );
    fileHeader[ 3] = (unsigned char)(fileSize >>  8);
    fileHeader[ 4] = (unsigned char)(fileSize >> 16);
    fileHeader[ 5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

unsigned char* createBitmapInfoHeader (int height, int width)
{
    static unsigned char infoHeader[] = {
        0,0,0,0, /// header size
        0,0,0,0, /// image width
        0,0,0,0, /// image height
        0,0,     /// number of color planes
        0,0,     /// bits per pixel
        0,0,0,0, /// compression
        0,0,0,0, /// image size
        0,0,0,0, /// horizontal resolution
        0,0,0,0, /// vertical resolution
        0,0,0,0, /// colors in color table
        0,0,0,0, /// important color count
    };

    infoHeader[ 0] = (unsigned char)(INFO_HEADER_SIZE);
    infoHeader[ 4] = (unsigned char)(width      );
    infoHeader[ 5] = (unsigned char)(width >>  8);
    infoHeader[ 6] = (unsigned char)(width >> 16);
    infoHeader[ 7] = (unsigned char)(width >> 24);
    infoHeader[ 8] = (unsigned char)(height      );
    infoHeader[ 9] = (unsigned char)(height >>  8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(BYTES_PER_PIXEL*8);

    return infoHeader;
}

////////////////////////////////////////////////////////
// Timer Functions to monitor performance
////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////
// Mandlebrot Set and Point functions
////////////////////////////////////////////////////////

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

    void colorOrbite0(char* color, int n, long double a, long double b) {
        //The color map from wikipedia Mandelbrot article
        if (n < this->MAX_ITER) {
			color[0] = (char)255;
			color[1] = (char)255;
			color[2] = (char)255;
		}
        else {
			color[0] = (char)0;
			color[1] = (char)0;
			color[2] = (char)0;
		}
	}

    void colorOrbite1(char* color, int n, long double a, long double b) {
        //The simplest one Black or White
        if (n < this->MAX_ITER) {
			int i = n % 16;
			switch ( i )
			{
			case 0:
				color[0] = (unsigned char)66;
				color[1] = (unsigned char)30;
				color[2] = (unsigned char)15;
				break;

			case 1:
				color[0] = (unsigned char)25;
				color[1] = (unsigned char)7;
				color[2] = (unsigned char)26;
				break;

			case 2:
				color[0] = (unsigned char)9;
				color[1] = (unsigned char)1;
				color[2] = (unsigned char)47;
				break;

			case 3:
				color[0] = (unsigned char)4;
				color[1] = (unsigned char)4;
				color[2] = (unsigned char)73;
				break;

			case 4:
				color[0] = (unsigned char)0;
				color[1] = (unsigned char)7;
				color[2] = (unsigned char)100;
				break;

			case 5:
				color[0] = (unsigned char)12;
				color[1] = (unsigned char)44;
				color[2] = (unsigned char)138;
				break;

			case 6:
				color[0] = (unsigned char)24;
				color[1] = (unsigned char)82;
				color[2] = (unsigned char)177;
				break;

			case 7:
				color[0] = (unsigned char)57;
				color[1] = (unsigned char)125;
				color[2] = (unsigned char)209;
				break;

			case 8:
				color[0] = (unsigned char)134;
				color[1] = (unsigned char)181;
				color[2] = (unsigned char)229;
				break;

			case 9:
				color[0] = (unsigned char)211;
				color[1] = (unsigned char)236;
				color[2] = (unsigned char)248;
				break;

			case 10:
				color[0] = (unsigned char)241;
				color[1] = (unsigned char)233;
				color[2] = (unsigned char)191;
				break;

			case 11:
				color[0] = (unsigned char)248;
				color[1] = (unsigned char)201;
				color[2] = (unsigned char)95;
				break;

			case 12:
				color[0] = (unsigned char)255;
				color[1] = (unsigned char)170;
				color[2] = (unsigned char)0;
				break;

			case 13:
				color[0] = (unsigned char)204;
				color[1] = (unsigned char)128;
				color[2] = (unsigned char)0;
				break;

			case 14:
				color[0] = (unsigned char)153;
				color[1] = (unsigned char)87;
				color[2] = (unsigned char)0;
				break;

			case 15:
				color[0] = (unsigned char)106;
				color[1] = (unsigned char)52;
				color[2] = (unsigned char)3;
				break;

			default:
				color[0] = (char)0;
				color[1] = (char)0;
				color[2] = (char)0;
				break;
			}
		}
        else {
			color[0] = (char)0;
			color[1] = (char)0;
			color[2] = (char)0;
		}
	}


	//Should return triplet of color code when finished
	void getColor(char* color) {
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
		colorOrbite1(color,i,a,b);
	}
};


class MandelbrotSet {
	long double _xmin;
	long double _xmax;
	long double _ymin;
	long double _ymax;
	int _width;
	int _height;
	unsigned char* _data;

public:

    //x,y coordinate of center of image, width and height the number of pixels to compute
    MandelbrotSet(long double xmin,long double xmax,long double ymin,long double ymax,int width=640,int height=480) {
        this->_xmin = xmin;
        this->_xmax = xmax;
        this->_ymin = ymin;
        this->_ymax = ymax;
        this->_width = width;
        this->_height = height;
        this->_data = new unsigned char[height*width*BYTES_PER_PIXEL];
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
                char color[3];
				m.getColor(color);
				this->_data[i*this->_height*BYTES_PER_PIXEL+j*BYTES_PER_PIXEL] = color[0];
				this->_data[i*this->_height*BYTES_PER_PIXEL+j*BYTES_PER_PIXEL+1] = color[1];
				this->_data[i*this->_height*BYTES_PER_PIXEL+j*BYTES_PER_PIXEL+2] = color[2];
                j+=1;
			}
            i+=1;
		}
	}

    void exportToImage(std::string filename) {
		generateBitmapImage((unsigned char*) this->_data, this->_height, this->_width, filename.c_str());
        //return Image.fromarray(this->_data)
	}

};

////////////////////////////////////////////////////////
// Main entry point
////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
   int counter;

    if(argc!=8) {
        std::cerr << "Usage: " << argv[0] << " Xmin Xmax Ymin Ymax Pixel_Width Pixel_Height Output_Filename" << std::endl;
        std::cerr << "Example: " << argv[0] << " -2 1 -1.5 1.5 1000 1000 mset.bmp" << std::endl;
	}
    if(argc==8)
    {
		long double xmin = stold(std::string(argv[1]));
		long double xmax = stold(std::string(argv[2]));
		long double ymin = stold(std::string(argv[3]));
		long double ymax = stold(std::string(argv[4]));
		int width = stoi(std::string(argv[5]));
		int height = stoi(std::string(argv[6]));
		std::string filename = std::string(argv[7]);

 		Timer t = Timer();
 		t.start();
 		//MandelbrotSet mSet = MandelbrotSet(-2,1,-1.5,1.5,1000,1000);
		MandelbrotSet mSet = MandelbrotSet(xmin,xmax,ymin,ymax,width,height);
 		mSet.compute();
 		t.stop();
		mSet.exportToImage(filename);
 		double e = t.getElapsed_time();
 		std::cout << "Mandelbrot Set generated in " << e << " seconds" << std::endl; 
		std::cout << "Area:       xmin=" << xmin << " xmax=" << xmax << std::endl; 
		std::cout << "            ymin=" << ymin << " ymax=" << ymax << std::endl; 
		std::cout << "Resolution: width=" << width << " height=" << height << std::endl; 
		std::cout << "Output:     file=" << filename << std::endl; 
    }

    return 0;
}

// int main(int argc, char *argv[]) {
// 	int NB_TEST = 10;
//     double total_time = 0;
// 	for (size_t n = 0; n < NB_TEST; n++)
// 	{	
// 		Timer t = Timer();
// 		t.start();
// 		MandelbrotSet mSet = MandelbrotSet(-2,1,-1.5,1.5,1000,10000);
// 		mSet.compute();
// 		t.stop();
// 		double e = t.getElapsed_time();
// 		std::cout << "Boucle " << n << " ran in " << e << "seconds" << std::endl;
// 		total_time += e;
// 	}
// 	std::cout << "Have run  " << NB_TEST << " Mandelbrot Set calculations in " << total_time/NB_TEST << " seconds average" << std::endl;
// 	return 0;
// }
