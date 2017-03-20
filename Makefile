all: CImg.h a2.cpp
	g++ a2.cpp -o a2 -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lpthread -I. -Isiftpp -O3 siftpp/sift.cpp

clean:
	rm a2