all:
	g++ -c -fPIC warpmodule.cpp -o warp_linux.o -I/usr/include/python3.8
	g++ -shared -o warp.so warp_linux.o
	nm -D ./warp.so
	python3 ./use_warp.py
clean:
