build: SHELL:=/bin/bash
build:
	-cd bamtools; mkdir build; cd build; cmake ..; make; cd ../lib; rm *.so; rm *.so.*; rm *.dylib;
	cd rMATS_C; make;
	cd rMATS_pipeline; python setup.py build_ext;
	cp `find ./rMATS_pipeline/build | grep so` .;
