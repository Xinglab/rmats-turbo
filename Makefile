build: SHELL:=/bin/bash
build:
	cd bamtools; mkdir -p build; cd build; cmake ..; make;
	-cd bamtools/lib; rm *.so *.so.* *.dylib
	cd rMATS_C; make;
	cd rMATS_pipeline; python setup.py build_ext;
	cp `find ./rMATS_pipeline/build | grep so` .;
