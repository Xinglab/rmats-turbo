build: SHELL:=/bin/bash
build:
	cd bamtools; mkdir -p build; cd build; cmake ..; make;
	# rm -f to ignore nonexistent files since *.dylib will only exist for mac
	cd bamtools/lib; rm -f *.so *.so.* *.dylib
	cd rMATS_C; make;
	cd rMATS_pipeline; python setup.py build_ext;
	cp `find ./rMATS_pipeline/build | grep so` .;
