xcode:
	mkdir build; cd build && cmake -G Xcode ..

clean:
	rm -rf build
	find . -name "*.so" -type f -delete
	find . -name "*.o"  -type f -delete

