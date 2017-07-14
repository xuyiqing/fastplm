xcode:
	mkdir build; cd build && cmake -G Xcode ..

clean:
	rm -rf build
	find . -name "*.so" -type f -delete
	find . -name "*.o"  -type f -delete

tar:
	tar cfv fastplm.tar \
		DESCRIPTION NAMESPACE \
		src/FixedEffect.h src/PlainModel.h \
		src/FastFESolver.h src/FastFESolver.cpp \
		src/SlowFESolver.h src/SlowFESolver.cpp \
		src/GPSolver.h src/GPSolver.cpp \
		src/MatrixReader.h src/MatrixReader.cpp \
		src/Common.h src/Wrapper.cpp
