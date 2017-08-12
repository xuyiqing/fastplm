xcode:
	mkdir build; cd build && cmake -G Xcode ..

clean:
	rm -rf build
	find . -name "*.so" -type f -delete
	find . -name "*.o"  -type f -delete

tar:
	tar cfv fastplm.tar \
		DESCRIPTION NAMESPACE LICENSE README.md \
		src/FixedEffect.h src/FixedEffect.cpp \
		src/PlainModel.h src/PlainModel.cpp \
		src/CrashQueue.h src/CrashQueue.cpp \
		src/FastFESolver.h src/FastFESolver.cpp \
		src/SlowFESolver.h src/SlowFESolver.cpp \
		src/GPSolver.h src/GPSolver.cpp \
		src/MatrixReader.h src/MatrixReader.cpp \
		src/Common.h src/Wrapper.cpp \
		man/solveFE.rd \
		tests/FEDataGenerator.r tests/GPDataGenerator.r \
		tests/SpeedResult.txt tests/SpeedTest.r \
		mathematics/GPanel.pdf mathematics/GPanel.tex \
		mathematics/GradientDescent.pdf mathematics/GradientDescent.tex \
		fastplm.Rproj

