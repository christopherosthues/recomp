CXX = g++
CXXFLAGS = -std=c++14 -O3 -DNODEBUG -W -Wall -Wno-deprecated -ftree-vectorize -march=native
LINKFLAGS = -lm

all: recomp

recomp: rlslp.cpp recompression.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) rlslp.cpp recompression.cpp -o recomp

clean:
	rm -f recomp *.o *~


