CXX = g++
CXXFLAGS = -std=c++14 -O3 -DNODEBUG -W -Wall -Wno-deprecated -ftree-vectorize -march=native
LINKFLAGS = -lm

all: recomp repair

recomp: rlslp.cpp recompression.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) rlslp.cpp recompression.cpp -o recomp

repair: repair.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) repair.cpp -o repair

clean:
	rm -f recomp repair *.o *~


