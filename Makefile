CXX = g++
CXXFLAGS = -std=c++1z -Ofast -DNODEBUG -W -Wall -Wno-deprecated
LINKFLAGS = -lm

all: recompression_pos

recompression_pos: recompression.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) recompression.cpp -o recompression_pos

clean:
	rm -f recompression_pos *.o *~


