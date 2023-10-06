CXXFLAGS=-DNDEBUG -Wall -O3 -std=c++17 -I. -I./KMC
# -Wno-sign-compare
# -Wno-unused-variable
# -Wno-char-subscripts
# -fopenmp
LIBS=-lz

.PHONY: all

all: main

main: main.o ./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
	rm -rf KMC/kmc_api/*.o
