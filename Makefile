.PHONY: all classic quantum evaluation clean

all: classic quantum evaluation plot 

classical/build/classical_eigensolver: $(wildcard classical/src/*) $(wildcard classical/include/*)
	mkdir -p classical/build
	cd classical/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

data/wavefunctions.dat: classical/build/classical_eigensolver
	cd classical/build && ./classical_eigensolver

classic: data/wavefunctions.dat

quantum: 
	cd quantum && python3 vqe.py

evaluation:
	cd quantum && python3 error_estimation.py

plot: 
	cd Plotting/ && make 

clean:
	cd Plotting/ && make clean 
	rm -rf classical/build
	rm -f data/*
	rm -f Plots/*