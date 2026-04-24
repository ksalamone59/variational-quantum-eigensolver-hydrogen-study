.PHONY: all classic clean

all: classic quantum evaluation plot 

classical/build/classical_eigensolver: $(wildcard classical/src/*) $(wildcard classical/include/*)
	mkdir -p classical/build
	cd classical/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

data/wavefunctions.dat: classical/build/classical_eigensolver
	cd classical/build && ./classical_eigensolver

classic: data/wavefunctions.dat

quantum: quantum/.vqe.done

quantum/.vqe.done: data/wavefunctions.dat quantum/vqe.py
	cd quantum && python3 vqe.py
	touch quantum/.vqe.done

evaluation: Plots/error_heatmap.png

Plots/error_heatmap.png: quantum/.vqe.done quantum/error_estimation.py
	cd quantum && python3 error_estimation.py

plot: 
	$(MAKE) -C Plotting/
	
clean:
	cd Plotting/ && make clean 
	rm -rf classical/build
	rm -f data/*
	rm -f Plots/*
	rm -f quantum/.vqe.done