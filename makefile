ph_cov: ph_cov.cpp
	g++ ph_cov.cpp -o ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./ph_cov.x
	python3 plot_ph_cov.py