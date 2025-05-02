ph_cov: ph_cov.cpp
	g++ ph_cov.cpp -o ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./ph_cov.x
	python3 plot_ph_cov.py

qubit: qubit.cpp
	g++ qubit.cpp -o qubit.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./qubit.x
	python3 plot_qubit.py