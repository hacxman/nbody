all:
	echo "With O3"
	g++ nb.cc -std=c++11 -O3 -lboost_timer -lboost_system
	./a.out
	echo "With O2"
	g++ nb.cc -std=c++11 -O2 -lboost_timer -lboost_system
	./a.out
	echo "Noopt"
	g++ nb.cc -std=c++11 -lboost_timer -lboost_system
	./a.out
