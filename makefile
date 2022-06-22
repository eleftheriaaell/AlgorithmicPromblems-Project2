executable1 = search
executable2 = cluster
executable3 = test

all : search cluster test

search : data.cpp Frechet_LSH.cpp LSH.cpp hypercube.cpp search_main.cc continuous/*.cpp
	g++ -o search data.cpp Frechet_LSH.cpp LSH.cpp hypercube.cpp search_main.cc continuous/*.cpp -O2

cluster: data.cpp LSH.cpp hypercube.cpp clustering.cpp Frechet_LSH.cpp Frechet_clustering.cpp cluster_main.cc continuous/*.cpp
	g++ -o cluster data.cpp LSH.cpp hypercube.cpp clustering.cpp Frechet_LSH.cpp Frechet_clustering.cpp cluster_main.cc continuous/*.cpp -O2

test: test_main.cc LSH.cpp data.cpp continuous/*.cpp
	g++ -o test test_main.cc LSH.cpp data.cpp continuous/*.cpp -lcunit

clean:
	rm -f $(executable1) $(executable2) $(executable3)