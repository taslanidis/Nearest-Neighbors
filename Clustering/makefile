all : cluster

cluster : Initializers.o Helper_Functions.o K_Means.o Clustering.o
	@echo "Compile Cluster ...";
	g++ -std=c++11 -I ./lib ./build/Initializers.o ./build/K_Means.o ./build/Clustering.o ./build/Helper_Functions.o -o ./build/cluster.x

Clustering.o :
	g++ -std=c++11 -I ./lib -c ./src/Clustering.cpp -o ./build/Clustering.o

K_Means.o :
	g++ -std=c++11 -I ./lib -c ./src/K_Means.cpp -o ./build/K_Means.o

Initializers.o :
	g++ -std=c++11 -I ./lib -c ./src/Initializers/Initializers.cpp -o ./build/Initializers.o

Helper_Functions.o :
	g++ -std=c++11 -I ./lib -c ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x