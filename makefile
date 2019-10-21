all : lsh hypercube curves_grid_lsh curves_projection_lsh

curves_projection_lsh : HashTable.o Helper_Functions.o Traversals.o Projections_Curves.o
	@echo "Compile projections curves ...";
	g++ -I ./lib ./build/LSH.o ./build/Projections_Curves.o ./build/Helper_Functions.o ./build/Traversals.o -o ./build/curves.x

Projections_Curves.o :
	g++ -I ./lib -c ./src/Projections_Curves/RandomProjections.cpp -o ./build/RandomProjections.o

curves_grid_lsh : HashTable.o Helper_Functions.o Traversals.o Grid_Curves.o
	@echo "Compile grid curves ...";
	g++ -I ./lib ./build/Grid_Curves.o ./build/LSH.o ./build/Helper_Functions.o ./build/Traversals.o -o ./build/curves_grid_lsh.x

Grid_Curves.o :
	g++ -I ./lib -c ./src/Grid_Curves/Curves.cpp -o ./build/Grid_Curves.o

Traversals.o :
	g++ -I ./lib -c ./src/Grid_Curves/Traversals.cpp -o ./build/Traversals.o

hypercube : HashTable.o Helper_Functions.o LSH_Functions.o BHC_Functions.o BHC.o BHC_main.o
	@echo "Compile hypercube ...";
	g++ -I ./lib ./build/BHC_main.o ./build/BHC.o ./build/BHC_Functions.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/hypercube.x

BHC_main.o :
	g++ -I ./lib -c ./src/Hypercube/BHC_main.cpp -o ./build/BHC_main.o

BHC.o :
	g++ -I ./lib -c ./src/Hypercube/BHC.cpp -o ./build/BHC.o

BHC_Functions.o :
	g++ -I ./lib -c ./src/Hypercube/BHC_Functions.cpp -o ./build/BHC_Functions.o

lsh : HashTable.o Helper_Functions.o LSH_Functions.o LSH.o LSH_main.o
	@echo "Compile lsh ...";
	g++ -I ./lib ./build/LSH_main.o ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/lsh.x

LSH_main.o :
	g++ -I ./lib -c ./src/LSH/LSH_main.cpp -o ./build/LSH_main.o

LSH.o :
	g++ -I ./lib -c ./src/LSH/LSH.cpp -o ./build/LSH.o

LSH_Functions.o :
	g++ -I ./lib -c ./src/LSH/LSH_Functions.cpp -o ./build/LSH_Functions.o

Helper_Functions.o :
	g++ -I ./lib -c ./src/Helper_Functions.cpp -o ./build/Helper_Functions.o

HashTable.o :
	g++ -I ./lib -c ./src/HashTable.cpp -o ./build/HashTable.o

clean:
	-rm -f ./build/*.o
	-rm -f ./build/*.x
