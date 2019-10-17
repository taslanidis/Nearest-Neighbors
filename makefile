all : lsh hypercube curves

curves : HashTable.o Helper_Functions.o Traversals.o Curves.o
	@echo "Compile curves ...";
	g++ -I ./lib ./build/Curves.o ./build/Helper_Functions.o ./build/Traversals.o -o ./build/curves.x

Curves.o :
	g++ -I ./lib -c ./src/Curves/Curves.cpp -o ./build/Curves.o

Traversals.o :
	g++ -I ./lib -c ./src/Curves/Traversals.cpp -o ./build/Traversals.o

hypercube : HashTable.o Helper_Functions.o LSH_Functions.o BHC_Functions.o BHC.o
	@echo "Compile hypercube ...";
	g++ -I ./lib ./build/BHC.o ./build/BHC_Functions.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/hypercube.x

BHC.o :
	g++ -I ./lib -c ./src/Hypercube/BHC.cpp -o ./build/BHC.o

BHC_Functions.o :
	g++ -I ./lib -c ./src/Hypercube/BHC_Functions.cpp -o ./build/BHC_Functions.o

lsh : HashTable.o Helper_Functions.o LSH_Functions.o LSH.o
	@echo "Compile lsh ...";
	g++ -I ./lib ./build/LSH.o ./build/LSH_Functions.o ./build/Helper_Functions.o ./build/HashTable.o -o ./build/lsh.x

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
