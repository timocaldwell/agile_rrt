# May need to change include directory
Inc = -I/opt/local/include/

all: 
	g++ -c -Os -std=c++11 $(Inc) main.cpp global.cpp interp.cpp system.cpp pendcart_3link.cpp treevertex.cpp tree.cpp
	g++ -o agilerrt $(Inc) main.o global.o interp.o system.o pendcart_3link.o treevertex.o tree.o -lgsl -lcblas -lm
# may need -lgslcblas

clean:
	rm -rf *o main.o global.o interp.o system.o pendcart_3link.o treevertex.o
