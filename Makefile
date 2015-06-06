Inc = -I/opt/local/include/
#Inc = -I/projects/caldwelt
#Inc2 = -I/curc/tools/x_86_64/rh6/boost/1.56/anaconda/2.0.0/openmpi/1.8.2/intel/13.0.0/include

all: 
	g++ -c -Os -std=c++11 $(Inc) main.cpp global.cpp interp.cpp system.cpp pendcart_3link.cpp treevertex.cpp tree.cpp
	g++ -o agilerrt $(Inc) main.o global.o interp.o system.o pendcart_3link.o treevertex.o tree.o -lgsl -lcblas -lm
#       icpc -c -Os -std=c++11 $(Inc) $(Inc) main.cpp global.cpp interp.cpp system.cpp pendcart_3link.cpp treevertex.cpp tree\
.cpp
#       icpc -o agilerrt $(Inc) $(Inc2) main.o global.o interp.o system.o pendcart_3link.o treevertex.o tree.o -lgsl -lcblas -\
lm

clean:
	rm -rf *o main.o global.o interp.o system.o pendcart_3link.o treevertex.o
