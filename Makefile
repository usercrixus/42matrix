OBJ = \
	srcs/main.o

HDR = \
	srcs/Matrix/VectorAlgebra.hpp
	srcs/Matrix/Matrix.hpp

all: tester

tester: $(OBJ) $(HDR)
	c++ $(OBJ) -o $@

%.o: %.cpp
	c++ -c $^ -o $@