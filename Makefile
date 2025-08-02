OBJ = \
	srcs/main.o

HDR = \
	srcs/Matrix/Matrix.hpp \
	srcs/Matrix/VectorAlgebra.hpp

all: tester

tester: $(OBJ)
	c++ $(OBJ) -o $@

%.o: %.cpp $(HDR)
	c++ -c $< -o $@