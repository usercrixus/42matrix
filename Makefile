OBJ = \
	srcs/main.o

HDR = \
	srcs/Matrix/Matrix.hpp \
	srcs/Matrix/VectorAlgebra.hpp

all: matrix

matrix: $(OBJ)
	c++ $(OBJ) -o $@

%.o: %.cpp $(HDR)
	c++ -Wall -Wextra -Werror -Werror -std=c++17 -c $< -o $@

clean:
	rm -fr $(OBJ)

fclean: clean
	rm -f matrix

re: fclean all

.PHONY: clean fclean re