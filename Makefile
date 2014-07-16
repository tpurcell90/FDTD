CC = g++

CFLAGS = -O2 -Wall -Wno-sign-compare -Wno-unused-function -Werror -std=c++11 -openmp

OBJ  =   obj/main.o  obj/Inputs.o obj/FDTDField.o

LIBS =  -lpthread -lm

HEADS =  src/source.hpp src/Grid.hpp src/FDTDField.hpp src/Inputs.hpp
BIN  =   FDTD

RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) -I./src $(OBJ) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) -c  src/main.cpp -o obj/main.o -I./src

obj/Inputs.o: src/Inputs.cpp
	$(CC) $(CFLAGS)  -c src/Inputs.cpp  -o obj/Inputs.o  -I./src

obj/FDTDField.o: src/FDTDField.cpp
	$(CC) $(CFLAGS)  -c src/FDTDField.cpp  -o obj/FDTDField.o  -I./src