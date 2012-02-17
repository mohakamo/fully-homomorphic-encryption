NTL = /Users/cs1/sw

VPATH = $(NTL)/lib:$(NTL)/include
objects = main.o R_Ring_Vector.o FHE.o

arithmetic_objects = R_Ring_Number.h R_Ring_Vector.h R_Ring_Matrix.h

SRCPATH = ./
LIBDIRS = -L$(NTL)/lib
INCDIRS = -I$(NTL)/include

LIBS = -lntl

CC = g++
CXXFLAGS = $(COMPILERFLAGS) -g $(INCDIRS) $(LIBDIRS)

chehoff : $(objects)
	$(CC) $(CXXFLAGS) -o chehoff $(objects) $(LIBS)

main.o : $(arithmetic_objects) Pair.h GLWE.h FHE.h LSS.h Tests.h NTL/ZZ.h main.cpp
	$(CC) $(CXXFLAGS) -c main.cpp
R_Ring_Vector.o : $(arithmetic_objects) NTL/ZZ.h R_Ring_Vector.cpp
	$(CC) $(CXXFLAGS) -c R_Ring_Vector.cpp
FHE.o : $(arithmetic_objects) Pair.h GLWE.h FHE.h NTL/ZZ.h FHE.cpp
	$(CC) $(CXXFLAGS) -c FHE.cpp

.PHONY : clean
clean :
	rm chehoff $(objects)