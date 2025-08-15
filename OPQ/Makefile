# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -Wall -O2 -Wno-deprecated-declarations

# Includes and libraries
INCLUDES = -I../AESN
LIBS = -lssl -lcrypto

# Source files
SRC = main.cpp OPQ.cpp PQ.cpp ../AESN/AMAES.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = test

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean
clean:
	rm -f $(TARGET) *.o ../AESN/*.o

.PHONY: all clean
