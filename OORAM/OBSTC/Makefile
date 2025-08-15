# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -Wall -O2 -Wno-deprecated-declarations

# Includes and libraries
INCLUDES = -I../AESN -I../OPQ
LIBS = -lssl -lcrypto

# Source files
SRC = main.cpp OORAM.cpp ../../OPQ/OPQ.cpp ../../AESN/AMAES.cpp Utils.cpp
TARGET = test

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean
clean:
	rm -f $(TARGET) *.o ../../AESN/*.o ../../OPQ/*.o

.PHONY: all clean
