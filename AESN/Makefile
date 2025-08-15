# 
CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -Wall -O2

#
INCLUDES = -I../AESN
LIBS = -lssl -lcrypto

# 
SRC = main.cpp ../AESN/AMAES.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = test

# Target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Clean
clean:
	rm -f $(TARGET) *.o ../AESN/*.o

.PHONY: all clean
