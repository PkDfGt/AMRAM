CXX := g++
CXXFLAGS := -std=c++17 -fopenmp -Wall -Wextra -Wno-deprecated-declarations
INCLUDES := -I. -I../AESN
LIBS := -lcrypto

SRC := main.cpp TC.cpp ../AESN/AMAES.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := TC

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ)
