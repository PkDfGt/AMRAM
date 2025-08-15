CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS = -lssl -lcrypto

# Source files
SRCS = main.cpp ../AESN/AMAES.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = test

# Build target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Clean target
clean:
	rm -f $(TARGET) $(OBJS)
