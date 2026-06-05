TARGET = main

CXX = mpicxx

INC_DIR = include
SRC_DIR = src

CXXFLAGS = -std=c++20 -Wall -Wextra -O3 -I$(INC_DIR) -fopenmp

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

DEPS = $(wildcard $(INC_DIR)/*.hpp)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)


.PHONY: all clean

.PHONY: doc
doc:
	doxygen Doxyfile


