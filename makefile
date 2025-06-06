# Compiler
CXX = g++
EIGENDIR = /opt/eigen-3.40/
# Compiler flags
CXXFLAGS = -I ./ -std=c++11 -Wall -I /usr/include/suitesparse -O3 -I $(EIGENDIR)
LDFLAGS = -lumfpack -lcholmod -lamd -lsuitesparseconfig
#-ffast-math -fno-finite-math-only 
# Source and object files
SRC_DIR = src
OBJ_DIR = obj
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Output executable
TARGET = simplex

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean target to remove object files and the executable
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

# Phony targets to avoid conflicts with files named 'clean', 'all', etc.
.PHONY: all clean

