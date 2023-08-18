

PROGRAM_NAME := md.out

# compiler
CC := gcc
CXX := g++
OPTS := -O3 -flto 
OPTS += -Wall -Wextra # -fopenmp
OPTS += -DDSFMT_MEXP=19937
CCFLAGS  := $(OPTS)
CXXFLAGS := -std=c++11 $(OPTS)

INCS := -I./dSFMT-src-2.2.3 -I$(HOME)/local/include #-I/opt/homebrew/include
LIBS := -L$(HOME)/local/lib #-L/opt/homebrew/lib -I$(HOME)/local/include
DYLIBS := -lyaml-cpp 


# directory
OUT_DIR := build
SOURCE_DIR := ./
HEADER_DIR := ./

SOURCES := $(wildcard $(SOURCE_DIR)/*.cpp) 
HEADERS := $(wildcard $(HEADER_DIR)/*.h)

PROGRAM_DIR := $(OUT_DIR)
OBJ_DIR := $(OUT_DIR)
DEPEND_DIR := $(OUT_DIR)

PROGRAM := $(PROGRAM_DIR)/$(PROGRAM_NAME)
SOURCE_NAMES = $(notdir $(SOURCES))
OBJS := $(addprefix $(OBJ_DIR)/,$(SOURCE_NAMES:.cpp=.o))
DEPENDS := $(addprefix $(DEPEND_DIR)/,$(SOURCE_NAMES:.cpp=.depend))

.PHONY: all
all: $(DEPENDS) $(PROGRAM)
$(PROGRAM): $(OBJS)
	@mkdir -p $(PROGRAM_DIR)
	$(CXX) $(CXXFLAGS) $(INCS) $(LIBS) $^ $(DYLIBS) -o $(PROGRAM)

$(DEPEND_DIR)/%.depend: $(SOURCE_DIR)/%.cpp $(HEADERS)
	@echo "generating $@"
	@mkdir -p $(DEPEND_DIR)
	@$(CXX) $(CXXFLAGS) $(INCS) -I$(HEADER_DIR) -MM $< > $@

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCS) -I$(HEADER_DIR) -c $^ -o $@


ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEPENDS)
endif

.PHONY : clean
clean:
	$(RM) -r $(OUT_DIR)
