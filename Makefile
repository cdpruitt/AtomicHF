################################### SETUP ######################################

# Define important locations
BIN = bin/
SOURCE = src/

# Define compilation flags
COMPILER = clang++ 
CFLAGS = -lstdc++ -I$(shell root-config --incdir) -std=c++11 -Wall -Wextra

# For integration with ROOT data analysis framework
LINKOPTION = $(shell root-config --libs)

############################### DEFINE TARGETS #################################

# List all targets
TARGETS = calculateHF

all: $(addprefix $(BIN), $(TARGETS))

# Build driver (main data analysis engine)
CALCULATE_HF_SOURCES = calculateHF.cpp buildHydrogenicWF.cpp mathFunctions.cpp
$(BIN)calculateHF: $(addprefix $(SOURCE), $(CALCULATE_HF_SOURCES))
	$(COMPILER) $(CFLAGS) -o $(BIN)calculateHF $(addprefix $(SOURCE), $(CALCULATE_HF_SOURCES)) $(LINKOPTION)

############################### DEFINE OBJECTS #################################

# Rule for building all targets
%.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@ $(LINKOPTION)

clean:
	rm -f *.o
