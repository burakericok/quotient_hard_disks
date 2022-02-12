#------------------------
# makefile for Tabu search
#------------------------
SRC_DIR = src
INCL_DIR_OPT = 

# LIBS = -lm -lopenblas
LIBS = -lopenblas -larmadillo 
LIBS_DIR_OPT = 

EXECUTABLE = calc_pair_dists
#LIBRARY_DIR = /usr/local/lib
#HEADERS_DIR = /usr/local/include/sbdesc
#EXAMPLE = sb_example

#------------------------
# Source code for project
#------------------------
SRCS = main.cpp generateNeighbors.cpp utility.cpp distances.cpp tabuSearch.cpp 
# EXAMPLE_SRCS = 

#------------------------
# Compiler and linker options
#------------------------
CC   	  = g++
FLAGS     = -fopenmp -DARMA_DONT_USE_WRAPPER
DBG_FLAGS = -g -Wall -Wextra -Wfatal-errors
FLAGS_OPT = 

OBJS = $(patsubst %.cpp, $(SRC_DIR)/%.po, $(SRCS))
DBG_OBJS = $(patsubst %.cpp, $(SRC_DIR)/%.do, $(SRCS))
# EXAMPLE_OBJS = $(patsubst %.cpp, $(SRC_DIR)/%.po, $(EXAMPLE_SRCS))

# HEADERS = $(wildcard $(SRC_DIR)/sb_*.h)

# LIBS_DIR_OPT += -L$(LIBRARY_DIR)
# LIB_PATH = $(LIBRARY_DIR)/$(LIBRARY)

.PHONY: default debug clean

#---------------------------
# make
#---------------------------
# default : FLAGS_OPT += -fPIC
default : $(OBJS)
	$(CC) -o $(EXECUTABLE) $(OBJS) $(LIBS) $(FLAGS)

#---------------------------
# make debug
#---------------------------
# debug : FLAGS_OPT += -fPIC
debug : $(DBG_OBJS)
	$(CC) -o $(EXECUTABLE) $(DBG_OBJS) $(LIBS)

#---------------------------
# make requirements
#---------------------------
%.po : %.cpp
	$(CC) -c $(FLAGS) $(FLAGS_OPT) $(INCL_DIR_OPT) $< -o $@

%.do : %.cpp
	$(CC) -c $(FLAGS) $(DBG_FLAGS) $(FLAGS_OPT) $(INCL_DIR_OPT) $< -o $@

#---------------------------
# make clean
#---------------------------
clean :
	rm -rf $(OBJS) $(DBG_OBJS) $(EXECUTABLE)

