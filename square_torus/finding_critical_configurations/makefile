#------------------------
# makefile for disk_search
#------------------------

#------------------------
# Paths to source files
#------------------------
SRC_DIR       = src
ARMADILLO_DIR = /usr/include/
NAUTY_DIR     = src/libraries/nauty


#------------------------
#  Source code for current project
#------------------------
SOURCES = canon_label.cpp config.cpp energy.cpp \
          hessian.cpp jacobian.cpp main.cpp \
          output.cpp parser.cpp radius.cpp \
          simulation.cpp utility.cpp


#------------------------
#  Source code for external libraries
#------------------------
NAUTY_SRCS = nauty.c nautil.c naugraph.c schreier.c naurng.c


#------------------------
#  Object files for current project
#------------------------
OBJS = $(patsubst %.cpp, $(SRC_DIR)/%.o, $(SOURCES) ) \
       $(patsubst %.c, $(NAUTY_DIR)/%.o, $(NAUTY_SRCS) )

OBJS_PAR = $(patsubst %.cpp, $(SRC_DIR)/%.po, $(SOURCES) ) \
           $(patsubst %.c, $(NAUTY_DIR)/%.po, $(NAUTY_SRCS) )

OBJS_DBG = $(patsubst %.cpp, $(SRC_DIR)/%.do, $(SOURCES) ) \
           $(patsubst %.c, $(NAUTY_DIR)/%.do, $(NAUTY_SRCS) )


#------------------------
#  Compiler options
#------------------------
CC = gcc
#CC  = clang
CXX = g++
#CXX = clang++ -std=c++11

LIBS = -larmadillo

FAST_FLAGS     = -O2 -Wno-write-strings
PARALLEL_FLAGS = -O2 -fopenmp -Wno-write-strings
#PARALLEL_FLAGS = -O2 -fopenmp=libgomp -Wno-write-strings
DEBUG_FLAGS    = -g -Wfatal-errors -Wno-write-strings

INCL_LOCATION_OPT = -I ${ARMADILLO_DIR} -I ${NAUTY_DIR}

EXEC = crit_search
.PHONY: clean debug fast parallel


#---------------------------
# make fast
#---------------------------
fast:$(OBJS)
	$(CXX) ${FAST_FLAGS} $(OBJS) -o $(EXEC) $(LIBS)

#---------------------------
# make parallel
#---------------------------
parallel: $(OBJS_PAR)
	$(CXX) ${PARALLEL_FLAGS} $(OBJS_PAR) -o $(EXEC) $(LIBS)

#---------------------------
# make debug
#---------------------------
debug:$(OBJS_DBG)
	$(CXX) $(DEBUG_FLAGS) $(OBJS_DBG) -o $(EXEC) $(LIBS)


#---------------------------
# make requirements
#---------------------------
%.o : %.cpp
	$(CXX) -c $(FAST_FLAGS) $(INCL_LOCATION_OPT) $< -o $@

%.po : %.cpp
	$(CXX) -c $(PARALLEL_FLAGS) $(INCL_LOCATION_OPT) $< -o $@

%.do : %.cpp
	$(CXX) -c $(DEBUG_FLAGS) $(INCL_LOCATION_OPT) $< -o $@

%.o : %.c
	$(CC) -c $(FAST_FLAGS) $(INCL_LOCATION_OPT) $< -o $@

%.po : %.c
	$(CC) -c $(PARALLEL_FLAGS) $(INCL_LOCATION_OPT) $< -o $@

%.do : %.c
	$(CC) -c $(DEBUG_FLAGS) $(INCL_LOCATION_OPT) $< -o $@


#---------------------------
# make clean
#---------------------------
clean:
	rm -f $(OBJS) $(OBJS_PAR) $(OBJS_DBG) $(EXEC) *~ 
