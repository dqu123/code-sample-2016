###############################################################################
# CS/CNS 171 Fall 2016
#
# David Qu
###############################################################################
CC = g++
FLAGS = -g -std=c++11

INCLUDE = -I./ -I/usr/X11R6/include -I/usr/include/GL -I/usr/include
LIBDIR = -L/usr/X11R6/lib -L/usr/local/lib
MAIN_FILES = main.cpp test.cpp
SOURCES = $(filter-out $(MAIN_FILES), $(wildcard *.h *.cpp))
LIBS = -lGLEW -lGL -lGLU -lglut -lm

EXENAME = opengl_renderer

exe: $(EXENAME)

all: $(EXENAME) test

$(EXENAME): $(SOURCES) main.cpp
	$(CC) $(FLAGS) -o $(EXENAME) $(INCLUDE) $(LIBDIR) $(SOURCES) main.cpp $(LIBS)

test: $(SOURCES) test.cpp
	$(CC) $(FLAGS) -o test $(INCLUDE) $(LIBDIR) $(SOURCES) test.cpp $(LIBS)

check: test
	./test

clean:
	rm -f *.o *~ test $(EXENAME)

.PHONY: all clean check

