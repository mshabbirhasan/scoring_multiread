ifndef CXXFLAGS
  CXXFLAGS = -O3
endif

CXXFLAGS+=-Iinclude

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  LIBS += -lrt
endif

CXX = g++
default: all
all: main

LIB_SRC = $(wildcard src/*.cpp)
LIB_OBJ = $(patsubst %.cpp, %.o, $(LIB_SRC))

-include $(pathsubst %.o, $(LIB_OBJ))

$(OBJS): %.o : %.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $< 

main: $(LIB_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

clean:
	rm -f $(LIB_OBJ) main

.phony: clean default
