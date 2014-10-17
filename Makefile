CXX = g++
GTEST_DIR = ../gtest-1.6.0
CPPFLAGS += -I$(GTEST_DIR)/include
# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra
#
# # All tests produced by this Makefile.  Remember to add new tests you
# # created to the list.
TESTS = indexmap_test spanningtree_test
#
ARFLAGS =-q -s
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

all : bdmst $(TESTS)
gtest-all.o : $(GTEST_SRCS_)
	        $(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		              $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	        $(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		              $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	        $(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	        $(AR) $(ARFLAGS) $@ $^

mersenne.o : randomc.h mersenne.cpp
	g++ $(CXXFLAGS) -c mersenne.cpp
graph.o : graph.h graph.cpp
	g++ $(CXXFLAGS) -c graph.cpp
indexset.o : indexset.h indexset.cpp indexmap.h
	g++ $(CXXFLAGS) -c indexset.cpp
indexmap.o : indexmap.h indexmap.cpp
	g++ $(CXXFLAGS) -c indexmap.cpp
spanningtree.o : spanningtree.h spanningtree.cpp indexset.h
	g++ $(CXXFLAGS) -c spanningtree.cpp
ant.o : ant.h ant.cpp graph.h
	g++ $(CXXFLAGS) -c ant.cpp
queue.o : queue.h queue.cpp
	g++ $(CXXFLAGS) -c queue.cpp
bdmst.o : bdmst.cpp ant.h randomc.h graph.h indexset.h parameters.h helper.h spanningtree.h
	g++ $(CXXFLAGS) -c bdmst.cpp
bdmst : bdmst.o graph.o indexset.o indexmap.o spanningtree.o ant.o mersenne.o queue.o
	g++ $(CXXFLAGS) -o bdmst bdmst.o graph.o indexset.o indexmap.o spanningtree.o ant.o mersenne.o queue.o -lgflags

indexmap_test.o : indexmap_test.cpp indexmap.h
	g++  $(CPPFLAGS) $(CXXFLAGS) -c indexmap_test.cpp 
indexmap_test : indexmap.o indexmap_test.o gtest_main.a
	g++  $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ -lpthread
spanningtree_test.o : spanningtree_test.cpp spanningtree.h
	g++  $(CPPFLAGS) $(CXXFLAGS) -c spanningtree_test.cpp 
spanningtree_test : spanningtree.o spanningtree_test.o indexset.o indexmap.o graph.o gtest_main.a
	g++  $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ -lpthread
clean : 
	rm *.o bdmst $(TESTS) gtest_main.a
