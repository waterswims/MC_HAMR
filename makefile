# Makefile for project

#################################################################
## Location of the HDF5 instalation
#################################################################
## IRIDIS BUILD - Using intel compilers
#HDFLIB=/local/software/gdf5/1.8.11/intel-par/lib
#HDFINC=/local/software/gdf5/1.8.11/intel-par/include

#################################################################
## Directories
#################################################################

INC_PATH = includes
LIB_PATH = lib
OBJ_PATH = obj
TEST_PATH = tests

#################################################################
## Files to use
#################################################################

TEST_FILES = $(wildcard $(TEST_PATH)/*.hpp)
NOMAIN_FILES = $(filter-out $(LIB_PATH)/main.cpp, $(wildcard $(LIB_PATH)/*.cpp))
## NON-INTEL
SOURCE_FILES = $(filter-out $(LIB_PATH)/mklrand.cpp, $(NOMAIN_FILES))
## INTEL
# SOURCE_FILES = $(filter-out $(LIB_PATH)/stdrand.cpp, $(NOMAIN_FILES))
OBJS = $(addprefix $(OBJ_PATH)/, $(notdir $(SOURCE_FILES:.cpp=.o)))

#################################################################
## Compile options
#################################################################
## INTEL
# CPPFLAGS = -std=c++11 -Ofast -qopenmp -DMKL_ILP64 -I${MKLROOT}/include -I${HDFINC} -ipo
## NON-INTEL
CPPFLAGS = -std=c++11 -Ofast -fopenmp -I${HDFINC} -Wno-narrowing -flto

#################################################################
## Link options
#################################################################
## INTEL
# LDFLAGS = -L${HDFLIB}  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lhdf5 -lz -ipo
## NON-INTEL
LDFLAGS = -L${HDFLIB} -lhdf5 -lz -flto

#################################################################
## Gtest options
#################################################################
GTEST_DIR=$(TEST_PATH)/googletest/googletest
GTEST_FLAGS=-isystem $(GTEST_DIR)/include
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(wildcard $(GTEST_DIR)/src/*.cc) $(wildcard $(GTEST_DIR)/src/*.h) $(GTEST_HEADERS)

#################################################################
## Build scripts
#################################################################

default: $(OBJS) $(OBJ_PATH)/main.o
	$(CC) $(OBJS) $(OBJ_PATH)/main.o -o run $(LDFLAGS)

$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

$(TEST_PATH)/test.o: $(TEST_FILES) $(TEST_PATH)/test.cpp $(TEST_PATH)/libs/gtest_main.a
	$(CC) tests/test.cpp $(GTEST_FLAGS) $(CPPFLAGS) -c -o $(TEST_PATH)/test.o

test: $(OBJS) $(TEST_PATH)/test.o $(TEST_PATH)/libs/gtest_main.a
	$(CC) $(OBJS) $(TEST_PATH)/test.o $(TEST_PATH)/libs/gtest_main.a $(GTEST_FLAGS) $(LDFLAGS) -o test

$(TEST_PATH)/libs/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) 	$(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ \
		$(GTEST_DIR)/src/gtest-all.cc

$(TEST_PATH)/libs/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) 	$(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ \
		$(GTEST_DIR)/src/gtest_main.cc

$(TEST_PATH)/libs/gtest.a : $(TEST_PATH)/libs/gtest-all.o
	$(AR) 	$(ARFLAGS) $@ $^

$(TEST_PATH)/libs/gtest_main.a : $(TEST_PATH)/libs/gtest-all.o $(TEST_PATH)/libs/gtest_main.o
	$(AR) 	$(ARFLAGS) $@ $^

clean:
	-rm -f $(OBJS)
	-rm -f $(OBJ_PATH)/main.o
	-rm -f run
	-rm -f $(TEST_PATH)/test.o
	-rm -f $(TEST_PATH)/libs/*
	-rm -f test
