#################################################################
########################   MakeVars   ###########################
#################################################################

#CXX=g++

LIB=../lib
HOME=./src
INCLUDE=./include
OBJ_DIR=./obj
EIGEN_DIR=../../eigen-3.4.0
XTC_DIR=../../xdrfile-1.1.4/install/include/xdrfile/
XTC_LIB=../../xdrfile-1.1.4/install/lib

NCXX_PREFIX=../../netcdf-cxx4-4.3.1/install
NCC_PREFIX=../../netcdf-c-4.9.2/install

NCXX_INCLUDE = $(NCXX_PREFIX)/include
NCXX_LIB     = $(NCXX_PREFIX)/lib64

NCC_INCLUDE  = $(NCC_PREFIX)/include
NCC_LIB      = $(NCC_PREFIX)/lib64

LDLIBS = -lnetcdf-cxx4 -lnetcdf -lstdc++fs -lxdrfile

CXX_OPT= -std=c++17 -I "$(INCLUDE)" \
          -I "$(EIGEN_DIR)"         \
          -I "$(LIB)"               \
          -I "$(XTC_DIR)"           \
          -L "$(XTC_LIB)"           \
          -I$(NCXX_INCLUDE)         \
          -I$(NCC_INCLUDE)          \
          -L$(NCXX_LIB)             \
          -L$(NCC_LIB)              \
          -Wall -Wextra -O2

MKDIR_OBJ=mkdir -p $(OBJ_DIR)

TARGET=xelym.x

SRC=$(HOME)/main.cpp         \
    $(HOME)/analyze.cpp      \
    $(HOME)/input.cpp        \
    $(LIB)/read_file.cpp     \
    $(LIB)/dcd.cpp           \
    $(LIB)/dcd_r.cpp         \
    $(LIB)/traj_tools.cpp    \
    $(LIB)/math_tools.cpp    \
    $(LIB)/pdb_tools.cpp     \
    $(LIB)/psf_tools.cpp     \
    $(LIB)/pos_tools.cpp     \
    $(LIB)/output.cpp        \
    $(LIB)/itp_tools.cpp     \
    $(LIB)/xtc_tools.cpp     \
    $(LIB)/extension.cpp     \
    $(LIB)/nc_tools.cpp      \
    $(LIB)/make_hist.cpp

OBJ=$(OBJ_DIR)/main.o          \
    $(OBJ_DIR)/analyze.o       \
    $(OBJ_DIR)/input.o         \
    $(OBJ_DIR)/read_file.o     \
    $(OBJ_DIR)/dcd.o           \
    $(OBJ_DIR)/dcd_r.o         \
    $(OBJ_DIR)/traj_tools.o    \
    $(OBJ_DIR)/math_tools.o    \
    $(OBJ_DIR)/pdb_tools.o     \
    $(OBJ_DIR)/psf_tools.o     \
    $(OBJ_DIR)/pos_tools.o     \
    $(OBJ_DIR)/output.o        \
    $(OBJ_DIR)/itp_tools.o     \
    $(OBJ_DIR)/xtc_tools.o     \
    $(OBJ_DIR)/extension.o     \
    $(OBJ_DIR)/nc_tools.o      \
    $(OBJ_DIR)/make_hist.o

#################################################################
########################   Makefile   ###########################
#################################################################

all: $(TARGET)
	@echo "Compilation Success"

$(TARGET): $(OBJ)
	$(CXX) $(CXX_OPT) $(OBJ) -o $@ $(LDLIBS)

$(OBJ_DIR)/%.o: $(HOME)/%.cpp
	@$(MKDIR_OBJ)
	$(CXX) $(CXX_OPT) -c $< -o $@

$(OBJ_DIR)/%.o: $(LIB)/%.cpp
	@$(MKDIR_OBJ)
	$(CXX) $(CXX_OPT) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ_DIR)/*.o
