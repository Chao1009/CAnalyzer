#############################################################################
# Makefile for building: PRadEventViewer
# Generated by qmake (3.0) (Qt 5.6.0)
# Project:  QMakefile.pro
# Template: app
# Command: /usr/bin/qmake-qt5 -o Makefile QMakefile.pro
#############################################################################

MAKEFILE      = Makefile
####### Compiler, tools and options

CC            = gcc
CXX           = g++
FORTRAN       = gfortran
CXXFLAGS      = -shared -pipe -std=c++11 -O2 -g -Wall -m64 -mtune=generic -fPIC
FFLAGS        = -fPIC -ffixed-line-length-none
INCPATH       = -Iinclude -I$(ROOTSYS)/include
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
LINK          = g++
LFLAGS        = -shared -Wl,-O1 -Wl,-z,relro
LIBS          = 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 


####### Files
TARGET        = libCAna.so

OBJECTS_DIR   = obj
SOURCES_DIR   = src

CXX_SUFFIX    = cpp
CXX_SOURCES   = canalib \
                ConfigParser \
                ConfigValue \
                ConfigObject \
                CExpData \
                CRadCorr \
                CAnalyzer \
                CMatrix \
                CEstimator \
                CElasTails \
                CHe3Elas \

F_SUFFIX      = f
F_SOURCES     = BostedFit \
                QFSFit \
                rtails \
				sub_poltail \
				sub_math


CXX_OBJECTS   = $(addprefix $(OBJECTS_DIR)/, $(CXX_SOURCES:=.cpp.o))
F_OBJECTS     = $(addprefix $(OBJECTS_DIR)/, $(F_SOURCES:=.f.o))
OBJECTS       = $(CXX_OBJECTS) $(F_OBJECTS)

####### Build rules
first: all

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

all: Makefile $(TARGET)

clean: cleanobj cleantgt

cleanobj:
	-$(DEL_FILE) $(OBJECTS)

cleantgt:
	-$(DEL_FILE) $(TARGET)

$(OBJECTS_DIR)/%.cpp.o: $(SOURCES_DIR)/%.$(CXX_SUFFIX)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

$(OBJECTS_DIR)/%.f.o: $(SOURCES_DIR)/%.$(F_SUFFIX)
	$(FORTRAN) -c $(FFLAGS) $(INCPATH) -o $@ $<

####### Install

#install:   FORCE

#uninstall:   FORCE

#FORCE:


