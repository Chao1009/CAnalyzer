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
CXXFLAGS      = -shared -pipe -std=c++11 -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -O2 -std=gnu++0x -Wall -W -D_REENTRANT -fPIC $(DEFINES)
INCPATH       = -Iinclude
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
LIBS          = $(SUBLIBS)
# -L$(ROOTSYS)/lib -lCore -lRint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui -lSpectrum -lGL -lpthread 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Files

TARGET        = libCAna.so

OBJECTS_DIR   = obj

SOURCES_DIR   = src

CXX_EXT       = cpp

CXX_OBJECTS   = $(OBJECTS_DIR)/ConfigParser.o \
                $(OBJECTS_DIR)/CAnalyzer.o \
                $(OBJECTS_DIR)/CMatrix.o

first: all
####### Build rules

$(TARGET):  $(CXX_OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(CXX_OBJECTS) $(OBJCOMP) $(LIBS)

all: Makefile $(TARGET)

clean: cleanobj cleantgt

cleanobj:
	-$(DEL_FILE) $(CXX_OBJECTS)

cleantgt:
	-$(DEL_FILE) $(TARGET)

$(OBJECTS_DIR)/%.o: $(SOURCES_DIR)/%.$(CXX_EXT)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

####### Install

#install:   FORCE

#uninstall:   FORCE

#FORCE:


