
ifeq ($(OS),Windows_NT)
	LBITS = $(shell echo %PROCESSOR_ARCHITECTURE%)
	LBITS64 = $(shell echo %PROCESSOR_ARCHITEW6432%)

	ifeq ($(LBITS),AMD64)
		GLUT_DLL_PATH = ./freeglut/bin/x64
		GLUT_DLL = ${GLUT_DLL_PATH}/freeglut.dl
	else ifeq ($(LBITS64),AMD64)
		GLUT_DLL_PATH = ./freeglut/bin/x64
		GLUT_DLL = ${GLUT_DLL_PATH}/freeglut.dl
	else
		GLUT_DLL_PATH = ./freeglut/bin
		GLUT_DLL = ${GLUT_DLL_PATH}/freeglut.dll
	endif

	INCLUDES = -I. -I./include -I./freeglut/include
	LIBS = -L./lib -lm -L${GLUT_DLL_PATH} -lglu32 -lopengl32 -lfreeglut
else
	UNAME = ${shell uname}

	ifeq ($(UNAME),Linux)
		INCLUDES = -I. -I./include
		LIBS = -L./lib -lm -lglut -lGLU -lGL
	else ifeq ($(UNAME),Darwin)
		INCLUDES = -I. -I./include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers/
		LIBS = -L./lib -lm -framework OpenGL -framework GLUT
	endif
endif

NVCC = nvcc
NVCCFLAGS = -DFOOD_ANGLE=${ANGLE} -DFOOD_DIST=${DIST} -D_FORCE_INLINES

CC = cc
CFLAGS = -O3 -DFOOD_ANGLE=${ANGLE} -DFOOD_DIST=${DIST} -D_FORCE_INLINES

CXX = g++
CXXFLAGS = -O3 -DFOOD_ANGLE=${ANGLE} -DFOOD_DIST=${DIST} -D_FORCE_INLINES

LD = nvcc
LDFLAGS = 

TARGET = ${ANGLE}.deg
CUSOURCES = main.cu IO.cu DataStructures.cu Variables.cu kernel.cu Display.cu
CSOURCES =
CXXSOURCES =
OBJECTS = ${CUSOURCES:.cu=.o} ${CSOURCES:.cu=.o} ${CXXSOURCES:.cpp=.o}

.SUFFIXES: .cu .c .cpp .o

all: ${TARGET}
	-cp ${GLUT_DLL} ./
	-rm -f ${OBJECTS}

${TARGET}: ${OBJECTS}
	${LD} ${LDFLAGS} $^ ${LIBS} -o $@

.cu.o :
	${NVCC} -dc ${NVCCFLAGS} ${INCLUDES} $<

.c.o :
	${CC} -c ${CFLAGS} ${INCLUDES} $<

.cpp.o :
	${CXX} -c ${CXXFLAGS} ${INCLUDES} $<

clean :
	-rm -f  ${OBJECTS} *.deg freeglut.dll
