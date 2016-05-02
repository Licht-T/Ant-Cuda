
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
NVCCFLAGS = -O3 -arch sm_35 -DMACRO_FOOD_ANGLE=${ANGLE} -DMACRO_FOOD_DIST=${DIST} -D_FORCE_INLINES

CC = cc
CFLAGS = -O3 -DMACRO_FOOD_ANGLE=${ANGLE} -DMACRO_FOOD_DIST=${DIST} -D_FORCE_INLINES

CXX = g++
CXXFLAGS = -O3 -DMACRO_FOOD_ANGLE=${ANGLE} -DMACRO_FOOD_DIST=${DIST} -D_FORCE_INLINES

LD = nvcc
LDFLAGS = -arch sm_35

TARGET = ${DIST}dist_${ANGLE}deg.exe
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
	-rm -f  ${OBJECTS} *.exe freeglut.dll
