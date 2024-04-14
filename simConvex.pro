QT -= core
QT -= gui

TARGET = simConvex
TEMPLATE = lib

DEFINES -= UNICODE
DEFINES += QT_COMPIL
DEFINES += SIM_MATH_DOUBLE
CONFIG += shared plugin
INCLUDEPATH += "../include"
INCLUDEPATH += "external/hacd"
INCLUDEPATH += "external/vhacd/inc"
INCLUDEPATH += "external/vhacd/public"
INCLUDEPATH += "external/vhacd/src"
INCLUDEPATH += "external/qHull"

*-msvc* {
    QMAKE_CXXFLAGS += -O2
    QMAKE_CXXFLAGS += -W3
}
*-g++* {
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -fvisibility=hidden
    QMAKE_CXXFLAGS += -Wno-unused-parameter
    QMAKE_CXXFLAGS += -Wno-strict-aliasing
    QMAKE_CXXFLAGS += -Wno-empty-body
    QMAKE_CXXFLAGS += -Wno-write-strings

    QMAKE_CXXFLAGS += -Wno-unused-but-set-variable
    QMAKE_CXXFLAGS += -Wno-unused-local-typedefs
    QMAKE_CXXFLAGS += -Wno-narrowing

    QMAKE_CXXFLAGS += -fpermissive

    QMAKE_CFLAGS += -O3
    QMAKE_CFLAGS += -Wall
    QMAKE_CFLAGS += -Wno-strict-aliasing
    QMAKE_CFLAGS += -Wno-unused-parameter
    QMAKE_CFLAGS += -Wno-unused-but-set-variable
    QMAKE_CFLAGS += -Wno-unused-local-typedefs
}

win32 {
    DEFINES += WIN_SIM
}

macx {
    DEFINES += MAC_SIM
}

unix:!macx {
    DEFINES += LIN_SIM
}

SOURCES += \
    sourceCode/simConvex.cpp \
    ../include/simLib/simLib.cpp \
    ../include/simStack/stackBool.cpp \
    ../include/simStack/stackNull.cpp \
    ../include/simStack/stackNumber.cpp \
    ../include/simStack/stackString.cpp \
    ../include/simStack/stackArray.cpp \
    ../include/simStack/stackMap.cpp \
    ../include/simStack/stackObject.cpp \
    ../include/simLib/scriptFunctionData.cpp \
    ../include/simLib/scriptFunctionDataItem.cpp \
    ../include/simMath/mathFuncs.cpp \
    ../include/simMath/3Vector.cpp \
    ../include/simMath/4Vector.cpp \
    ../include/simMath/7Vector.cpp \
    ../include/simMath/3X3Matrix.cpp \
    ../include/simMath/4X4Matrix.cpp \
    ../include/simMath/mXnMatrix.cpp \
    external/hacd/hacdGraph.cpp \
    external/hacd/hacdHACD.cpp \
    external/hacd/hacdICHull.cpp \
    external/hacd/hacdManifoldMesh.cpp \
    external/hacd/hacdMeshDecimator.cpp \
    external/hacd/hacdMicroAllocator.cpp \
    external/hacd/hacdRaycastMesh.cpp \
    external/vhacd/src/btAlignedAllocator.cpp \
    external/vhacd/src/btConvexHullComputer.cpp \
    external/vhacd/src/FloatMath.cpp \
    external/vhacd/src/VHACD.cpp \
    external/vhacd/src/VHACD-ASYNC.cpp \
    external/vhacd/src/vhacdICHull.cpp \
    external/vhacd/src/vhacdManifoldMesh.cpp \
    external/vhacd/src/vhacdMesh.cpp \
    external/vhacd/src/vhacdRaycastMesh.cpp \
    external/vhacd/src/vhacdVolume.cpp \
    external/qHull/userprintf_rbox.c \
    external/qHull/userprintf.c \
    external/qHull/usermem.c \
    external/qHull/user.c \
    external/qHull/stat.c \
    external/qHull/rboxlib.c \
    external/qHull/random.c \
    external/qHull/qset.c \
    external/qHull/poly2.c \
    external/qHull/poly.c \
    external/qHull/merge.c \
    external/qHull/mem.c \
    external/qHull/libqhull.c \
    external/qHull/io.c \
    external/qHull/global.c \
    external/qHull/geom2.c \
    external/qHull/geom.c \

HEADERS +=\
    sourceCode/simConvex.h \
    ../include/simLib/simLib.h \
    ../include/simStack/stackBool.h \
    ../include/simStack/stackNull.h \
    ../include/simStack/stackNumber.h \
    ../include/simStack/stackString.h \
    ../include/simStack/stackArray.h \
    ../include/simStack/stackMap.h \
    ../include/simStack/stackObject.h \
    ../include/simLib/scriptFunctionData.h \
    ../include/simLib/scriptFunctionDataItem.h \
    ../include/simMath/mathFuncs.h \
    ../include/simMath/mathDefines.h \
    ../include/simMath/3Vector.h \
    ../include/simMath/4Vector.h \
    ../include/simMath/7Vector.h \
    ../include/simMath/3X3Matrix.h \
    ../include/simMath/4X4Matrix.h \
    ../include/simMath/mXnMatrix.h \
    external/hacd/hacdCircularList.h \
    external/hacd/hacdGraph.h \
    external/hacd/hacdHACD.h \
    external/hacd/hacdICHull.h \
    external/hacd/hacdManifoldMesh.h \
    external/hacd/hacdMeshDecimator.h \
    external/hacd/hacdMicroAllocator.h \
    external/hacd/hacdRaycastMesh.h \
    external/hacd/hacdSArray.h \
    external/hacd/hacdVector.h \
    external/hacd/hacdVersion.h \
    external/vhacd/public/VHACD.h \
    external/vhacd/inc/btAlignedAllocator.h \
    external/vhacd/inc/btAlignedObjectArray.h \
    external/vhacd/inc/btConvexHullComputer.h \
    external/vhacd/inc/btMinMax.h \
    external/vhacd/inc/btScalar.h \
    external/vhacd/inc/btVector3.h \
    external/vhacd/inc/FloatMath.h \
    external/vhacd/inc/vhacdCircularList.h \
    external/vhacd/inc/vhacdlCHull.h \
    external/vhacd/inc/vhacdManifoldMesh.h \
    external/vhacd/inc/vhacdMesh.h \
    external/vhacd/inc/vhacdMutex.h \
    external/vhacd/inc/vhacdRaycastMesh.h \
    external/vhacd/inc/vhacdSArray.h \
    external/vhacd/inc/vhacdTimer.h \
    external/vhacd/inc/vhacdVector.h \
    external/vhacd/inc/vhacdVHACD.h \
    external/vhacd/inc/vhacdVolume.h \
    external/qHull/user.h \
    external/qHull/stat.h \
    external/qHull/random.h \
    external/qHull/qset.h \
    external/qHull/qhull_a.h \
    external/qHull/poly.h \
    external/qHull/merge.h \
    external/qHull/mem.h \
    external/qHull/libqhull.h \
    external/qHull/io.h \
    external/qHull/geom.h \

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
