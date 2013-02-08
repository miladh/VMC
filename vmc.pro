TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lconfig++

SOURCES += src/main.cpp \
    src/vmcsolver.cpp \
    src/includes/lib.cpp \
    src/Wavefunction/wavefunction.cpp

HEADERS += \
    src/vmcsolver.h \
    src/includes/lib.h \
    src/Wavefunction/wavefunction.h

