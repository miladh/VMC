TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lconfig++

SOURCES += src/main.cpp \
    src/vmcsolver.cpp \
    src/includes/lib.cpp \
    src/Wavefunction/wavefunction.cpp \
    src/Wavefunction/basicwavefunction.cpp \
    src/Wavefunction/jastrowwavefunction.cpp \
    src/Potential/potential.cpp \
    src/Potential/coulomb_potential.cpp \
    src/Kinetic/kinetic.cpp \
    src/Kinetic/numericalkinetic.cpp \
    src/Kinetic/closedformkinetic.cpp

HEADERS += \
    src/vmcsolver.h \
    src/includes/lib.h \
    src/Wavefunction/wavefunction.h \
    src/Wavefunction/basicwavefunction.h \
    src/Wavefunction/jastrowwavefunction.h \
    src/Potential/potential.h \
    src/Potential/coulomb_potential.h \
    src/Kinetic/kinetic.h \
    src/Kinetic/numericalkinetic.h \
    src/Kinetic/closedformkinetic.h

