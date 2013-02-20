TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lconfig++

TARGET = vmc

QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x

SOURCES += src/main.cpp \
    src/includes/lib.cpp \
    src/Wavefunction/wavefunction.cpp \
    src/Wavefunction/basicwavefunction.cpp \
    src/Wavefunction/jastrowwavefunction.cpp \
    src/Potential/potential.cpp \
    src/Potential/coulomb_potential.cpp \
    src/Kinetic/kinetic.cpp \
    src/Kinetic/numericalkinetic.cpp \
    src/Kinetic/closedformkinetic.cpp \
    src/Solver/solver.cpp \
    src/Solver/mcbf.cpp \
    src/Hamiltonian/hamiltonian.cpp \
    src/VMCApp/vmcapp.cpp \
    src/Minimizer/minimizer.cpp \
    src/Wavefunction/hydrogenicwavefunction.cpp \
    src/Solver/mcis.cpp

HEADERS += \
    src/includes/lib.h \
    src/Wavefunction/wavefunction.h \
    src/Wavefunction/basicwavefunction.h \
    src/Wavefunction/jastrowwavefunction.h \
    src/Potential/potential.h \
    src/Potential/coulomb_potential.h \
    src/Kinetic/kinetic.h \
    src/Kinetic/numericalkinetic.h \
    src/Kinetic/closedformkinetic.h \
    src/Solver/solver.h \
    src/Solver/mcbf.h \
    src/Hamiltonian/hamiltonian.h \
    src/VMCApp/vmcapp.h \
    src/Minimizer/minimizer.h \
    src/Wavefunction/hydrogenicwavefunction.h \
    src/Solver/mcis.h


release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}




# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
