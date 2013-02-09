TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lconfig++

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
    src/Minimizer/minimizer.cpp

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
    src/Minimizer/minimizer.h

