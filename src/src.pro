TEMPLATE = app
CONFIG += console
CONFIG -= qt

include(../defaults.pri)

TARGET = vmc

SOURCES += main.cpp \
    includes/lib.cpp \
    Wavefunction/wavefunction.cpp \
    Potential/potential.cpp \
    Potential/coulombPotential.cpp \
    Kinetic/kinetic.cpp \
    Solver/solver.cpp \
    Solver/mcbf.cpp \
    Hamiltonian/hamiltonian.cpp \
    VMCApp/vmcapp.cpp \
    Minimizer/minimizer.cpp \
    Solver/mcis.cpp \
    Orbitals/orbitals.cpp \
    Orbitals/hydrogenic.cpp \
    slater/slater.cpp \
    Jastrow/padejastrow.cpp \
    Jastrow/jastrow.cpp \
    Jastrow/nojastrow.cpp \
    Wavefunction/hlikewavefunction.cpp \
    OnebodyDensity/onebodydensity.cpp \
    Blocking/blocking.cpp \
    electronInteraction/electroninteraction.cpp \
    electronInteraction/coulombinteraction.cpp \
    electronInteraction/nointeraction.cpp \
    Minimizer/bfminimizer.cpp \
    Minimizer/steepestdescent.cpp \
    Observables/observables.cpp \
    ConfigurationParser/configurationparser.cpp

HEADERS += \
    includes/lib.h \
    Wavefunction/wavefunction.h \
    Potential/potential.h \
    Kinetic/kinetic.h \
    Solver/solver.h \
    Solver/mcbf.h \
    Hamiltonian/hamiltonian.h \
    VMCApp/vmcapp.h \
    Minimizer/minimizer.h \
    Solver/mcis.h \
    includes/Defines.h \
    Orbitals/orbitals.h \
    Orbitals/hydrogenic.h \
    slater/slater.h \
    Jastrow/padejastrow.h \
    Jastrow/jastrow.h \
    Jastrow/nojastrow.h \
    Wavefunction/hlikewavefunction.h \
    OnebodyDensity/onebodydensity.h \
    Blocking/blocking.h \
    electronInteraction/electroninteraction.h \
    electronInteraction/coulombinteraction.h \
    electronInteraction/nointeraction.h \
    Potential/coulombPotential.h \
    Minimizer/bfminimizer.h \
    Minimizer/steepestdescent.h \
    Observables/observables.h \
    ConfigurationParser/configurationparser.h



OTHER_FILES += \
    config.cfg
