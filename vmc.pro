TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += src/main.cpp \
    src/vmcsolver.cpp \
    src/includes/lib.cpp

HEADERS += \
    src/vmcsolver.h \
    src/includes/lib.h

