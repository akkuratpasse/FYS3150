TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += /usr/local/include/
LIBS += -L/usr/local/lib -larmadillo

SOURCES += main.cpp \
    ../libraries/lib.cpp
HEADERS += \
    ../libraries/lib.h
