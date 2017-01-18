TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ec.cpp \
    point.cpp

HEADERS += \
    ec.hpp \
    test.hpp \
    point.hpp

LIBS += -lgcrypt -lgmp

QMAKE_CXXFLAGS += -std=c++1y
