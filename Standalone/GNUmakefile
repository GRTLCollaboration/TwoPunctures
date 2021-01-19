# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)
# Also need to set the GRCHOMBO_SOURCE environment variable
ifndef GRCHOMBO_SOURCE
    $(error Please define GRCHOMBO_SOURCE.)
endif

# Note this code requires linking with the GNU Scientific Library
# In addition to the line below, the include paths to the relevant header files
# might also need to be passed to the preprocessor/compiler
XTRALIBFLAGS = -lgsl

TWOPUNCTURES_SOURCE = ../Source

ebase := Main_TwoPunctures

LibNames :=

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(TWOPUNCTURES_SOURCE)


include $(CHOMBO_HOME)/mk/Make.test
