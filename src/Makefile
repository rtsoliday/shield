.DEFAULT_GOAL := all
SDDS_REPO = $(firstword $(wildcard ../../SDDS ../../../../../epics/extensions/src/SDDS))
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the shield repository)
endif

ENABLE_OPENMP ?= 1
BUNDLE_OPENMP ?= 1
include ../Makefile.rules
PROD = shield

shield_SRC = shield.c

ifneq ($(wildcard $(DD)../.svn $(DD)../../../../.svn),)
  SVN_VERSION=$(shell svnversion .)
  ifeq ($(SVN_VERSION),)
    SVN_VERSION=unknown
  endif
  ifeq ($(SVN_VERSION),Unversioned directory)
    SVN_VERSION=unknown
  endif
else
  SVN_VERSION=GIT_$(shell git log --oneline | wc -l | sed 's/^[[:space:]]*//')
endif

CFLAGS += -I$(SDDS_REPO)/include -I$(OBJ_DIR) -DSVN_VERSION=\"$(SVN_VERSION)\"

ifeq ($(OS), Linux)
  CFLAGS += -Wno-unused-variable
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = -lmdbcommon -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Darwin)
  CFLAGS += -Wno-unused-variable
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = -lmdbcommon -lmatlib -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Windows)
  CFLAGS += -I$(SDDS_REPO)/lzma -wd4244
  PROD_LIBS += gsl.lib gslcblas.lib
  PROD_LIBS_SDDS = mdbcommon.lib matlib.lib SDDS1.lib namelist.lib rpnlib.lib mdbmth.lib mdblib.lib lzma.lib z.lib
  LIB_LINK_DIRS += -LIBPATH:$(SDDS_REPO)/lib/$(OS)-$(ARCH)
endif

CUSTOM_LINK_PRODUCTS += shield
BUILD_GENERATED += $(patsubst %.nl,$(OBJ_DIR)/%.h,$(wildcard *.nl))
include ../Makefile.build
$(BUILD_GENERATED): | $(OBJ_DIR)
$(OBJ_DIR)/shield$(EXEEXT): $(shield_OBJS) $(PROD_DEPS)
	$(LINKEXE) $(OUTPUTEXE) $(shield_OBJS) $(LDFLAGS) $(LIB_LINK_DIRS) $(PROD_LIBS) $(PROD_LIBS_SDDS) $(PROD_SYS_LIBS)
	cp -f $@ $(BIN_DIR)/
	@if [ -n "$(EPICS_BIN_DIR)" ]; then echo cp -f $@ $(EPICS_BIN_DIR)/; fi
	@if [ -n "$(EPICS_BIN_DIR)" ]; then cp -f $@ $(EPICS_BIN_DIR)/; fi

$(OBJ_DIR)/shield.$(OBJEXT): shield.c $(OBJ_DIR)/shield.h

$(OBJ_DIR)/shield.h: shield.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp shield.nl $(OBJ_DIR)/shield.h

