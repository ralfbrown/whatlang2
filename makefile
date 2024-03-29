## makefile for LangIdent

FRAMEPAC=./framepac
DESTDIR=/usr/bin
DBDIR=/usr/share/langident

SHAREDLIB=

OBJS = 	build/langid.o \
	build/mtrie.o \
	build/prepfile.o \
	build/ptrie.o \
	build/scan_langid.o \
	build/smooth.o \
	build/trie.o \
	build/trigram.o

EXES =	bin/mklangid \
	bin/romanize \
	bin/subsample \
	bin/whatlang

DISTFILES = COPYING README makefile manual.txt *.C *.h \
	mklangid romanize whatlang

LIBRARY=langident.a

#########################################################################
## define the compiler

CC=g++ --std=c++11
CCLINK=$(CC)
CPUTYPE = -D__386__ -D__586__ -D__886__

#########################################################################
## define compilation options

ifndef BUILD_DBG
### compile fully optimized for distribution
BUILD_DBG=0
### compile with debugging info
#BUILD_DBG=1
### compile with debugging info and all optimizations off
#BUILD_DBG=2
endif

# build statically-linked executable (1=yes, 0=no)
#STATIC?=1
STATIC?=0

# enable multi-threading? (1=yes, 0=no)
THREADS?=1
#THREADS?=0

ifndef NODEBUG
#NODEBUG=-DNDEBUG
endif

ifndef GDB
#GDB = -ggdb3
endif

ifeq ($(DO_PROFILE),1)
PROFILE=-pg
#PROFILE=-DPURIFY
else
ifndef PROFILE
PROFILE=
endif
endif

ifeq ($(SANE),1)
SANITIZE=-fsanitize=thread -fPIE -DHELGRIND
LINKSAN=-pie
else ifeq ($(SANE),2)
SANITIZE=-fsanitize=address -fno-omit-frame-pointer -DPURIFY
else ifeq ($(SANE),3)
SANITIZE=-fsanitize=leak -DPURIFY
else ifeq ($(SANE),4)
SANITIZE=-fsanitize=memory -fno-omit-frame-pointer
else ifeq ($(SANE),5)
SANITIZE=-fsanitize=undefined
endif

ifndef CPU
## Uncomment the appropriate CPU type
### 486
#CPU=4
### Pentium
#CPU=5
### PentiumPro or higher
#CPU=6
### AMD Athlon; not supported by GCC 2.x
#CPU=7
### AMD64/x86_64 CPUs in 64-bit mode; not supported by GCC 2.x
###    (AMD K8 [Opteron/Athlon64], newest PentiumIV with EM64t)
#CPU=8
### AMD64 "PhenomII" (K10) or newer
#CPU=10
### Let GCC auto-determine CPU type, but assume at least CPU=8 capabilities
CPU=99
endif

ifndef BITS
#BITS=32
BITS=64
endif

ifeq ($(NOTHREADS),1)
  PTHREAD=-DFrSINGLE_THREADED
else
  PTHREAD=-pthread
#-fopenmp
endif

ifeq ($(STATIC),1)
  LINKTYPE=-static -z muldefs
else
  LINKTYPE=
endif

ifdef MAKE_SHAREDLIB
SHAREDLIB=-fPIC -DSHARED
endif

ifeq ($(NOICONV),1)
ICONV=-DNO_ICONV
else
ICONV=
endif

ifndef RELEASE
RELPATH=LangIdent
ZIPNAME=langident.zip
else
RELPATH=LangIdent-$(RELEASE)
ZIPNAME=langident-$(RELEASE).zip
endif

WARN=-Wall -Wextra -Wno-deprecated -Wshadow -Wcast-align -Wmissing-noreturn -Wmissing-format-attribute
WARN += -Wunused-result -Wno-multichar -Wpacked -Wdisabled-optimization
#WARN += -Wpadded

LINKBITS=-m$(BITS)
CPUDEF=-march=native -D__BITS__=$(BITS)

CC = g++ -std=c++11
CCLINK = $(CC)
CFLAGS = $(WARN)
CFLAGS +=$(CPUDEF)

CFLAGS +=$(PTHREAD)
CFLAGS +=$(PROFILE)
CFLAGS +=$(ICONV)
CFLAGS +=$(NODEBUG)
CFLAGS +=$(LINKBITS) -pipe
CFLAGS +=$(EXTRAINC)
CFLAGS +=$(SANITIZE)
CFLAGS +=$(INCLUDEDIRS)
CFLAGS +=$(SHAREDLIB)
CFLAGS +=$(COMPILE_OPTS)
CFLAGS +=-I$(FRAMEPAC) -I..
CFLAGEXE = -L$(LIBINSTDIR) $(PROFILE) -o $@
LINKFLAGS =$(LINKBITS)
LINKFLAGS +=$(LINKTYPE)
LINKFLAGS +=$(PTHREAD)
LINKFLAGS +=$(SANITIZE)
LINKFLAGS +=$(LINKSAN)
LINKFLAGS +=-L$(FRAMEPAC)

ifeq ($(BUILD_DBG),2)
 CFLAGS += -ggdb3 -O0 -fno-inline -g3 '-DDBDIR="$(DBDIR)"'
 CFLAGSLOOP=
 LINKFLAGS += -ggdb3
else ifeq ($(BUILD_DBG),1)
 CFLAGS +=-Og -fno-inline -ggdb3 '-DDBDIR="$(DBDIR)"'
 CFLAGSLOOP=
 LINKFLAGS += -ggdb3
else ifeq ($(BUILD_DBG),-1)
 CFLAGS +=-O3 -fexpensive-optimizations -ggdb3 '-DDBDIR="$(DBDIR)"'
 CFLAGSLOOP=-fno-align-labels -fno-align-loops
 LINKFLAGS += -ggdb3
else
 CFLAGS +=-O3 -fexpensive-optimizations -fmerge-constants '-DDBDIR="$(DBDIR)"'
 CFLAGSLOOP=-fno-align-labels -fno-align-loops
# CFLAGS += -fweb -ftracer -fgcse-sm -fgcse-las -fno-math-errno
endif

CFLAGSLOOP += $(CFLAGS)

CFLAGEXE=-o $@ -m$(BITS) $(PROFILE)

RM=rm -f
CP=cp -p
NULL=/dev/null

#########################################################################
## define the programs needed to create the target library

LIBRARIAN=ar
LIBFLAGS=rucl
LIBINDEXER = ranlib
LIBIDXFLAGS = $(LIBRARY)

#########################################################################
## standard targets

all:  framepac $(EXES)

clean:
	-$(RM) $(LIBRARY) build/*.o $(EXES)

allclean: clean
	-rmdir bin build >$(NULL) ; true

tags:
	etags --c++ *.h *.C

install:
	$(CP) bin/mklangid bin/romanize bin/whatlang bin/subsample $(DESTDIR)

zip: 	$(EXES)
#	-strip $(EXES)
	-$(RM) $(ZIPNAME)
	mkdir $(RELPATH)
	-$(CP) -i --parents $(DISTFILES) $(RELPATH)
	-( cd $(RELPATH) ; md5sum ${DISTFILES} >MD5SUM )
	-( cd $(RELPATH) ; sha1sum ${DISTFILES} >SHA1SUM )
	zip -mro9 $(ZIPNAME) $(RELPATH)/

lib:	$(LIBRARY)

#########################################################################
## executables

bin/mklangid: build/mklangid.o $(LIBRARY) $(FRAMEPAC)/framepacng.a
	@mkdir -p bin
	$(CCLINK) $(LINKFLAGS) $(CFLAGEXE) -o $@ $^

bin/romanize: build/romanize.o $(LIBRARY) $(FRAMEPAC)/framepacng.a
	@mkdir -p bin
	$(CCLINK) $(LINKFLAGS) $(CFLAGEXE) -o $@ $^

bin/whatlang: build/whatlang.o $(LIBRARY) $(FRAMEPAC)/framepacng.a
	@mkdir -p bin
	$(CCLINK) $(LINKFLAGS) $(CFLAGEXE) -o $@ $^

bin/subsample: build/subsample.o $(FRAMEPAC)/framepacng.a
	@mkdir -p bin
	$(CCLINK) $(LINKFLAGS) $(CFLAGEXE) -o $@ $^

#########################################################################
## object modules

build/langid.o: langid.C langid.h
	$(CC) $(CFLAGSLOOP) -c -o $@ $<

build/mklangid.o: mklangid.C langid.h prepfile.h trie.h mtrie.h ptrie.h

build/whatlang.o: whatlang.C langid.h

build/scan_langid.o: scan_langid.C langid.h

build/mtrie.o: mtrie.C mtrie.h

build/prepfile.o: prepfile.C prepfile.h

build/ptrie.o: ptrie.C ptrie.h mtrie.h

build/romanize.o: romanize.C framepac/framepac/romanize.h

build/smooth.o: smooth.C langid.h

build/subsample.o: subsample.C

build/trie.o: trie.C trie.h

build/trigram.o: trigram.C langid.h trie.h

#########################################################################
## header files -- touching to ensure proper recompilation

langid.h:	mtrie.h ptrie.h
	touch langid.h

#########################################################################
## submodule initialization

framepac:
	@[ -e ../framepac-ng/framepac ] && echo "Linking to local install of FramepaC-ng" && ln -s ../framepac-ng framepac ; true
	@[ -e ../framepac/framepac ] && echo "Linking to local install of FramepaC-ng" && ln -s ../framepac framepac ; true
	@[ ! -e framepac ] && [ -e .git ] && echo "Fetching FramepaC-ng" && git submodule add -f ../framepac-ng.git framepac && git submodule update --init
	@[ -e framepac ] || (echo "Please install FramepaC-ng in subdir 'framepac'" ;exit 1)

#########################################################################
## library

$(LIBRARY): $(OBJS)
	-$(RM) $(LIBRARY)
	$(LIBRARIAN) $(LIBFLAGS) $(LIBRARY) $(OBJS)
	$(LIBINDEXER) $(LIBIDXFLAGS)

$(FRAMEPAC)/framepacng.a:
	( cd $(FRAMEPAC) ; make lib )

#########################################################################
## default compilation rule

.C.o: ; $(CC) $(CFLAGS) $(CPUTYPE) -c $<

build/%.o : %.C
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPUTYPE) -c -o $@ $<
