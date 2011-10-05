# common SWIG language build makefile
# Author: Karl Wette, 2011

# first rule
.PHONY : all
all :

# link this makefile to Makefile.in for automatic remaking
$(srcdir)/Makefile.in : $(top_srcdir)/swig/swig-common.mk

# the library SWIG wrappings are being generated for
swig_lib = $(PACKAGE_NAME)

# name of the output SWIG wrapping module
swig_wrapname = swig$(swig_lib)

# name of the SWIG interface file
swig_ifacename = swig-$(swig_lib)
swig_ifacefile = $(swig_ifacename).i
swig_iface = $(top_builddir)/swig/$(swig_ifacefile)

# SWIG interface file dependencies
swig_iface_deps = .$(swig_ifacename).deps

# name of the SWIG wrapping code
swig_wrapfile = $(swig_wrapname).cpp

# name of the test script for the SWIG wrapping module
swig_wrapcheck = $(srcdir)/check-$(swig_wrapname)

# cleanup files
swig_cleanfiles = $(swig_wrapfile) \
                  $(swig_iface) \
                  $(swig_iface_deps)

# name of scripting language being targeted
# if this is defined, this is a language Makefile
# which has been configured for use
ifdef swig_language

# name of an output generated by compiling the
# SWIG wrapping code, used for make dependencies
ifndef swig_wrapobj
$(error 'swig_wrapobj' must be defined)
endif

# build rules
all : $(swig_wrapobj)
$(swig_wrapobj) : $(swig_wrapfile)

# install rule
install-exec-local : $(swig_wrapobj)

# SWIG interface file preamble
swig_iface_pre = $(top_srcdir)/swig/$(swig_ifacename).swg

# SWIG interface header include files
swig_iface_headers = $(wildcard $(top_builddir)/swig/swig-iface-header-*.swg)

# build SWIG interface file
$(swig_iface) : $(swig_iface_pre) $(swig_iface_headers)
	@cat $+ > $@

# generate dummy dependencies file
$(swig_iface_deps) : $(swig_iface)
	@echo "#$(swig_iface_deps)" > $@
include $(swig_iface_deps)

# script which checks headers for correct SWIGLAL constructs
swig_check_headers = $(top_srcdir)/swig/swig-check-headers.pl

# path where SWIG should look for SWIG interface / LAL header files
swig_inclpath = $(top_srcdir)/swig \
                $(top_builddir)/include \
                $(wildcard $(sort $(SWIG_INCLPATH)))

# generate SWIG wrapping code
$(swig_wrapfile) : $(swig_iface_deps) $(swig_check_headers)
	@cd $(top_builddir)/swig && $(MAKE) iface-links
	CPP='$(CPP)' $(PERL) $(swig_check_headers) \
	--include $(top_builddir)/include --interface $(swig_iface)
	$(SWIG) $(addprefix -D,$(SWIG_SWIG_DEFINES)) \
        -module $(swig_wrapname) -$(swig_language) -c++ $(swig_language_opts) \
	-MD -MF $(swig_iface_deps) -MT $(swig_wrapfile) \
	$(addprefix -I,$(swig_inclpath)) \
	-o $(swig_wrapfile) $(swig_iface)

# symbols to define when compiling SWIG wrapping code
swig_defines = $(SWIG_CXX_DEFINES) \
               $(subst -D,,$(DEFS))

# path where SWIG should look for LAL libraries while linking
swig_libpath = $(wildcard $(sort $(SWIG_PREINST_LIBPATH) $(SWIG_LIBPATH)))

# directory where LAL libraries will eventually be installed
swig_libdir = $(libdir)

# libraries that SWIG wrapping code should link against
swig_libs = $(sort lal lalsupport $(swig_lib))

# export variables to SWIG wrapping library build script
export build_vendor
export PACKAGE_NAME
export PACKAGE_VERSION
export swig_wrapname
export swig_defines
export swig_inclpath
export swig_libpath
export swig_libdir
export swig_libs
export swig_wrapfile

# build colon-separated path by repeated concatenation
swig_makepath = $(if $(word 2,$1),$(word 1,$1):$(call swig_makepath,$(wordlist 2,$(words $1),$1)),$(word 1,$1))

# set library load path when running check scripts prior to installation
ifeq "$(build_vendor)" "apple"
swig_ldlibpathname = DYLD_LIBRARY_PATH
else
swig_ldlibpathname = LD_LIBRARY_PATH
endif
swig_ldlibpath = $(call swig_makepath,$(SWIG_PREINST_LIBPATH))

endif # ifdef swig_language
