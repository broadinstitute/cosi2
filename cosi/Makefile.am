#
# Makefile: cosi makefile
#
# Describes how to build the cosi library, libcosi, and the main cosi program, 'coalescent'.
#
# Run "make doc" to build HTML documentation.
#
# Most of cosi code is compiled into a library, libcosi, with which the main
# executable is then linked.
#

#SUBDIRS=libsigc++-2.3.1

ACLOCAL_AMFLAGS = -I m4
BOOST_CPPFLAGS=-isystem $(srcdir)/boost
AM_CPPFLAGS = -I$(srcdir)/.. $(BOOST_CPPFLAGS) -DNDEBUG 
AM_CXXFLAGS = $(CXX11_CXXFLAGS) -O3
COSI_CPPFLAGS = $(AM_CPPFLAGS) -DCOSI_TRACK_LAST_COAL #-I$(builddir)/libsigc++-2.3.1 -I$(srcdir)/libsigc++-2.3.1 
COSI_CXXFLAGS= $(AM_CXXFLAGS)

EXTRA_DIST=
lib_LTLIBRARIES=
bin_PROGRAMS=

AM_LDFLAGS = -L$(srcdir)/../cosi_rand
AM_DEFAULT_SOURCE_EXT = .cc

####################################################################################################
#
# Section: Building Boost libraries
#
# Build Boost libraries that we use and that are not header-only.
# Building them as part of our build reduces dependence on the build environment and ensures
# that they are built with the same flags as the code that uses them, preventing link-time issues.
#
####################################################################################################

BOOST_PROGRAM_OPTIONS_LIBS = $(builddir)/libprogopt.la
BOOST_REGEX_LIBS = $(builddir)/libcosiregex.la
BOOST_FILESYSTEM_LIBS = $(builddir)/libcosifilesystem.la
BOOST_SYSTEM_LIBS = $(builddir)/libcosisystem.la
BOOST_CHRONO_LIBS = $(builddir)/libcosichrono.la
BOOST_TIMER_LIBS = $(builddir)/libcositimer.la
BOOST_SERIALIZATION_LIBS = $(builddir)/libcosiserialization.la

lib_LTLIBRARIES += libprogopt.la libcosiregex.la libcosifilesystem.la libcosisystem.la libcosichrono.la libcositimer.la libcosiserialization.la

PO=$(srcdir)/boost/libs/program_options/src
libprogopt_la_SOURCES = $(PO)/cmdline.cpp $(PO)/config_file.cpp $(PO)/convert.cpp $(PO)/options_description.cpp $(PO)/parsers.cpp $(PO)/positional_options.cpp $(PO)/split.cpp $(PO)/utf8_codecvt_facet.cpp $(PO)/value_semantic.cpp $(PO)/variables_map.cpp $(PO)/winmain.cpp

RE=$(srcdir)/boost/libs/regex/src
libcosiregex_la_SOURCES = $(RE)/cpp_regex_traits.cpp $(RE)/cregex.cpp $(RE)/c_regex_traits.cpp $(RE)/fileiter.cpp $(RE)/icu.cpp \
	$(RE)/instances.cpp $(RE)/internals.hpp $(RE)/posix_api.cpp $(RE)/regex.cpp $(RE)/regex_debug.cpp $(RE)/regex_raw_buffer.cpp \
	$(RE)/regex_traits_defaults.cpp $(RE)/static_mutex.cpp $(RE)/usinstances.cpp $(RE)/w32_regex_traits.cpp $(RE)/wc_regex_traits.cpp \
	$(RE)/wide_posix_api.cpp $(RE)/winstances.cpp
libcosiregex_la_CXXFLAGS = $(AM_CXXFLAGS) -Wno-strict-overflow

FS=$(srcdir)/boost/libs/filesystem/src
SY=$(srcdir)/boost/libs/system/src
libcosifilesystem_la_SOURCES = $(FS)/codecvt_error_category.cpp $(FS)/operations.cpp $(FS)/path.cpp $(FS)/path_traits.cpp $(FS)/portability.cpp $(FS)/unique_path.cpp $(FS)/utf8_codecvt_facet.cpp $(FS)/windows_file_codecvt.cpp $(SY)/error_code.cpp
libcosifilesystem_la_CPPFLAGS=$(AM_CPPFLAGS)

libcosisystem_la_SOURCES = $(SY)/error_code.cpp $(SY)/local_free_on_destruction.hpp
libcosisystem_la_CPPFLAGS=$(AM_CPPFLAGS)

CR=$(srcdir)/boost/libs/chrono/src
libcosichrono_la_SOURCES = $(CR)/chrono.cpp $(CR)/process_cpu_clocks.cpp $(CR)/thread_clock.cpp
libcosichrono_la_CPPFLAGS = $(AM_CPPFLAGS)
#libcosichrono_la_LDFLAGS = -lrt

TI=$(srcdir)/boost/libs/timer/src
libcositimer_la_SOURCES = $(TI)/auto_timers_construction.cpp $(TI)/cpu_timer.cpp
libcositimer_la_CPPFLAGS = $(AM_CPPFLAGS)
libcositimer_la_LIBADD = $(BOOST_CHRONO_LIBS)

SR=$(srcdir)/boost/libs/serialization/src
libcosiserialization_la_SOURCES = $(SR)/archive_exception.cpp $(SR)/basic_archive.cpp $(SR)/basic_iarchive.cpp $(SR)/basic_iserializer.cpp $(SR)/basic_oarchive.cpp $(SR)/basic_oserializer.cpp \
 $(SR)/basic_pointer_iserializer.cpp $(SR)/basic_pointer_oserializer.cpp $(SR)/basic_serializer_map.cpp $(SR)/basic_text_iprimitive.cpp $(SR)/basic_text_oprimitive.cpp \
 $(SR)/basic_text_wiprimitive.cpp $(SR)/basic_text_woprimitive.cpp $(SR)/basic_xml_archive.cpp $(SR)/binary_iarchive.cpp $(SR)/binary_oarchive.cpp \
 $(SR)/binary_wiarchive.cpp $(SR)/binary_woarchive.cpp $(SR)/codecvt_null.cpp $(SR)/extended_type_info.cpp $(SR)/extended_type_info_no_rtti.cpp $(SR)/extended_type_info_typeid.cpp \
 $(SR)/polymorphic_iarchive.cpp $(SR)/polymorphic_oarchive.cpp $(SR)/shared_ptr_helper.cpp $(SR)/stl_port.cpp $(SR)/text_iarchive.cpp $(SR)/text_oarchive.cpp $(SR)/text_wiarchive.cpp \
 $(SR)/text_woarchive.cpp $(SR)/utf8_codecvt_facet.cpp $(SR)/void_cast.cpp $(SR)/xml_archive_exception.cpp $(SR)/xml_grammar.cpp $(SR)/xml_iarchive.cpp $(SR)/xml_oarchive.cpp \
 $(SR)/xml_wgrammar.cpp $(SR)/xml_wiarchive.cpp $(SR)/xml_woarchive.cpp
libcosiserialization_la_CPPFLAGS = $(AM_CPPFLAGS)

EXTRA_DIST += $(srcdir)/boost

################################################
# End section: Building Boost libraries
################################################

################################################
#
# Section: Building cosi
#
################################################

lib_LTLIBRARIES += libcosiutil.la libcosi.la 
bin_PROGRAMS += coalescent

if DEVOPTS
COSI_CXXFLAGS += --param max-inline-insns-single=1000 \
	-Werror -Wall -Wextra -Wno-unused-function -Wuninitialized \
	-Wredundant-decls -Wpointer-arith -Wcast-qual -Wcast-align -Wformat \
	-Woverloaded-virtual -Wstrict-null-sentinel -Wlogical-op \
	-Wmissing-declarations -Wdouble-promotion -Winit-self

endif  # if DEVOPTS

#if COALAPX
# Turn on support for approximating the coalescent by restricting coalescence to nodes whose
# genetic material overlaps or nearly-overlaps.
COSI_CPPFLAGS += -DCOSI_SUPPORT_COALAPX
#endif

#if CONDSNP
# Turn on support for conditioning the generated simulations on those that have a SNP at a
# specified location with specified sample frequencies.
#COSI_CPPFLAGS += -DCOSI_CONDSNP
#endif

#if EHH
# Turn on support for computing EHH and related statistics.
#COSI_CPPFLAGS += -DCOSI_EHH -DCOSI_LEAFSET_BITSET
#endif

COSIUTILHEADERS = utils.h arrival.h generalmath.h cosirand.h gauss_legendre.h valgrind.h evtgrp.h typedval.h \
	order_statistics.hpp random_access_skip_list.h skip_list_detail.h utildefs.h drd.h memcheck.h mempool.h 
libcosiutil_la_SOURCES = utils.cc arrival.cc cosirand.cc gauss_legendre.c mempool.cc $(COSIUTILHEADERS)

COSISOURCES = $(srcdir)/coalesce.cc $(srcdir)/coalescent.cc $(srcdir)/demography.cc $(srcdir)/file.cc $(srcdir)/geneconversion.cc \
	$(srcdir)/historical.cc \
	$(srcdir)/migrate.cc $(srcdir)/mutate.cc $(srcdir)/mutlist.cc $(srcdir)/node.cc $(srcdir)/output.cc $(srcdir)/pop.cc $(srcdir)/recomb.cc \
	$(srcdir)/seglist.cc $(srcdir)/segsumm.cc $(srcdir)/simulator.cc \
	$(srcdir)/sweep.cc $(srcdir)/sweep1.cc $(srcdir)/sweep2.cc $(srcdir)/sweep3.cc $(srcdir)/genmap.cc $(srcdir)/stats.cc $(srcdir)/hooks.cc \
	$(srcdir)/traj.cc $(srcdir)/cositop.cc $(srcdir)/mutcontext.cc \
	$(srcdir)/condsnp.cc $(srcdir)/hullmgr.cc $(srcdir)/utils.cc $(srcdir)/arrival.cc $(srcdir)/cosirand.cc $(srcdir)/gauss_legendre.c \
	$(srcdir)/mempool.cc $(srcdir)/leafset-tree.cc $(srcdir)/leafset-counts.cc

libcosi_la_SOURCES = coalesce.cc coalescent.cc demography.cc file.cc geneconversion.cc historical.cc \
	migrate.cc mutate.cc mutlist.cc node.cc output.cc pop.cc recomb.cc seglist.cc segsumm.cc simulator.cc \
	sweep.cc sweep1.cc sweep2.cc sweep3.cc genmap.cc stats.cc hooks.cc traj.cc cositop.cc \
	callgrind.h defs.h coalesce.h demography.h ehh.h ehh.cc file.h geneconversion.h historical.h leafset.h \
	leafset-bitset.h leafset-tree.h leafset-sizeonly.h leafset-sizeonly.cc leafset-counts.h mutcontext.cc \
	migrate.h mutate.h node.h nodefwd.h nodelist.h output.h pop.h popfwd.h recomb.h seglist.h seglistfwd.h \
	seglistfwd2.h segsumm.h simulator.h cositop.h \
	sweep.h sweep1.h sweep2.h sweep3.h recomb.h genmap.h stats.h traj.h decls.h hooks.h \
	coalrate.h condsnp.h condsnp.cc \
	mutcontext.h hullmgr.h hullmgr.cc \
	cosicfg.h $(COSIUTILHEADERS)

libcosi_la_CPPFLAGS = $(COSI_CPPFLAGS)
libcosi_la_CXXFLAGS = $(COSI_CXXFLAGS)
libcosi_la_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
libcosi_la_LIBADD = $(builddir)/../cosi_rand/libcosi_rand.la $(builddir)/libcosiutil.la $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LIBS) $(BOOST_TIMER_LIBS) #$(builddir)/libsigc++-2.3.1/sigc++/libsigc-2.0.la

coalescent_SOURCES = cosimain.cc
coalescent_CPPFLAGS = $(COSI_CPPFLAGS)
coalescent_CXXFLAGS = $(COSI_CXXFLAGS)
coalescent_LDADD = $(builddir)/libcosi.la $(builddir)/../cosi_rand/libcosi_rand.la $(BOOST_PROGRAM_OPTIONS_LIBS)	\
	$(BOOST_FILESYSTEM_LIBS) $(BOOST_TIMER_LIBS) $(BOOST_SERIALIZATION_LIBS) #$(builddir)/libsigc++-2.3.1/sigc++/libsigc-2.0.la
coalescent_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)

subpkgincludedir = $(includedir)/cosi
subpkginclude_HEADERS = coalescent.h mutlist.h defs.h utildefs.h cositop.h

#if EHH
#libcosi_la_SOURCES += leafset-bitset.cc
#else
libcosi_la_SOURCES += leafset-tree.cc leafset-counts.cc
#endif

#COSI_CXX_FLAGS_DEBUG += -O0 -UNDEBUG -D_GLIBCXX_DEBUG -fno-inline-functions -ggdb3

COSISTATICSOURCES=$(COSISOURCES) $(srcdir)/cosimain.cc $(libprogopt_la_SOURCES) $(libcosifilesystem_la_SOURCES) $(libcosiserialization_la_SOURCES) $(libcositimer_la_SOURCES) $(libcosichrono_la_SOURCES) $(srcdir)/../cosi_rand/mtwist.c $(srcdir)/../cosi_rand/poisson.c $(srcdir)/../cosi_rand/ranbinom.c $(srcdir)/../cosi_rand/gammln.c $(srcdir)/../cosi_rand/rng.c $(srcdir)/../cosi_rand/expdev.c $(srcdir)/../cosi_rand/gamma.c

coalescent_static: $(COSISTATICSOURCES)
	g++ -v -static -o $@ -I$(srcdir)/.. -I$(srcdir)/boost $(COSI_CPPFLAGS) $(CPPFLAGS) $(COSI_CXXFLAGS) $(CXXFLAGS)  -Wno-error -pipe $(COSISTATICSOURCES) -lrt -lpthread

# End section: Building cosi


###########################################################################################################################

#
# Section: Tests
#

check_PROGRAMS =
check_SCRIPTS =

include cositests.am

check_PROGRAMS += seglist_test funcs_test test_funcs test_arrproc
seglist_test_SOURCES = seglist_test.cc seglist.cc seglist.h leafset-tree.h leafset-tree.cc leafset-counts.cc leafset-counts.h leafset.h segsumm.h segsumm.cc \
	utils.h utils.cc defs.h utildefs.h typedval.h mempool.h mempool.cc
seglist_test_CPPFLAGS = $(COSI_CPPFLAGS) -UNDEBUG -D_GLIBCXX_DEBUG
seglist_test_CXXFLAGS = $(COSI_CXXFLAGS)
seglist_test_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
seglist_test_LDADD = $(builddir)/../cosi_rand/libcosi_rand_noinst.la $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_FILESYSTEM_LIBS)

funcs_test_SOURCES = funcs_test.cc generalmath.h utils.h utils.cc
funcs_test_CPPFLAGS = $(COSI_CPPFLAGS) -UNDEBUG -D_GLIBCXX_DEBUG
funcs_test_LDADD = $(BOOST_SYSTEM_LIBS)

test_funcs_SOURCES = test_funcs.cc generalmath.h utils.h utils.cc
test_funcs_CPPFLAGS = $(COSI_CPPFLAGS) -UNDEBUG -D_GLIBCXX_DEBUG -DCOSI_DEV_PRINT
test_funcs_LDADD = $(BOOST_SYSTEM_LIBS)

test_arrproc_SOURCES = test_arrproc.cc generalmath.h utils.h utils.cc arrproc.h typeutil.h
test_arrproc_CPPFLAGS = $(COSI_CPPFLAGS) -DCOSI_DEV_PRINT -UNDEBUG
test_arrproc_LDADD = $(BOOST_SYSTEM_LIBS) #$(builddir)/libsigc++-2.3.1/sigc++/libsigc-2.0.la

EXTRA_DIST += $(TESTS) $(COSI_TESTDATA) $(srcdir)/make_cosi_tests.py $(srcdir)/update_cosi_tests.py $(srcdir)/runtest.py \
	$(srcdir)/make_cosi_tests2.py $(srcdir)/update_tests.py $(srcdir)/cmpdist.py $(srcdir)/timecmd.sh

TESTS += seglist_test test_funcs

if HAVE_PYTHON
if HAVE_CMPDIST
check_SCRIPTS += cmpdist.py runtest.py timecmd.sh
endif
endif

#bin_PROGRAMS += test_random_access
#test_random_access_SOURCES = test_random_access.cpp random_access_skip_list.h skip_list_detail.h test_types.h catch.hpp 
#test_random_access_CXXFLAGS = -Wno-unknown-pragmas

check_PROGRAMS += ost_test
ost_test_SOURCES = ost_test.cpp order_statistics.hpp utils.h utils.cc
ost_test_LDADD = $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(builddir)/../cosi_rand/libcosi_rand.la
ost_test_CPPFLAGS = $(AM_CPPFLAGS) -UNDEBUG

if HAVE_GTEST
check_PROGRAMS += ost_test_main
ost_test_main_SOURCES = ost_test_main.cpp order_statistics.hpp utils.h utils.cc
ost_test_main_CPPFLAGS = $(AM_CPPFLAGS) -UNDEBUG -I$(GTEST_DIR)/include
ost_test_main_LDADD = $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(builddir)/../cosi_rand/libcosi_rand.la
ost_test_main_LDFLAGS = $(AM_LDFLAGS) -L$(GTEST_DIR)/lib -lgtest
endif  # if HAVE_GTEST

#
# Section: Building documentation
#

EXTRA_DIST += NaturalDocs ndproj

doc:
	$(mkdir_p) $(builddir)/doc/html
	$(srcdir)/NaturalDocs/NaturalDocs -p $(srcdir)/ndproj -i $(srcdir) -xi $(srcdir)/boost -xi $(srcdir)/NaturalDocs -o FramedHTML $(builddir)/doc/html

.PHONY: doc


#
# Section: Development/experimental options
#

if HAVEROOT
bin_PROGRAMS += cosiroot ldhist
endif

bin_PROGRAMS += sample_stats_extra


sample_stats_extra_SOURCES = sample_stats_extra.cc tajd.cc
sample_stats_extra_CPPFLAGS = $(AM_CPPFLAGS) #-I @LIBSEQUENCE_INC_DIR@
sample_stats_extra_CXXFLAGS = $(COSI_CXXFLAGS)
sample_stats_extra_LDFLAGS = $(BOOST_REGEX_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 
sample_stats_extra_LDADD = $(BOOST_REGEX_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) 

if HAVEROOT
cosiroot_SOURCES = cosiroot.cc
cosiroot_CPPFLAGS = $(AM_CPPFLAGS) $(ROOTCFLAGS)
cosiroot_LDADD = $(coalescent_LDADD) -L$(ROOTLIBDIR) $(ROOTLIBS)

ldhist_SOURCES = ldhist.cc
ldhist_CPPFLAGS = $(AM_CPPFLAGS) $(ROOTCFLAGS)
ldhist_LDADD = $(coalescent_LDADD) -L$(ROOTLIBDIR) $(ROOTLIBS)

endif

COSIFLAGS = $(builddir)/coalescent.flags.txt

$(COSIFLAGS): $(builddir)/coalescent
	echo "coalescent " $(coalescent_CPPFLAGS) $(CPPFLAGS) $(coalescent_CXXFLAGS) $(CXXFLAGS) $(coalescent_LDFLAGS) $(LDFLAGS) > $@
	echo "libcosi " $(libcosi_la_CPPFLAGS) $(CPPFLAGS) $(libcosi_la_CXXFLAGS) $(CXXFLAGS) $(libcosi_la_LDFLAGS) $(LDFLAGS) >> $@
	$(CXX) --version >> $@
if GITSTATE
	-git status $(srcdir) >> $@
	-git branch --no-abbrev -v
	-git diff $(srcdir)	>> $@
endif


all-local: $(COSIFLAGS)

MOSTLYCLEANFILES = $(COSIFLAGS)


