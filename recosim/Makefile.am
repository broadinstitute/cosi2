

bin_PROGRAMS += recosim get_recomap recomap_hapmap2

BOOSTRAND=$(srcdir)/ext/boost/libs/random/src/random_device.cpp $(srcdir)/ext/boost/boost/random/random_device.hpp

recosim_SOURCES = $(srcdir)/recosim/recosim.cc $(srcdir)/recosim/gamma.hpp $(BOOSTRAND)
recosim_CPPFLAGS = -I$(srcdir) -I$(srcdir)/ext/boost -DBOOST_ERROR_CODE_HEADER_ONLY
recosim_LDFLAGS = -lm

# The programs below are used to extract a genetic map from a random region of the
# 1000 Genomes genetic map.

get_recomap_SOURCES = $(srcdir)/recosim/get_recomap.cc $(BOOSTRAND)
get_recomap_CPPFLAGS = -I$(srcdir)/ext/boost -DBOOST_ERROR_CODE_HEADER_ONLY
get_recomap_LDFLAGS = -lm

recomap_hapmap2_SOURCES = $(srcdir)/recosim/recomap_hapmap2.cc $(BOOSTRAND)
recomap_hapmap2_CPPFLAGS = -I$(srcdir)/ext/boost -DBOOST_ERROR_CODE_HEADER_ONLY
recomap_hapmap2_LDFLAGS = -lm






