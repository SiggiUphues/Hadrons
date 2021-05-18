AM_LDFLAGS += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)

bin_PROGRAMS = \
  MDWF_TCorr \
	MDWF_SCorr

MDWF_TCorr_SOURCES = MDWF_TCorr.cpp
MDWF_TCorr_LDADD   = -lHadrons -lGrid

MDWF_SCorr_SOURCES = MDWF_SCorr.cpp
MDWF_SCorr_LDADD   = -lHadrons -lGrid
