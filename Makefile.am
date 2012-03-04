INCLUDES = -Wall -g
LFLAGS = `gsl-config --libs` -g

include_HEADERS = \
DescriptiveStatistics.h  DiscreteDistribution.h  FileReader.h  \
LowCoverageFinder.h  MappingReader.h  \
SAMAlignment.h GeneralFeature.h \
CoverageArray.h SAMReader.h CoverageDistribution.h EvidenceFinder.h \
SequenceSet.h \
Tool.h Option.h  Utility.h

bin_PROGRAMS = chopsticks

chopsticks_SOURCES = \
chopsticks.cc \
Option.cc Utility.cc FileReader.cc SAMAlignment.cc SequenceSet.cc GeneralFeature.cc \
SAMReader.cc CoverageArray.cc CoverageDistribution.cc EvidenceFinder.cc \
Tool.cc
chopsticks_LDFLAGS = $(LFLAGS)
