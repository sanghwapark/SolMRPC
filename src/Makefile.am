AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(ROOTSYS)/include

AM_LDFLAGS = \
	-L$(libdir) \
	-L$(ROOTSYS)/lib

lib_LTLIBRARIES = \
  libsolmrpc.la

libsolmrpc_la_SOURCES = \
  utlSolidBase.C \
  SolMRPCdigit.C \
  SolMRPCStrip.C

libsolmrpc_la_LDFLAGS = \
	-L$(ROOTSYS)/lib `root-config --libs`

include_HEADERs = \
  utlSolidBase.h \
  SolidConst.h \
  SolMRPCdigit.h \
  SolMRPCStrip.h

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libsolmrpc.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# This is where the dictionary file is generated
%_Dict.C: \
	utlSolidBase.h \
        SolMRPCdigit.h \
	SolMRPCStrip.h \
	utlSolidBaseLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
