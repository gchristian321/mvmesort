

ANALYSIS_LIBS =
ANALYSIS_DEPS = $(EXP_LIB)
ALL_TARGETS += analysis.so
EXTRA_CLEAN = analysis.so analysis.o

analysis.o: analysis.cxx $(ANALYSIS_DEPS)
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -c $<

analysis.so: analysis.o
	$(CXX) $(LDFLAGS) $(ROOTLIBS) $(MVMELIBS) $(ANALYSIS_LIBS) -Wl,-rpath,. $(EXP_LINK) -shared -o $@ $<

