
# Do not modify this file. It was generated by mvme_root_client and will be
# overwritten when next running the client.

ROOTCFLAGS	=	$(shell root-config --cflags)
ROOTLIBS	=	$(shell root-config --libs)
MVMELIBS        =       -L$(PWD)  -lmvme_root_event
#INCLUDES    =	

CXXFLAGS   +=	-std=c++14 -fPIC
CXXFLAGS   +=   $(INCLUDES)

EXP_LIB		= libExperiment_mvme.so
EXP_DICT	= libExperiment_mvme_dict.cxx
EXP_LINK    = -L. -lExperiment_mvme
EXP_LINKDEF = Experiment_mvme_LinkDef.h

EVENT_OBJECTS_HEADERS = $(current_dir)mvme_root_event_objects.h $(current_dir)mvme_root_event_objects_LinkDef.h
EVENT_OBJECTS_DEPS = $(current_dir)mvme_root_event_objects.cxx $(EVENT_OBJECTS_HEADERS)

ALL_TARGETS	= libmvme_root_event.so $(EXP_LIB) mvmesort
EXTRA_CLEAN =

.DEFAULT_GOAL = all
.PHONY: all clean

$(EXP_DICT): Experiment_mvme.cxx Experiment_mvme.h $(EXP_LINKDEF)
	rootcling -f $@ -rml $(EXP_LIB) -rmf libExperiment_mvme.rootmap $(INCLUDES) Experiment_mvme.h $(EXP_LINKDEF)

$(EXP_LIB): Experiment_mvme.cxx Experiment_mvme.h $(EXP_DICT)
	$(CXX) $(LDFLAGS) $(ROOTCFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(MVMELIBS) -shared -fPIC -o $@ $< $(EXP_DICT)

mvme_root_event_rdict.cxx: $(EVENT_OBJECTS_DEPS)
	rootcling -f $@ -rml libmvme_root_event.so -rmf mvme_root_event.rootmap $(EVENT_OBJECTS_HEADERS)

libmvme_root_event.so: $(EVENT_OBJECTS_DEPS) mvme_root_event_rdict.cxx
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) $(ROOTLIBS) -shared -fPIC -o $@ $< mvme_root_event_rdict.cxx

MVMESort.o: MVMESort.cxx MVMESort.h 
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -c $<

mvmesort: mvmesort.cxx $(EXP_LIB) MVMESort.o
	$(CXX) $(LDFLAGS) $(ROOTCFLAGS) $(CXXFLAGS) MVMESort.o $< '-Wl,-rpath,$$ORIGIN/:$$ORIGIN/../lib' $(EXP_LINK) -lmvme_root_event $(ROOTLIBS) -o $@

all: $(ALL_TARGETS)

clean:
	-rm -f $(EXP_LIB) $(EXP_DICT) $(EXTRA_CLEAN) *.so *.o *.pcm *.rootmap mvme_root_event_rdict.cxx
