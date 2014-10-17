PVFMM_DIR=./../../ext/pvfmm
-include $(PVFMM_DIR)/MakeVariables

ifndef CXXFLAGS_PVFMM
$(error Cannot find file: MakeVariables)
endif

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_BIN = \
       $(BINDIR)/semilag_tree

CXXFLAGS := $(CXXFLAGS_PVFMM) -std=c++0x

all : $(TARGET_BIN)

ifeq ($(INTEL_OFFLOAD_OK),yes)

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS) -no-offload         $^       $(LDFLAGS_PVFMM) -o $@
	$(CXX_PVFMM) $(CXXFLAGS)                     $^_async $(LDFLAGS_PVFMM) -o $@_async
	$(CXX_PVFMM) $(CXXFLAGS) -D__DEVICE_SYNC__=1 $^_mic   $(LDFLAGS_PVFMM) -o $@_mic

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS) -no-offload         -I$(INCDIR) -c $^ -o $@
	$(CXX_PVFMM) $(CXXFLAGS)                     -I$(INCDIR) -c $^ -o $@_async
	$(CXX_PVFMM) $(CXXFLAGS) -D__DEVICE_SYNC__=1 -I$(INCDIR) -c $^ -o $@_mic

else

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS)                  $^       $(LDFLAGS_PVFMM) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS)                  -I$(INCDIR) -c $^ -o $@

endif

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~
