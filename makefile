SRC_DIR := src
SUPER_DIR := $(SRC_DIR)/super
FC := gfortran
NBITS ?= 64
DLEXT ?= so
LIB_DIR ?= .

FFLAGS := -O3 -fPIC -fopenmp
ifeq ($(NBITS), 64)
    FFLAGS += -fdefault-integer-8
endif

OBJS := $(SRC_DIR)/cfs.o $(SRC_DIR)/bs.o $(SRC_DIR)/op.o $(SRC_DIR)/diag.o $(SRC_DIR)/diag_re.o $(SRC_DIR)/ent.o
OBJS_SUPER := $(SUPER_DIR)/scfs.o $(SUPER_DIR)/sbs.o $(SUPER_DIR)/sop.o $(SUPER_DIR)/sent.o
MODS := cfs.mod bs.mod op.mod diag.mod diag_re.mod ent.mod scfs.mod sbs.mod sop.mod sent.mod

.PHONY: all clean cleantmp libfuzzified libfuzzifino

all: libfuzzified libfuzzifino cleantmp

$(SRC_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(SUPER_DIR)/%.o: $(SUPER_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

libfuzzified: $(OBJS)
	$(FC) $(FFLAGS) -shared -o $(LIB_DIR)/libfuzzified.$(DLEXT) $^ -L$(LIB_DIR) -larpack

libfuzzifino: $(OBJS_SUPER)
	$(FC) $(FFLAGS) -shared -o $(LIB_DIR)/libfuzzifino.$(DLEXT) $^

cleantmp:
	rm -f $(OBJS) $(OBJS_SUPER) $(MODS)

clean:
	rm -f $(OBJS) $(OBJS_SUPER) $(MODS) $(LIB_DIR)/libfuzzified.$(DLEXT) $(LIB_DIR)/libfuzzifino.$(DLEXT)
