VERSION = 3.02
CC      = /usr/bin/g++
CFLAGS  = -Wall -DVERSION=\"$(VERSION)\" -DJULIA_ENABLE_THREADING=1

JL_SHARE = $(shell julia -e 'print(joinpath(JULIA_HOME,Base.DATAROOTDIR,"julia"))')
CFLAGS   += $(shell $(JL_SHARE)/julia-config.jl --cflags)
CXXFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)

SRC_DIR = ./src

PACKAGE = $(SRC_DIR)/core/QMDDpackage.o $(SRC_DIR)/core/qcost.o $(SRC_DIR)/core/textFileUtilities.o $(SRC_DIR)/core/timing.o $(SRC_DIR)/core/QMDDcircuit.o $(SRC_DIR)/core/QMDDcomplexD.o  $(SRC_DIR)/core/QMDDreorder.o

all: qmdd_cmd qmdd_eq_check qmdd_minimize

qmdd_cmd: package $(SRC_DIR)/applications/qmdd_cmd.o 
	$(CC) $(CFLAGS) -o qmdd_cmd $(SRC_DIR)/applications/qmdd_cmd.o $(PACKAGE) $(LDFLAGS) $(LDLIBS)

qmdd_eq_check: package $(SRC_DIR)/applications/qmdd_equivalence_check.o
	$(CC) $(CFLAGS) -o qmdd_equivalence_check $(SRC_DIR)/applications/qmdd_equivalence_check.o $(PACKAGE) $(LDFLAGS) $(LDLIBS)

qmdd_minimize: package $(SRC_DIR)/applications/qmdd_minimize.o
	$(CC) $(CFLAGS) -o qmdd_minimize $(SRC_DIR)/applications/qmdd_minimize.o $(PACKAGE) $(LDFLAGS) $(LDLIBS)
	
test3: package test3.o
	$(CC) $(CFLAGS) -o test3 test3.o $(PACKAGE) $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

package: $(PACKAGE)
	ar -q libQMDD.a $(PACKAGE)
	
clean:
	rm -rf $(SRC_DIR)/core/*.o
	rm -rf $(SRC_DIR)/applications/*.o
	rm -rf qmdd_*
	rm -rf libQMDD.a
