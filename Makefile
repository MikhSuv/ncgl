# Build the demo program and the tests without fpm.
#
#   make        build the demo program ./nclg
#   make test   build the test programs and run them
#   make clean  remove the object files, the module files and the programs
#
# Everything is built in the current directory, so no -J and no -I are needed:
# gfortran writes the .mod files where it is called from and finds them there.

FC     := gfortran
FFLAGS := -O3 -DNDEBUG

# Library modules, in the order the compiler has to see them.
LIB_SRC := src/precision_mod.f90 src/function_interfaces.f90 \
           src/test_functions.f90 src/nclg.f90
LIB_OBJ := $(patsubst src/%.f90,%.o,$(LIB_SRC))

# The program built by the default target. Its object is main.o: the name nclg
# belongs to the module of the same name.
APP     := nclg
APP_OBJ := main.o

# The programs run by "make test".
TESTS    := test_grad test_cg_min
TEST_OBJ := $(TESTS:%=%.o)

ALL_OBJ := $(LIB_OBJ) $(APP_OBJ) $(TEST_OBJ)

# The sources are looked up by file name in these three directories.
vpath %.f90 src test app

.PHONY: all test clean
.DEFAULT_GOAL := all

all: $(APP)

# The library is compiled in a single call, so that the modules end up in the
# order they use each other. Its objects are all older than any change of the
# library, which is what makes the programs below be recompiled as well.
$(LIB_OBJ): $(LIB_SRC)
	$(FC) $(FFLAGS) -c $(LIB_SRC)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# A program needs the module files of the library, so its object is compiled
# after the objects that write them.
$(APP_OBJ) $(TEST_OBJ): $(LIB_OBJ)

$(APP): $(APP_OBJ) $(LIB_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(TESTS): %: %.o $(LIB_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# Run every test program. The whole loop runs in one shell, so the exit status
# of the last failure is the exit status of the target.
test: $(TESTS)
	@status=0; \
	for t in $(TESTS); do \
	  echo "--- ./$$t ---"; \
	  ./$$t || status=1; \
	done; \
	exit $$status

clean:
	rm -f $(ALL_OBJ) *.mod $(APP) $(TESTS)
