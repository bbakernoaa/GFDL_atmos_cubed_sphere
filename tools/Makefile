# Makefile for FV3 Tracer Transport Optimization Tests

# Compiler settings
FC = mpif90
FFLAGS = -O3 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
MPIFC = mpif90

# Directories
SRCDIR = ../model
TOOLDIR = .

# Source files
SOURCES = $(SRCDIR)/fv_arrays.F90 \
          $(SRCDIR)/fv_grid_utils.F90 \
          $(SRCDIR)/fv_tracer2d.F90 \
          $(SRCDIR)/fv_tracer2d_optimized.F90 \
          $(SRCDIR)/tp_core.F90 \
          $(SRCDIR)/tp_core_optimized.F90 \
          $(TOOLDIR)/performance_measurement.F90 \
          $(TOOLDIR)/tracer_test_suite.F90 \
          $(TOOLDIR)/benchmark_test_suite.F90 \
          $(TOOLDIR)/regression_test_suite.F90 \
          $(TOOLDIR)/tracer_transport_integrator.F90

# Object files
OBJECTS = $(SOURCES:.F90=.o)

# Executables
TEST_EXE = run_comprehensive_tests.exe
BENCH_EXE = run_benchmark_tests.exe
REPORT_EXE = generate_performance_report.exe

# Default target
all: $(TEST_EXE) $(BENCH_EXE) $(REPORT_EXE)

# Compile the comprehensive test executable
$(TEST_EXE): $(OBJECTS) $(TOOLDIR)/run_comprehensive_tests.F90
	$(MPIFC) $(FFLAGS) -o $(TEST_EXE) $(TOOLDIR)/run_comprehensive_tests.F90 $(OBJECTS)

# Compile the benchmark test executable
$(BENCH_EXE): $(OBJECTS) $(TOOLDIR)/run_benchmark_tests.F90
	$(MPIFC) $(FFLAGS) -o $(BENCH_EXE) $(TOOLDIR)/run_benchmark_tests.F90 $(OBJECTS)

# Compile the performance report executable
$(REPORT_EXE): $(OBJECTS) $(TOOLDIR)/generate_performance_report.F90
	$(MPIFC) $(FFLAGS) -o $(REPORT_EXE) $(TOOLDIR)/generate_performance_report.F90 $(OBJECTS)

# Compile source files to object files
%.o: %.F90
	$(MPIFC) $(FFLAGS) -c $< -o $@

# Clean object files
clean_obj:
	rm -f *.o

# Clean executables
clean_exe:
	rm -f $(TEST_EXE) $(BENCH_EXE) $(REPORT_EXE)

# Clean all
clean: clean_obj clean_exe
	rm -f *.mod

# Run tests
test: $(TEST_EXE)
	./$(TEST_EXE)

# Run benchmarks
benchmark: $(BENCH_EXE)
	./$(BENCH_EXE)

# Generate report
report: $(REPORT_EXE)
	./$(REPORT_EXE)

.PHONY: all clean clean_obj clean_exe test benchmark report