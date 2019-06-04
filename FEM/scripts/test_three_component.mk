EXE = ../../build/bin/test_three_component
# MODEL = ../data/beam/beam.vtk
# FIXED = ../data/beam/beam.csv
# OUTDIR = ../results/beam
RUN_COMMAND = $(EXE) # $(MODEL) $(FIXED) $(OUTDIR)

run: $(EXE)
	$(RUN_COMMAND)	


