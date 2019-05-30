EXE = ../../build/bin/test_elas
MODEL = ../data/beam/beam.vtk
FIXED = ../data/beam/beam.csv
OUTDIR = ../results/beam
RUN_COMMAND = $(EXE) $(MODEL) $(FIXED) $(OUTDIR)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


