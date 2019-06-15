EXE = ../../build/bin/test_elas
MODEL = ../data/beam/beam.vtk
FIXED = ../data/beam/beam.csv
OUTDIR = ../results/beam
ROTATE = ../data/beam/rotate.csv


RUN_COMMAND = $(EXE) $(MODEL) $(FIXED) $(OUTDIR) $(ROTATE)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


