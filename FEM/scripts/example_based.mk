EXE = ../../build/bin/example_based
MODEL = ../data/beam/beam.vtk
FIXED = ../data/beam/beam.csv
OUTDIR = ../results/example_based/beam_origin
EXAMPLE = ../data/beam/example.vtk


RUN_COMMAND = $(EXE) $(MODEL) $(FIXED) $(OUTDIR) $(EXAMPLE)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


