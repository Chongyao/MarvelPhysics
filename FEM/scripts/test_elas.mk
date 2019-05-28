EXE = ../../build/bin/test_elas
MODEL = ../data/beam/beam.vtk
FIXED = ../data/beam/beam.csv
OUTDIR = ../results/beam

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(EXE) $(MODEL) $(FIXED) $(OUTDIR)


