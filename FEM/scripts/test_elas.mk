EXE = ../../build/bin/test_elas

MODEL_NAME = beam
MODEL = ../data/$(MODEL_NAME)/$(MODEL_NAME).vtk
FIXED = ../data/$(MODEL_NAME)/$(MODEL_NAME).csv
OUTDIR = ../results/$(MODEL_NAME)



RUN_COMMAND = $(EXE) $(MODEL) $(FIXED) $(OUTDIR)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


