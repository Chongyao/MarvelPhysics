EXE = ../../build/bin/pipline_Z


MODEL_NAME = cube_hex
TYPE = hex

FIXED = ./data/$(MODEL_NAME)/$(MODEL_NAME).csv
MODEL = ./data/$(MODEL_NAME)/$(MODEL_NAME).vtk
OUTDIR = results/$(MODEL_NAME)

RUN_COMMAND = $(EXE) $(MODEL) $(FIXED) $(OUTDIR) $(TYPE)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


