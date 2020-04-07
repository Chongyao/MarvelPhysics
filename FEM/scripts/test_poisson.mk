EXE = ../../build/bin/test_poisson

MODEL_NAME = cube_hex
TYPE = hex
MODEL = ../data/$(MODEL_NAME)/$(MODEL_NAME).vtk
TOP_FIXED = ../data/$(MODEL_NAME)/$(MODEL_NAME)_top.csv
BOTTOM_FIXED = ../data/$(MODEL_NAME)/$(MODEL_NAME)_bottom.csv
OUTDIR = ../results/POI/$(MODEL_NAME)




RUN_COMMAND = $(EXE) $(MODEL) $(TOP_FIXED) $(BOTTOM_FIXED) $(OUTDIR) $(TYPE)

run: $(FIXED) $(MODEL) $(EXE)
	@mkdir -p $(OUTDIR)
	$(RUN_COMMAND)


