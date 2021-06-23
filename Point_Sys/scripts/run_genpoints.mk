EXE = ../../build/bin/test_genpoints

model_name = dragon
IN_DIR = ../../Point_Sys/data
OUT_DIR = ../../Point_Sys/result


surf = $(IN_DIR)/$(model_name).obj
points_out = $(OUT_DIR)/$(model_name).vtk
num_in_axis=30
test_genpoints: $(EXE)
	$(EXE) surf=$(surf) points_out=$(points_out) num_in_axis=$(num_in_axis)

