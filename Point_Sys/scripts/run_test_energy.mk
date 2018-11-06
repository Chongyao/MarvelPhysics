EXE = ../../build/bin/test_energy

model_name = 53092_sf-G_G0_subdivide937772130937772130
IN_DIR = ../../Point_Sys/data
OUT_DIR = ../../Point_Sys/result
poi = 0.45
Young = 100

surf = $(IN_DIR)/$(model_name).obj
points_out = $(OUT_DIR)/$(model_name).vtk
num_in_axis=5

test_genpoints: $(EXE)
	$(EXE) surf=$(surf) points_out=$(points_out) num_in_axis=$(num_in_axis) rho=5 Poission=$(poi) Young=$(Young) res=$(OUT_DIR)/$(model_name) ker_cof=8

