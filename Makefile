FC = mpif90
FLAGS = -fallow-argument-mismatch -O3 -ffast-math -march=native -fopenmp
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
OUT_DIR = ./outs
EXAMPLES_DIR = ./examples/

MODULE_SOURCES = $(SRC_DIR)/particle_class.f90   
MODULE_OBJECTS = $(MODULE_SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

SOURCES = $(filter-out $(MODULE_SOURCES), $(wildcard $(SRC_DIR)/*.f90))  
OBJECTS = $(SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

directories:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OUT_DIR)

mpi: directories $(BIN_DIR)/program_mpi

$(BIN_DIR)/program_mpi: $(EXAMPLES_DIR)/example_mpi.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

mpigrid: directories $(BIN_DIR)/program_mpigrid

$(BIN_DIR)/program_mpigrid: $(EXAMPLES_DIR)/example_mpigrid.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

grid: directories $(BIN_DIR)/program_grid

$(BIN_DIR)/program_grid: $(EXAMPLES_DIR)/example_grid.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

omp: directories $(BIN_DIR)/program_omp

$(BIN_DIR)/program_omp: $(EXAMPLES_DIR)/example_omp.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

tt: directories $(BIN_DIR)/program_twotime

$(BIN_DIR)/program_twotime: $(EXAMPLES_DIR)/twotime_mpi.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

dtime: directories $(BIN_DIR)/program_dtime

$(BIN_DIR)/program_dtime: $(EXAMPLES_DIR)/D_time.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

mtime: directories $(BIN_DIR)/program_mtime

$(BIN_DIR)/program_mtime: $(EXAMPLES_DIR)/multitime.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@

infer: directories $(BIN_DIR)/program_infer

$(BIN_DIR)/program_infer: $(EXAMPLES_DIR)/inference.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@
.SILENT: clean
clean:
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(BIN_DIR)/program*
