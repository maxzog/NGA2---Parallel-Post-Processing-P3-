FC = gfortran
FLAGS = -fallow-argument-mismatch -O3 -ffast-math -fopenmp
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
EXAMPLES_DIR = ./examples/
MODULE_SOURCES = $(SRC_DIR)/particle_class.f90   
MODULE_OBJECTS = $(MODULE_SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
SOURCES = $(filter-out $(MODULE_SOURCES), $(wildcard $(SRC_DIR)/*.f90))  
OBJECTS = $(SOURCES:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
all: $(BIN_DIR)/program
$(BIN_DIR)/program: $(EXAMPLES_DIR)/temp_omp.f90 $(MODULE_OBJECTS) $(OBJECTS)  
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $^
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -c $< -o $@
.SILENT: clean
clean:
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(BIN_DIR)/program
