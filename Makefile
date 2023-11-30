FC = gfortran    
FLAGS = -O2
LD_FLAGS = 
SRC_DIR = ./src
EXAMPLES_DIR = ./examples
OBJ_DIR = ./obj
BIN_DIR = ./bin

OBJ_FILES = $(addprefix $(OBJ_DIR)/, $(notdir $(patsubst %.f90,%.o,$(wildcard $(SRC_DIR)/*.f90))))
TEST_DRIVER = $(BIN_DIR)/testing

all: directories $(TEST_DRIVER)

directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

$(TEST_DRIVER): $(EXAMPLES_DIR)/temp_driver.f90 $(OBJ_FILES)
	$(FC) $(FLAGS) -I$(SRC_DIR) -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) -J$(OBJ_DIR) -o $@ -c $<

clean:
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
