# INCLUDE FOLDER
INC_DIR+= -I./inc -I./inc/fp$(BITLENGTH_OF_P)/

# REQUIRED FOR TESTS
FILES_REQUIRED_IN_CSIDH=./lib/rng.c \
			./lib/fp$(BITLENGTH_OF_P).S \
			./lib/point_arith.c ./lib/isogenies.c \
			./lib/action_simba_$(shell echo $(TYPE) | tr A-Z a-z).c \
			./main/csidh.c

OUTPUT_CSIDH=./bin/csidh
CFLAGS_CSIDH=-O3 -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -DFP_$(BITLENGTH_OF_P) -D$(TYPE)

# REQUIRED FOR COSTS
FILES_REQUIRED_IN_ACTION=./lib/rng.c \
			./lib/fp$(BITLENGTH_OF_P).S \
			./lib/point_arith.c ./lib/isogenies.c \
			./lib/action_simba_$(shell echo $(TYPE) | tr A-Z a-z).c \
			./main/action_cost.c

OUTPUT_ACTION=./bin/action_cost
CFLAGS_ACTION=-O3 -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -DFP_$(BITLENGTH_OF_P) -D$(TYPE) -lm

# REQUIRED FOR CLOCK CYCLES
FILES_REQUIRED_IN_ACTION_CC=./lib/rng.c \
			./lib/fp$(BITLENGTH_OF_P).S \
			./lib/point_arith.c ./lib/isogenies.c \
			./lib/action_simba_$(shell echo $(TYPE) | tr A-Z a-z).c \
			./main/action_timing.c

OUTPUT_ACTION_CC=./bin/action_timing
CFLAGS_ACTION_CC=-O3 -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -DFP_$(BITLENGTH_OF_P) -D$(TYPE) -lm

# GLOBAL FLAGS
CFLAGS_ALWAYS=

# COMPILER
CC=gcc

help:
	@echo "\nusage: make csidh BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make action_cost BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make action_timing BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make clean\n"
	@echo "In addition, you can use an specific compiler by setting the variable CC with the "
	@echo "compiler name.\n\t\tCC=[any version of gcc compiler]"

csidh:
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_CSIDH) -o $(OUTPUT_CSIDH) $(CFLAGS_CSIDH) $(CFLAGS_ALWAYS)

action_cost: 
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_ACTION) -o $(OUTPUT_ACTION) $(CFLAGS_ACTION) $(CFLAGS_ALWAYS)

action_timing: 
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_ACTION_CC) -o $(OUTPUT_ACTION_CC) $(CFLAGS_ACTION_CC) $(CFLAGS_ALWAYS)

clean:
	rm -f ./bin/*


