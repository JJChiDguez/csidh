# INCLUDE FOLDER
BITLENGTH_OF_P?=512
BITS?=$(BITLENGTH_OF_P)
TYPE?=DUMMYFREE
INC_DIR+= -I./inc -I./inc/fp$(BITLENGTH_OF_P)/
# GLOBAL FLAGS
CFLAGS_ALWAYS?=-fcommon
# COMPILER
CC?=gcc-10

# REQUIRED FOR TESTS
FILES_REQUIRED_IN_CSIDH=./lib/rng.c \
			./lib/fp$(BITLENGTH_OF_P).S \
			./lib/point_arith.c ./lib/isogenies.c \
			./lib/action_simba_$(shell echo $(TYPE) | tr A-Z a-z).c \
			./main/csidh.c

OUTPUT_CSIDH=./bin/csidh
CFLAGS_CSIDH=-O3 -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -DFP_$(BITLENGTH_OF_P) -D$(TYPE)

FILES_REQUIRED_IN_CSIDH_UTIL=./lib/rng.c \
			./lib/fp$(BITLENGTH_OF_P).S \
			./lib/point_arith.c ./lib/isogenies.c \
			./lib/action_simba_$(shell echo $(TYPE) | tr A-Z a-z).c \
			./main/csidh_util.c
OUTPUT_CSIDH_UTIL=./bin/csidh-p$(BITS)-util
CFLAGS_CSIDH_UTIL=-O3 -funroll-loops -fomit-frame-pointer -m64 -mbmi2 -DFP_$(BITLENGTH_OF_P) -D$(TYPE) -DBITS=$(BITS)

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

help:
	@echo "\nusage: make csidh BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make csidh_util BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make util_test"
	@echo "usage: make regenerate_test_vectors"
	@echo "usage: make action_cost BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make action_timing BITLENGTH_OF_P=[512] TYPE=[WITHDUMMY_1/WITHDUMMY_2/DUMMYFREE]"
	@echo "usage: make clean\n"
	@echo "In addition, you can use an specific compiler by setting the variable CC with the "
	@echo "compiler name.\n\t\tCC=[any version of gcc compiler]"

util: csidh_util
csidh_util:
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_CSIDH_UTIL) -o $(OUTPUT_CSIDH_UTIL) $(CFLAGS_CSIDH_UTIL) $(CFLAGS_ALWAYS)

regenerate_test_vectors:
	./bin/csidh-p512-util -g -p sample-keys/1.montgomery.le.pk -s sample-keys/1.montgomery.le.sk
	./bin/csidh-p512-util -g -p sample-keys/2.montgomery.le.pk -s sample-keys/2.montgomery.le.sk
	./bin/csidh-p512-util -g -p sample-keys/3.montgomery.le.pk -s sample-keys/3.montgomery.le.sk
	./bin/csidh-p512-util -g -p sample-keys/4.montgomery.le.pk -s sample-keys/4.montgomery.le.sk
	./bin/csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/1.montgomery.le.sk > sample-keys/1-2.ss
	./bin/csidh-p512-util -d -p sample-keys/1.montgomery.le.pk -s sample-keys/2.montgomery.le.sk > sample-keys/2-1.ss
	./bin/csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/3.montgomery.le.sk > sample-keys/3-2.ss
	./bin/csidh-p512-util -d -p sample-keys/3.montgomery.le.pk -s sample-keys/4.montgomery.le.sk > sample-keys/4-3.ss

util_test: util
	echo "BEGIN util-test"
	./bin/csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/1.montgomery.le.sk > sample-keys/1-2.ss.test_result
	diff sample-keys/1-2.ss.test_result sample-keys/1-2.ss
	./bin/csidh-p512-util -d -p sample-keys/1.montgomery.le.pk -s sample-keys/2.montgomery.le.sk > sample-keys/2-1.ss.test_result
	diff sample-keys/2-1.ss.test_result sample-keys/2-1.ss
	./bin/csidh-p512-util -d -p sample-keys/2.montgomery.le.pk -s sample-keys/3.montgomery.le.sk > sample-keys/3-2.ss.test_result
	diff sample-keys/3-2.ss.test_result sample-keys/3-2.ss
	./bin/csidh-p512-util -d -p sample-keys/3.montgomery.le.pk -s sample-keys/4.montgomery.le.sk > sample-keys/4-3.ss.test_result
	diff sample-keys/4-3.ss.test_result sample-keys/4-3.ss
	rm sample-keys/*.test_result
	echo "END util-test"

csidh:
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_CSIDH) -o $(OUTPUT_CSIDH) $(CFLAGS_CSIDH) $(CFLAGS_ALWAYS)

action_cost: 
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_ACTION) -o $(OUTPUT_ACTION) $(CFLAGS_ACTION) $(CFLAGS_ALWAYS)

action_timing: 
	$(CC) $(INC_DIR) $(FILES_REQUIRED_IN_ACTION_CC) -o $(OUTPUT_ACTION_CC) $(CFLAGS_ACTION_CC) $(CFLAGS_ALWAYS)

clean:
	rm -f ./bin/* sample-keys/*.test_result

