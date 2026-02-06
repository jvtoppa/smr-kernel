# Compiler flags
OPTIONS = -O3 -ffast-math -std=c++17
SEQAN_INCLUDE = /home/joao/Documentos/seqan/include
SDSL_INCLUDE = sdsl-lite/include
SDSL_LIB = sdsl-lite/lib

# Ensure build directory exists
BUILD_DIR = build
$(shell mkdir -p $(BUILD_DIR))

# Targets
all: sdsl main gen_test kernelize

debug: OPTIONS = -O0 -g -DDEBUG
debug: all

# Build SDSL
sdsl:
	@if [ ! -d "sdsl-lite" ]; then \
		git clone https://github.com/simongog/sdsl-lite.git; \
	fi
	@if [ ! -f "$(SDSL_LIB)/libsdsl.a" ]; then \
		cd sdsl-lite && ./install.sh $(PWD)/sdsl-lite; \
	fi

# libsais object
$(BUILD_DIR)/libsais.o: libsais/src/libsais.c | $(BUILD_DIR)
	gcc -O3 -c -I libsais/include $< -o $@

# Main executable
main: sdsl $(BUILD_DIR)/libsais.o src/main.cpp | $(BUILD_DIR)
	g++ $(OPTIONS) -I$(SEQAN_INCLUDE) -I$(SDSL_INCLUDE) -Ilibsais/include \
		src/main.cpp $(BUILD_DIR)/libsais.o -o $(BUILD_DIR)/smr \
		-L$(SDSL_LIB) -lsdsl -ldivsufsort -ldivsufsort64

# Test generator
gen_test: src/gen_test.cpp | $(BUILD_DIR)
	g++ $(OPTIONS) -I$(SEQAN_INCLUDE) \
		src/gen_test.cpp -o $(BUILD_DIR)/gen_test

# Kernelizer
kernelize: sdsl $(BUILD_DIR)/libsais.o src/kernelize.cpp | $(BUILD_DIR)
	g++ $(OPTIONS) -I$(SEQAN_INCLUDE) -I$(SDSL_INCLUDE) -Ilibsais/include \
		src/kernelize.cpp $(BUILD_DIR)/libsais.o -o $(BUILD_DIR)/kernelize \
		-L$(SDSL_LIB) -lsdsl -ldivsufsort -ldivsufsort64

# Phony targets
.PHONY: all debug clean sdsl

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)/*
