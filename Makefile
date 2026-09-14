CC = cc
PYTHON = python3
CPPFLAGS += -Iinclude
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
LDLIBS += -lm

.PHONY: all test clean
all: build/alcubierre_sim

build:
	mkdir -p build

build/alcubierre.o: src/alcubierre.c include/alcubierre.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

build/main.o: src/main.c include/alcubierre.h | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

build/alcubierre_sim: build/main.o build/alcubierre.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/test_model: tests/test_model.c build/alcubierre.o include/alcubierre.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) tests/test_model.c build/alcubierre.o $(LDLIBS) -o $@

test: build/test_model build/alcubierre_sim
	./build/test_model
	$(PYTHON) tests/test_cli.py ./build/alcubierre_sim

clean:
	rm -rf build
