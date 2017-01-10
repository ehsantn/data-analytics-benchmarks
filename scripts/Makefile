CC=mpiicpc
CFLAGS=-std=c++11 -O3 -wd2593
LDFLAGS=

CPP_DIR=../src/main/cpp/

all: linear_regression_gen logistic_regression_gen kmeans_gen kernelscore_gen

linear_regression_gen: $(CPP_DIR)/linear_regression_gen_main.cc $(CPP_DIR)/linear_regression_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

logistic_regression_gen: $(CPP_DIR)/logistic_regression_gen_main.cc $(CPP_DIR)/logistic_regression_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

kmeans_gen: $(CPP_DIR)/kmeans_gen_main.cc $(CPP_DIR)/kmeans_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

kernelscore_gen: $(CPP_DIR)/kernelscore_gen_main.cc $(CPP_DIR)/kernelscore_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

