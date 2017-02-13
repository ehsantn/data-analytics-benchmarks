CC=mpiicpc
CFLAGS=-std=c++11 -O3 -xHost -wd2593
LDFLAGS=

CPP_DIR=../src/main/cpp/

all: linear_regression_gen logistic_regression logistic_regression_gen kmeans_gen kernelscore_gen

manual: logistic_regression_manual_gen logistic_regression_manual kernelscore_gen_manual linear_regression_manual_gen kmeans_manual_gen

linear_regression_gen: $(CPP_DIR)/linear_regression_gen_main.cc $(CPP_DIR)/linear_regression_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

logistic_regression_gen: $(CPP_DIR)/logistic_regression_gen_main.cc $(CPP_DIR)/logistic_regression_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

logistic_regression: $(CPP_DIR)/logistic_regression_main.cc $(CPP_DIR)/logistic_regression.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< -lhdf5

kmeans_gen: $(CPP_DIR)/kmeans_gen_main.cc $(CPP_DIR)/kmeans_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

kernelscore_gen: $(CPP_DIR)/kernelscore_gen_main.cc $(CPP_DIR)/kernelscore_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

logistic_regression_manual_gen: $(CPP_DIR)/logistic_regression_manual_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

logistic_regression_manual: $(CPP_DIR)/logistic_regression_manual.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< -lhdf5

kernelscore_gen_manual: $(CPP_DIR)/kernelscore_gen_manual.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

linear_regression_manual_gen: $(CPP_DIR)/linear_regression_manual_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

kmeans_manual_gen: $(CPP_DIR)/kmeans_manual_gen.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

