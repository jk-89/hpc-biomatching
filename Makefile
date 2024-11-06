CXX := g++
NVCC := /usr/local/cuda/bin/nvcc
CXXFLAGS := -O3 -Wall -std=c++20
NVCCFLAGS := -O3 -std=c++20

SRC_DIR := src
OBJ_DIR := obj

COMMON_CPP_FILES := $(SRC_DIR)/common.cpp
CPU_SRC_FILES := $(SRC_DIR)/cpugenv.cpp
GPU_SRC_FILES := $(SRC_DIR)/gpugenv.cu
GPU_CONTINUOUS_SRC_FILES := $(SRC_DIR)/gpugenv_continuous.cu
GPU_STREAM_SRC_FILES := $(SRC_DIR)/gpugenv_stream.cu

COMMON_OBJ_FILES := $(COMMON_CPP_FILES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
CPU_OBJ_FILES := $(CPU_SRC_FILES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
GPU_OBJ_FILES := $(GPU_SRC_FILES:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)
GPU_CONTINUOUS_OBJ_FILES := $(GPU_CONTINUOUS_SRC_FILES:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)
GPU_STREAM_OBJ_FILES := $(GPU_STREAM_SRC_FILES:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)

CPU_EXE := cpugenv
GPU_EXE := gpugenv
GPU_CONTINUOUS_EXE := gpugenv_continuous
GPU_STREAM_EXE := gpugenv_stream

$(shell mkdir -p $(OBJ_DIR))

all: $(CPU_EXE) $(GPU_EXE) $(GPU_CONTINUOUS_EXE) $(GPU_STREAM_EXE)

$(CPU_EXE): $(COMMON_OBJ_FILES) $(CPU_OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(GPU_EXE): $(COMMON_OBJ_FILES) $(GPU_OBJ_FILES)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(GPU_CONTINUOUS_EXE): $(COMMON_OBJ_FILES) $(GPU_CONTINUOUS_OBJ_FILES)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(GPU_STREAM_EXE): $(COMMON_OBJ_FILES) $(GPU_STREAM_OBJ_FILES)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

clean:
	rm -f $(COMMON_OBJ_FILES) $(CPU_OBJ_FILES) $(GPU_OBJ_FILES) $(GPU_CONTINUOUS_OBJ_FILES) $(GPU_STREAM_OBJ_FILES) $(CPU_EXE) $(GPU_EXE) $(GPU_CONTINUOUS_EXE) $(GPU_STREAM_EXE)
	rm -rf $(OBJ_DIR)
