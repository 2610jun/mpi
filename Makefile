# Compiler and flags
CXX = mpicxx
CXXFLAGS = -std=c++11 -O2

# Files
TARGET = ecc
SRCS = ecc.cpp polar.cpp
RM = rm -f

# Run parameters
NUM_PROCESSORS = 200
# N, K, ENCODER, DECODER, CRC_LEN, ...

# Simulation parameter

RM_METHOD = POLAR_5G_STANDARD
#RM_METHOD = HUAWEI_APPROX
ENCODER_TYPE = POLAR
DECODER_TYPE = CA_SCL

N = 128
K = 64
NUM_LIST = 256
CRC_TYPE = 7

EBN0_BEGIN = 0
EBN0_END = 3
EBN0_INTERVAL = 0.5

N_MAX_BLK = 10000000
TARGET_FRAME_ERR = 250



# Default target
all:
	$(CXX) $(CXXFLAGS) -o  $(TARGET) $(SRCS)

# Run target
run: all
	mpiexec -f ojh_hostfile -n $(NUM_PROCESSORS) ./$(TARGET) $(ENCODER_TYPE) $(DECODER_TYPE) $(N) $(K) $(NUM_LIST) $(CRC_TYPE) $(EBN0_BEGIN) $(EBN0_END) $(EBN0_INTERVAL) $(N_MAX_BLK) $(TARGET_FRAME_ERR) $(RM_METHOD)
# Clean target
clean:
	$(RM) $(TARGET)
