# 后缀规则
.SUFFIXES :
.SUFFIXES : .o .cpp

# 编译器选择
CXX = mpic++    # 使用 OpenMPI 进行 MPI 编译
#CXX = mpicxx    # 可选 MPICH
#CXX = mpiicpx   # 可选 Intel OneAPI
#CXX = icpx      # 可选 Intel C++ 编译器
#CXX = g++       # 可选 GNU C++ 编译器

# 输出目标
PROG = ./makinput.exe

# 编译标志
CXXFLAGS = -O1 -std=c++17 -I/home/pan/petsc-3.24.5/include -I/home/pan/petsc-3.24.5/arch-linux-c-release/include
# 可选调试标志（取消注释以启用）
#CXXFLAGS = -O3 -march=core-avx2 -g -Wall    # 调试和警告
#CXXFLAGS = -O3 -fopenmp    # 启用 OpenMP 并行
#CXXFLAGS = -O3 -mavx       # 启用 AVX 指令
CXXLIB = -Wl,-rpath,/home/pan/petsc-3.24.5/arch-linux-c-release/lib -L/home/pan/petsc-3.24.5/arch-linux-c-release/lib -lpetsc -lHYPRE -lf2clapack -lf2cblas -lm -lstdc++ -lquadmath

# 源文件目录
DATA_DIR = .
FSI_DIR = $(DATA_DIR)/fsi
FLUID_DIR = $(DATA_DIR)/fluid
SOLID_DIR = $(DATA_DIR)/solid
MODULE_DIR = ../module

# 自动检测 data 目录下的 .cpp 文件（仅包含根目录和指定子目录）
# 当前默认使用 fluid，用户可手动改为 solid
DATA_SUBDIR = $(FSI_DIR)
# 查找 data 根目录下的 .cpp 文件
DATA_ROOT_SRC = $(wildcard $(DATA_DIR)/*.cpp)
# 查找 DATA_SUBDIR 下的 .cpp 文件
DATA_SUB_SRC = $(shell find $(DATA_SUBDIR) -type f -name '*.cpp')

# 合并 data 根目录和指定子目录的源文件
DATA_SRC = $(DATA_ROOT_SRC) $(DATA_SUB_SRC)

# 自动发现 module 目录下的所有源文件（排除未更新的 MPM fluid solver 和 explicit solid）
MODULE_SRC = $(shell find $(MODULE_DIR) -type f -name '*.cpp' -not -path '$(MODULE_DIR)/fluid/MPM/*' -not -path '$(MODULE_DIR)/solid/explicit/*')

# 所有源文件
SRC = $(DATA_SRC) $(MODULE_SRC)

# 生成目标文件
OBJS = $(SRC:%.cpp=%.o)

# 默认目标
all: $(PROG)

# 编译规则：将 .cpp 编译为 .o 文件
.cpp.o:
	$(CXX) $(CXXFLAGS) $(CXXLIB) -c $< -o $@

# 链接生成可执行文件
$(PROG): $(OBJS)
	$(CXX) $(OBJS) $(CXXFLAGS) $(CXXLIB) -o $(PROG)

# 清理生成文件
clean:
	@rm -f $(OBJS) $(PROG) *~

# 声明伪目标
.PHONY: all clean
