# 后缀规则
.SUFFIXES :
.SUFFIXES : .o .cpp

# 编译器选择
#CXX = mpic++    # 使用 OpenMPI 进行 MPI 编译
CXX = mpicxx    # 可选 MPICH
#CXX = mpiicpx   # 可选 Intel OneAPI
#CXX = icpx      # 可选 Intel C++ 编译器
#CXX = g++       # 可选 GNU C++ 编译器

# 输出目标
PROG = ./divide.exe

# 编译标志
CXXFLAGS = -O3 -std=c++17 -I/home/pan/petsc-3.24.5/include -I/home/pan/petsc-3.24.5/arch-linux-c-release/include
#CXXLIB =                   # 如果需要添加库（例如 -lm）

# 源文件目录
DIVIDE_DIR = .
DIVIDE_FLUID_DIR = $(DIVIDE_DIR)/divide_fluid
DIVIDE_SOLID_DIR = $(DIVIDE_DIR)/divide_solid
DIVIDE_FSI_DIR = $(DIVIDE_DIR)/divide_fsi
MODULE_DIR = ../../module

# 当前默认使用 divide_fluid，用户可手动改
DIVIDE_SRC_DIR = $(DIVIDE_FSI_DIR)
DIVIDE_SRC = $(shell find $(DIVIDE_SRC_DIR) -type f -name '*.cpp') \
             $(wildcard $(DIVIDE_DIR)/*.cpp)

# 手动指定 module 目录下的源文件
MODULE_SRC = \
    $(MODULE_DIR)/mesh.cpp \
    $(MODULE_DIR)/bc.cpp \
    $(MODULE_DIR)/data_io.cpp \
    $(MODULE_DIR)/material_point.cpp \
    $(MODULE_DIR)/solid/solid_material_point.cpp \
    $(MODULE_DIR)/solid/solid_data_io.cpp \
    $(MODULE_DIR)/fluid/FEM/fem_data_io.cpp \
    $(MODULE_DIR)/shape_function.cpp

# 所有源文件
SRC = $(DIVIDE_SRC) $(MODULE_SRC)

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
