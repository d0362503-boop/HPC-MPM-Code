# 编译器设置
#CXX = mpic++   # --- for MPICH/OpenMPI ---
CXX = mpicxx   # --- for MPICH/OpenMPI ---
#CXX = mpiicpx  # --- for Intel OneAPI ---

# PETSc 设置
PETSC_DIR  = /home/pan/petsc-3.24.5
PETSC_ARCH = arch-linux-c-release
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables

CXXFLAGS = -O3 -march=native -std=c++17 $(PETSC_CC_INCLUDES) #-g -Wall -Wextra -Wshadow #-fsanitize=address,undefined 

# 目录设置
BUILD_DIR = build
MODULE_DIR = module
MODULE_SOLVER_DIR = $(MODULE_DIR)/solver
MODULE_FLUID_DIR = $(MODULE_DIR)/fluid
MODULE_SOLID_DIR = $(MODULE_DIR)/solid

# 手动切换 src 目录（src_fluid, src_solid, src_fsi）
USE_SRC_FLUID = false
USE_SRC_SOLID = false #true
USE_SRC_FSI = true

# 手动切换是否编译 module 下的子目录（true 或 false）
USE_SOLVER = true
USE_FLUID = true
USE_SOLID = true

# 手动切换 fluid 中的 FEM 和 MPM 模式（true 或 false）
USE_FEM = true
USE_MPM = false #true

# 手动切换 solid 中的 explicit 和 implicit 模式（true 或 false）
USE_EXPLICIT_SOLID = false #true
USE_IMPLICIT_SOLID = true

# 自动处理 implicit 优先逻辑
SOLID_MODE = $(if $(filter true,$(USE_IMPLICIT_SOLID)),implicit,$(if $(filter true,$(USE_EXPLICIT_SOLID)),explicit,none))
MODULE_SOLID_EXPLICIT_DIR = $(MODULE_SOLID_DIR)/explicit
MODULE_SOLID_IMPLICIT_DIR = $(MODULE_SOLID_DIR)/implicit
SRC_SOLID_EXPLICIT_DIR = src_solid/explicit
SRC_SOLID_IMPLICIT_DIR = src_solid/implicit

# 流体子目录设置
MODULE_FLUID_FEM_DIR = $(MODULE_FLUID_DIR)/FEM
MODULE_FLUID_MPM_DIR = $(MODULE_FLUID_DIR)/MPM
SRC_FLUID_FEM_DIR = src_fluid/FEM
SRC_FLUID_MPM_DIR = src_fluid/MPM

# 主程序文件
MAIN_SRC = MPM_main.cpp
MAIN_OBJ = $(BUILD_DIR)/$(MAIN_SRC:.cpp=.o)

# 自动查找源文件
# module 目录下的文件（不包含 solver, fluid, solid 子目录）
MODULE_BASE_SRCS = $(shell find $(MODULE_DIR) -type f -name '*.cpp' -not -path '$(MODULE_SOLVER_DIR)/*' -not -path '$(MODULE_FLUID_DIR)/*' -not -path '$(MODULE_SOLID_DIR)/*')
# solver 子目录下的文件（动态决定是否包含）
MODULE_SOLVER_SRCS = $(if $(filter true,$(USE_SOLVER)),$(shell find $(MODULE_SOLVER_DIR) -type f -name '*.cpp'))
# fluid 子目录下的文件（动态决定是否包含 FEM 和 MPM）
MODULE_FLUID_FEM_SRCS = $(if $(filter true,$(USE_FLUID)),$(if $(filter true,$(USE_FEM)),$(shell find $(MODULE_FLUID_FEM_DIR) -type f -name '*.cpp')))
MODULE_FLUID_MPM_SRCS = $(if $(filter true,$(USE_FLUID)),$(if $(filter true,$(USE_MPM)),$(shell find $(MODULE_FLUID_MPM_DIR) -type f -name '*.cpp')))
MODULE_FLUID_SRCS = $(MODULE_FLUID_FEM_SRCS) $(MODULE_FLUID_MPM_SRCS)
# solid 子目录下的文件（动态决定是否包含，并根据模式选择）
MODULE_SOLID_SUBDIR_SRCS = $(if $(filter true,$(USE_SOLID)),$(shell find $(MODULE_SOLID_$(shell echo $(SOLID_MODE) | tr '[:lower:]' '[:upper:]')_DIR) -type f -name '*.cpp'))
# solid 目录下直接放置的文件（动态决定是否包含）
MODULE_SOLID_BASE_SRCS = $(if $(filter true,$(USE_SOLID)),$(shell find $(MODULE_SOLID_DIR) -maxdepth 1 -type f -name '*.cpp'))
# 合并 solid 的源文件
MODULE_SOLID_SRCS = $(MODULE_SOLID_BASE_SRCS) $(MODULE_SOLID_SUBDIR_SRCS)
# src 目录下的文件（动态决定是否包含 FEM 和 MPM）
SRC_FLUID_FEM_SRCS = $(if $(filter true,$(USE_SRC_FLUID)),$(if $(filter true,$(USE_FEM)),$(shell find $(SRC_FLUID_FEM_DIR) -type f -name '*.cpp')))
SRC_FLUID_MPM_SRCS = $(if $(filter true,$(USE_SRC_FLUID)),$(if $(filter true,$(USE_MPM)),$(shell find $(SRC_FLUID_MPM_DIR) -type f -name '*.cpp')))
SRC_FLUID_SRCS = $(SRC_FLUID_FEM_SRCS) $(SRC_FLUID_MPM_SRCS)
SRC_SOLID_SRCS = $(if $(filter true,$(USE_SRC_SOLID)),$(shell find $(SRC_SOLID_$(shell echo $(SOLID_MODE) | tr '[:lower:]' '[:upper:]')_DIR) -type f -name '*.cpp'))
SRC_FSI_SRCS = $(if $(filter true,$(USE_SRC_FSI)),$(shell find src_fsi -type f -name '*.cpp'))
# 合并 module 和 src 的源文件
MODULE_SRCS = $(MODULE_BASE_SRCS) $(MODULE_SOLVER_SRCS) $(MODULE_FLUID_SRCS) $(MODULE_SOLID_SRCS)
SRC_SRCS = $(SRC_FLUID_SRCS) $(SRC_SOLID_SRCS) $(SRC_FSI_SRCS)

# 所有源文件
SRCS = $(MODULE_SRCS) $(SRC_SRCS) $(MAIN_SRC)

# 生成目标对象路径（保持原始目录结构）
OBJS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# 默认目标
all: $(BUILD_DIR)/MPM.exe

# 创建 build 下的中间目录
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 链接最终可执行文件
$(BUILD_DIR)/MPM.exe: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(PETSC_WITH_EXTERNAL_LIB)

# 清理
clean:
	@find $(BUILD_DIR) -type f -name '*.o' -delete
	@find $(BUILD_DIR) -type d -empty -not -path $(BUILD_DIR) -delete

.PHONY: all clean

parallel:
	@echo "Use 'make -jN' to compile with N parallel jobs, e.g., 'make -j8'"
