# CMake 重构清单

> 面向当前 `MPM-Code` 的可执行重构 checklist。目标不是“一次性推翻重写”，而是在不破坏现有求解器工作流的前提下，逐步把构建系统整理成可维护、可移植、可诊断的状态。

## 1. 当前问题概览

基于目前仓库中的以下入口文件整理：

- `CMakeLists.txt`
- `work/CMakeLists.txt`
- `module/CMakeLists.txt`
- `cmake/BuildHDF5.cmake`
- `data/CMakeLists.txt` 及其子目录

当前 CMake 配置能工作，但存在几个明显的结构性问题：

- 构建模式选择依赖手工注释/取消注释 `add_subdirectory(...)`，不利于批量测试和自动化。
- 第三方依赖策略不统一：MPI/ZLIB 走 `find_package`，HDF5/PETSc 部分是手工脚本式管理。
- `OBJECT` library 用得比较多，但 target 之间的职责边界还不够清晰。
- 编译器兼容逻辑此前分散，虽然已经开始收口，但还没有完全形成统一约束层。
- 根目录对构建目录名、优化选项、外部依赖路径做了较强假设，不利于集群和多环境切换。
- 缺少标准化的“配置验证”和“失败诊断”路径，出错时经常只能靠人工排查。

## 2. 重构目标

建议把这次 CMake 重构的终态目标定为：

- 本地和集群都能用一致命令完成配置与编译。
- Solver 模式选择改为显式 CMake option，而不是改源码注释。
- 外部依赖的来源、优先级、失败提示都清晰可控。
- 所有 target 的 include、link、compile definitions 都通过 target 传播，而不是靠全局变量或隐式副作用。
- 新增源文件时，维护者能快速判断要改哪里，不容易漏改。
- 后续如果继续拆分 `src_solid` / `src_fluid` / `src_fsi`，构建系统不会成为阻力。

## 3. 重构原则

- 优先“整理结构”，而不是先追求功能扩展。
- 每一阶段都应可独立提交，并保持可编译。
- 尽量保留现有目录组织，不做无收益的大搬家。
- 先把“模式选择、依赖策略、target 边界”理顺，再考虑安装、导出、打包。

## 4. 分阶段 Checklist

## Phase 0: 建立基线

- [ ] 在 `docs/` 中补一份构建矩阵说明，至少覆盖：
  - Ubuntu/WSL 本地 GCC 11+
  - 集群 GCC 8.1
  - MPI on/off
  - HDF5 on/off
  - PETSc 外部已有安装 / 内部构建
- [ ] 固化当前推荐命令，统一写成两套：
  - 现代写法：`cmake -S . -B build`
  - 兼容老环境写法：`mkdir build && cd build && cmake ..`
- [ ] 明确当前“成功配置”的最小判据：
  - 能找到 MPI
  - 能找到或构建 PETSc
  - 能找到或构建 HDF5
  - 能成功生成 `build/MPM`
- [ ] 记录现有 target 拓扑：
  - `MPM`
  - `mpm_modules`
  - `mpm_src_solid`
  - `mpm_src_fluid`
  - `mpm_src_fsi`
  - `makinput_*`
  - `makdivide_*`

验证：

- [ ] 在本地执行一次全新配置和编译。
- [ ] 在集群执行一次全新配置和编译。
- [ ] 保存成功日志，后续作为回归对照。

## Phase 1: 根 CMake 入口收口

重点文件：

- `CMakeLists.txt`
- `cmake/options.cmake`

Checklist：

- [x] 把根目录里“项目声明、编译标准、全局策略、选项入口、子目录入口、主程序链接”这几类职责明确分段。
- [ ] 把当前对 `build` 目录名的硬性限制改成可选检查，而不是强制失败。
- [x] 把 `-march=native` 改成 option 控制，例如 `MPM_ENABLE_NATIVE_ARCH`。
- [x] 明确区分：
  - 用户配置项
  - 自动探测结果
  - 内部兼容性辅助 target
- [ ] 检查所有 `set(... CACHE ... FORCE)` 的使用场景，避免不必要地覆盖用户传入值。
- [x] 为常见失败场景补更清晰的 `message(FATAL_ERROR ...)` 提示。（配置摘要已加入）

建议结果：

- 根 `CMakeLists.txt` 只负责 orchestration，不承载具体依赖构建细节。

验证：

- [ ] `cmake ..` 和 `cmake -S . -B build` 都能工作。
- [ ] 切换 `Release/Debug` 不需要手工清缓存。
- [ ] 关闭 `native arch` 后，集群不会因为 CPU 差异出现潜在兼容问题。

## Phase 2: Solver 模式选择去手工注释化

重点文件：

- `work/CMakeLists.txt`
- `MPM_main.cpp`

Checklist：

- [ ] 用明确 option 替代手工注释，例如：
  - `MPM_BUILD_SOLID`
  - `MPM_BUILD_FLUID`
  - `MPM_BUILD_FSI`
- [ ] 保证至少有一个模式被启用，否则配置阶段直接报错。
- [ ] 如果多个模式同时启用，定义清楚行为：
  - 只构建多个库但主程序选一个
  - 或者明确禁止多选
- [ ] 把 `MPM_main.cpp` 中与模式选择强耦合的部分同步整理，避免 “CMake 选的是 A，主程序调用的是 B”。
- [ ] 若有必要，引入单独的配置头文件或 compile definition，用于显式传递当前构建模式。

建议结果：

- 切换构建模式不再需要修改源码注释。

验证：

- [ ] `SOLID` 单独配置通过。
- [ ] `FLUID` 单独配置通过。
- [ ] `FSI` 单独配置通过。
- [ ] 错误组合能在配置阶段给出明确报错。

## Phase 3: target 依赖传播整理

重点文件：

- `module/CMakeLists.txt`
- `module/solver/*`
- `module/fluid/*`
- `module/solid/*`
- `work/src_*/*`
- `data/*/CMakeLists.txt`

Checklist：

- [x] 重新审视 `mpm_modules` 是否继续使用 `OBJECT` library：
  - 若继续使用，明确它是“对象聚合层”
  - 若后续条件成熟，可评估改为 `STATIC` library
- [x] 统一所有公共编译兼容逻辑挂到专门的 interface target 上。
- [x] 检查所有 target 的 `PUBLIC/PRIVATE/INTERFACE` 使用是否准确。
- [x] 清理对子目录中重复 link 相同系统库的写法，避免多点散落。
- [x] 明确 `data/generate` 与 `data/divide` 中公共对象层和最终可执行层的职责。
- [x] 为 `makinput_*` / `makdivide_*` / `MPM` 建立一致的链接策略，避免“某个工具补了兼容库、另一个忘补”的情况再次出现。

建议结果：

- 谁依赖谁，由 `target_link_libraries()` 明确表达，而不是靠“碰巧在别处已经 link 过”。

验证：

- [ ] GCC 8.1 下所有使用 `std::filesystem` 的 target 都能链接通过。
- [ ] 新增一个使用 `std::filesystem` 的可执行 target 时，不需要再单独手补 `stdc++fs`。

## Phase 4: 第三方依赖策略统一

重点文件：

- `cmake/BuildHDF5.cmake`
- `cmake/BuildPETSc.cmake`
- `cmake/options.cmake`
- 根 `CMakeLists.txt`

Checklist：

- [ ] 明确依赖优先级策略：
  - 优先系统安装
  - 找不到再使用仓库内 tarball / external 构建
  - 或完全反过来
- [ ] 为 HDF5 增加“tarball 完整性检查”或至少更清晰的失败提示。
- [ ] 针对 Git LFS pointer 文件场景，补专门报错提示。
- [ ] 检查 `BuildHDF5.cmake` 中 `execute_process` 的错误输出是否足够保真。
- [ ] 评估是否用 `ExternalProject_Add` 替代手写解压-配置-编译流程。
- [ ] 明确 `PETSC_DIR` / `PETSC_ARCH` 的优先级与回退规则，并写进注释和文档。
- [x] HDF5、PETSc、MPI、ZLIB 的发现结果在配置阶段统一打印摘要。

建议结果：

- 第三方依赖失败时，能快速分辨是：
  - 路径错误
  - 版本错误
  - tarball 损坏
  - Git LFS pointer
  - 编译器不兼容

验证：

- [ ] 正常 tarball 能自动构建 HDF5。
- [ ] pointer 文件会给出可读报错，而不是只看到 `not in gzip format`。
- [ ] 集群已有 PETSc 时，不会误走错误路径。

## Phase 5: 配置选项与模式文档化

重点文件：

- `cmake/options.cmake`
- `docs/`

Checklist：

- [ ] 给每个 option 补一句话职责说明。
- [ ] 区分“用户应该手动设置的变量”和“内部中间变量”。
- [ ] 为常用场景给出推荐命令模板：
  - 本地开发
  - 集群编译
  - 只编译数据工具
  - 调试 HDF5/PETSc 发现问题
- [ ] 补一份 “最小构建手册”，集中说明：
  - 需要哪些依赖
  - 哪些目录需要 Git LFS 完整文件
  - 哪些环境变量会影响配置

建议结果：

- 新人不看源码也能知道该传哪些 `-D` 选项。

验证：

- [ ] 根据文档在空环境里复现实验，确认文档不是过期的。

## Phase 6: data 子工程整理

重点文件：

- `data/CMakeLists.txt`
- `data/generate/CMakeLists.txt`
- `data/divide/CMakeLists.txt`
- `data/generate/*/CMakeLists.txt`
- `data/divide/*/CMakeLists.txt`

Checklist：

- [ ] 统一 `generate` 和 `divide` 两条工具链的命名风格。
- [ ] 统一公共 object target 的命名和职责边界。
- [ ] 检查是否存在“只有某个子工具手工 link 某库”的不对称写法。
- [ ] 明确工具程序是否应该继承 `mpm_modules` 的全部依赖，还是只依赖局部子集。
- [ ] 若工具程序彼此结构重复，抽一层小的 CMake helper function，减少复制粘贴。

建议结果：

- 数据工具构建逻辑与主求解器构建逻辑风格一致。

验证：

- [ ] `makinput_*` 全部可编译。
- [ ] `makdivide_*` 全部可编译。
- [ ] GCC 8.1 环境下工具链不再遗漏 `filesystem` 兼容链接。

## Phase 7: 测试与诊断入口补齐

重点文件：

- 根 `CMakeLists.txt`
- `docs/`

Checklist：

- [ ] 至少补一个 `cmake --build build --target help` 可读的目标组织。
- [ ] 如果暂时不引入单元测试，也至少补“配置 smoke test”的文档和命令。
- [x] 可考虑增加简单自检 target，例如打印关键依赖路径和选项摘要。
- [ ] 统一常见错误的排查顺序：
  - 编译器版本
  - MPI
  - PETSc
  - HDF5 tarball
  - Git LFS
  - CMake 缓存污染

建议结果：

- 出问题时先看配置摘要和文档，而不是直接陷入源码排查。

验证：

- [ ] 让一个不熟悉当前构建系统的人按文档排查一次，确认路径清晰。

## Phase 8: Legacy 构建路径清理（已完成）

决策：根 `Makefile` 及 `data/` 下的旧 Makefile 已彻底废弃，不再维护。

重点文件：

- 根 `Makefile`（保留文件本体作为历史参考，但文档中不再描述其用法）
- `data/` 下旧 Makefile
- `docs/`、`README.md`、`AGENTS.md`、各 `README.md`

Checklist：

- [x] 明确 legacy Makefile 不再支持。
- [x] 从所有 Markdown 文档中移除 Makefile 用法说明。
- [x] 文档统一以 CMake 为唯一官方构建路径。
- [x] 新增 `.cpp` 文件只需更新 `CMakeLists.txt`。

建议结果：

- 团队成员知道 CMake 是唯一官方主线，不存在“双主线”误导。

验证：

- [x] 文档和实际行为一致；除本清理记录外，所有使用说明性质的 `*.md` 均已移除 `make -j`、`make clean`、`Makefile` 等旧构建说明。

## 5. 优先级建议

如果只分三轮推进，我建议顺序如下：

### 第一轮：先把最容易踩坑的结构问题收口

- [ ] Phase 1
- [ ] Phase 2
- [ ] Phase 3

### 第二轮：把外部依赖问题压下去

- [ ] Phase 4
- [ ] Phase 5
- [ ] Phase 6

### 第三轮：补长期维护能力

- [ ] Phase 7
- [ ] Phase 8

## 6. 建议的完成判据

满足下面这些条件时，可以认为这轮 CMake 重构基本达标：

- [ ] 不需要修改源码注释就能切换 `SOLID` / `FLUID` / `FSI` 构建。
- [ ] GCC 8.1 和 GCC 11+ 都能稳定编译通过。
- [ ] `std::filesystem` 兼容逻辑只维护一处。
- [ ] HDF5 tarball 错误能给出明确诊断。
- [ ] PETSc/HDF5/MPI 的发现路径和优先级在文档中说清楚。
- [ ] `MPM`、`makinput_*`、`makdivide_*` 的链接行为一致且可解释。
- [ ] 新增一个可执行 target 时，维护者能很快知道应该 link 哪些公共 target。

## 7. 我对当前框架的建议结论

这套 CMake 不是“不能用”，但确实还处于“能编起来、但维护成本偏高”的阶段。最值得先动手的不是大改目录，而是把以下三件事做实：

- Solver 模式选择 option 化
- 第三方依赖策略统一化
- target 依赖传播规范化

如果这三件事先落地，后面无论你继续扩展 `SOLID MPM`、`FSI` 还是数据工具，构建层都会顺很多。
