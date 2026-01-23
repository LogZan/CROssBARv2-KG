# PPI/STRING 性能瓶颈分析 - 使用指南

## 快速开始

如果你想了解为什么 STRING 数据处理需要 50+ 小时，请阅读:

📄 **[STRING_BOTTLENECK_SUMMARY.md](./STRING_BOTTLENECK_SUMMARY.md)** (5分钟阅读)

这个文档回答了:
- ✓ 为什么有缓存还需要 15秒/物种？
- ✓ 12,535 个物种需要多久？
- ✓ 主要时间花在哪里？
- ✓ 这是正常的吗？

## 详细分析

如果你想深入了解技术细节:

📄 **[STRING_BOTTLENECK_DEEP_ANALYSIS.md](./STRING_BOTTLENECK_DEEP_ANALYSIS.md)** (15分钟阅读)

包含:
- PyPath 源码分析
- 性能测试数据
- 时间分解细节
- 缓存机制说明

## 其他报告

最初的分析 (后来发现理解有误):

- **[PPI_BOTTLENECK_ANALYSIS.md](./PPI_BOTTLENECK_ANALYSIS.md)** - 初步分析 (假设是配置错误)
- **[PPI_BOTTLENECK_SUMMARY_CN.md](./PPI_BOTTLENECK_SUMMARY_CN.md)** - 中文总结 (配置角度)

## 运行测试

所有测试脚本位于 `test/` 目录:

### 1. 快速性能测试 (1分钟)
```bash
python test/profile_ppi_simple.py
```
测试:
- normalize_curie 性能
- 缓存优化效果
- 不同实现方法对比

### 2. STRING 瓶颈分析 (3-5分钟)
```bash
python test/analyze_string_bottleneck.py
```
测试:
- Physical interactions 读取时间
- Links 读取时间
- 完整流程耗时

### 3. 配置和时间估算 (几秒)
```bash
python test/analyze_real_bottleneck.py
```
分析:
- 当前配置
- 物种数量
- 总时间估算

### 4. 精确复现 (可能失败)
```bash
python test/profile_string_exact.py
```
精确复现 ppi_adapter.py 的调用流程

### 5. 完整 PPI 分析 (需要较长时间)
```bash
python test/profile_ppi_bottleneck.py
```
完整的 PPI adapter 性能测试 (使用 test_mode)

## 关键发现总结

### 问题
为什么处理 12,535 个物种的 STRING 数据需要 55+ 小时？

### 答案
**这是正常的！**

每个物种需要约 15 秒:
- 解压 .gz 文件 (33%)
- 文本解析 (47%)
- 构建字典 (7%)
- 过滤查询 (13%)

12,535 × 15秒 = 52.2 小时

### 原因
1. **缓存避免了网络下载，但无法避免**:
   - 解压缩 (CPU)
   - 文本解析 (CPU)
   - 数据结构构建 (CPU)

2. **PyPath 设计问题**:
   - 为了获取 physical_combined_score
   - 必须读取 ALL physical interactions
   - 即使用户只要 high_confidence

3. **数据量巨大**:
   - 人类: 1.4M physical + 11M links
   - 每个物种都要独立处理

## 文件清单

### 测试脚本 (test/)
- ✓ profile_ppi_simple.py (7.6K) - PPI边生成性能测试
- ✓ analyze_string_bottleneck.py (7.6K) - STRING读取深度分析
- ✓ profile_string_exact.py (6.9K) - 精确复现实际调用
- ✓ analyze_real_bottleneck.py (5.7K) - 配置和时间分析
- ✓ profile_ppi_bottleneck.py (7.7K) - 完整PPI性能分析

### 分析报告 (docs/)
- ✓ STRING_BOTTLENECK_SUMMARY.md (4.1K) - **推荐阅读** 🌟
- ✓ STRING_BOTTLENECK_DEEP_ANALYSIS.md (9.4K) - 详细技术分析
- ✓ PPI_BOTTLENECK_ANALYSIS.md (6.1K) - 初步分析
- ✓ PPI_BOTTLENECK_SUMMARY_CN.md (6.9K) - 中文总结

## 结论

对于处理所有物种的 STRING 数据:
- ✓ 52+ 小时是正常的处理时间
- ✓ 不是 bug 或配置错误
- ✓ 这是数据处理的固有成本

如果想加速:
- 考虑并行处理
- 优化 PyPath 源码
- 预处理为二进制格式

---

📅 分析日期: 2026-01-23  
🔬 分析方法: 源码分析 + 性能测试 + 数据验证  
✅ 状态: 已完成
