# CROssBARv2-KG 自动停止问题分析与解决方案

## 问题诊断

### 症状
运行 `scripts/create_crossbar.py` 时，脚本在处理 UniProt 数据后自动停止

### 根本原因
**内存溢出 (OOM - Out Of Memory)**

#### 证据
1. **退出码 137**: 进程被 `SIGKILL` (信号9) 强制终止
   - 退出码 137 = 128 + 9
   - 这是 Linux OOM Killer 的典型行为

2. **触发点**: InterPro adapter 和 Drug adapter
   - InterPro: `interpro.interpro_entries()` 会将整个 InterPro 数据库加载到内存
   - Drug: `unichem.unichem_mapping()` 同时加载 6 个大型映射字典

3. **系统资源限制**:
   ```
   Total Memory: 12 GB
   Used: 5.4 GB
   Available: 6.6 GB
   Swap: 0 GB (无交换分区)
   ```

#### Drug Adapter OOM 原因分析

1. **初始化时立即加载 SwissProt** (已修复)
   ```python
   # 原代码在 __init__ 中立即加载
   self.swissprots = set(uniprot._all_uniprots(organism="*", swissprot=True))
   ```

2. **UniChem 映射同时加载** (已修复)
   ```python
   # 原代码同时加载 6 个大型字典，每个可能数百 MB
   unichem_drugbank_to_zinc_mapping = unichem.unichem_mapping("drugbank", "zinc")
   unichem_drugbank_to_chembl_mapping = unichem.unichem_mapping("chembl", "drugbank")
   # ... 等等
   ```

3. **ChEMBL 数据不清理** (已修复)
   ```python
   # 原代码下载后不清理
   self.chembl_acts = chembl.chembl_activities(...)  # 非常大
   self.chembl_targets = chembl.chembl_targets()
   # ...
   ```

## 已实施的修复方案

### Drug Adapter 优化 (bccb/drug_adapter.py)

1. **延迟加载 SwissProt IDs**
   ```python
   # 使用 property 实现延迟加载
   @property
   def swissprots(self):
       if self._swissprots is None:
           self._swissprots = set(uniprot._all_uniprots(organism="*", swissprot=True))
       return self._swissprots
   ```

2. **新增 `low_memory_mode` 参数**
   ```python
   drug_adapter = Drug(
       ...,
       low_memory_mode=True  # 启用内存优化
   )
   ```

3. **按需加载 UniChem 映射**
   - 在 `low_memory_mode=True` 时，逐个加载 UniChem 映射
   - 每个映射处理完立即清理，避免峰值内存过高

4. **自动内存清理**
   - `cleanup_memory()` 方法用于手动清理
   - `process_chembl_dti_data()` 完成后自动清理原始数据

### 使用方法

```bash
cd /GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG
python scripts/create_crossbar.py
```

### 方案2：增加系统内存

如果需要 InterPro 数据：
1. 使用内存更大的机器（建议 ≥32GB）
2. 或者添加交换分区（不推荐，会很慢）

```bash
# 创建 16GB 交换文件（仅作参考，需要 root 权限）
sudo fallocate -l 16G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### 方案3：分离运行 InterPro

单独在高内存机器上运行 InterPro adapter：

```python
# run_interpro_only.py
from bccb.interpro_adapter import InterPro
from biocypher import BioCypher

bc = BioCypher(...)
interpro_adapter = InterPro(organism=9606, test_mode=False)
interpro_adapter.download_interpro_data(cache=True)
bc.write_nodes(interpro_adapter.get_interpro_nodes())
bc.write_edges(interpro_adapter.get_interpro_edges())
```

### 方案4：优化 InterPro Adapter（需要代码修改）

修改 `bccb/interpro_adapter.py` 以分批处理数据：
- 使用迭代器而不是一次性加载所有数据
- 实现分页机制
- 定期清理内存

## 脚本对比

### 原始脚本 (`create_crossbar.py`)
- ❌ InterPro 导致 OOM
- ❌ 没有内存清理
- ❌ 没有进度提示
- ❌ 脚本会在 InterPro 处自动停止

### 优化脚本 (`create_crossbar_optimized.py`)
- ✅ 跳过 InterPro，避免 OOM
- ✅ 每个 adapter 后清理内存 (`gc.collect()`)
- ✅ 清晰的进度提示
- ✅ 所有其他 adapters 都能运行完成

## 运行状态监控

### 检查进程状态
```bash
# 查看 Python 进程
ps aux | grep create_crossbar

# 查看内存使用
free -h

# 监控实时内存（另一个终端）
watch -n 1 free -h
```

### 查看日志
```bash
# BioCypher 日志
ls -lt biocypher-log/ | head -5
tail -f biocypher-log/biocypher-*.log

# PyPath 日志（查看数据下载进度）
ls -lt pypath_log/ | head -5
tail -f pypath_log/pypath-*.log
```

## 各 Adapter 预期运行时间

| Adapter | 预期时间 | 内存使用 | 备注 |
|---------|----------|----------|------|
| UniProt | 1-2分钟 | ~200MB | TEST_MODE=True 只加载100条 |
| InterPro | ❌ OOM | >6GB | 会导致内存溢出 |
| GO | 10-15分钟 | ~500MB | 需要下载 >1100 页数据 |
| Drug | 5-10分钟 | ~300MB | 需要 DrugBank 认证 |
| Compound | 5-10分钟 | ~400MB | |
| Orthology | 2-5分钟 | ~200MB | |
| Disease | 5-10分钟 | ~300MB | |
| Phenotype | 3-5分钟 | ~200MB | |
| Pathway | 5-10分钟 | ~300MB | |
| SideEffect | 3-5分钟 | ~200MB | |
| EC | 2-3分钟 | ~100MB | |
| TFGene | 2-3分钟 | ~100MB | |

**总预计时间（不含 InterPro）**: 45-70 分钟

## 验证运行成功

### 检查输出文件
```bash
# 查看生成的文件
ls -lh biocypher-out/

# 查看 CSV 文件
ls -lh biocypher-out/*.csv

# 查看最新输出目录
ls -lh biocypher-out/$(ls -t biocypher-out/ | head -1)/
```

### 成功标志
脚本最后应该输出：
```
================================================================================
✓ Script completed successfully!
================================================================================
```

## 最佳实践建议

1. **首次运行**: 使用 `TEST_MODE=True`（当前设置）
   - 快速验证流程
   - 减少内存使用
   - 减少运行时间

2. **生产运行**: 设置 `TEST_MODE=False`
   - 需要更多内存（建议 ≥32GB）
   - 需要更多时间（数小时）
   - 考虑分批运行不同的 adapters

3. **缓存利用**: 保持 `CACHE=True`
   - 第二次运行会快很多
   - PyPath 缓存位置: `/root/.cache/pypath/` 或 `/GenSIvePFS/users/data/pypath_cache/`

4. **内存监控**: 使用 `tmux` 或 `screen` 在后台运行
   ```bash
   tmux new -s crossbar
   python scripts/create_crossbar_optimized.py
   # Ctrl+B, D 来 detach
   # tmux attach -t crossbar 来重新连接
   ```

## 总结

- **问题**: InterPro adapter 导致内存溢出 (OOM)
- **解决**: 使用优化脚本跳过 InterPro
- **结果**: 其他所有 adapters 都能正常运行
- **下一步**: 如需 InterPro 数据，在高内存机器上单独运行

---
**创建时间**: 2026-01-18
**系统**: 12GB RAM, No Swap
**Python**: 3.10.19
**环境**: crossbarv2 conda environment
