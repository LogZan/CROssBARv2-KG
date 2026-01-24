# PPI Adapter 并行化实现方案讨论

## 📋 当前代码分析

### 关键方法：`download_string_data()`

**位置**: `bccb/ppi_adapter.py` line 651-759

**当前流程**:
```python
def download_string_data(self, cache=False):
    # 1. 获取物种列表 (line 661-668)
    if self.organism is None:
        string_species = string.string_species()
        self.tax_ids = list(string_species.keys())  # 12,535 物种
    else:
        self.tax_ids = [self.organism]
    
    # 2. 构建 string_to_uniprot 映射 (line 674-684)
    uniprot_to_string = uniprot.uniprot_data(...)
    self.string_to_uniprot = collections.defaultdict(list)
    for k, v in uniprot_to_string.items():
        # 构建映射...
    
    # 3. 串行处理每个物种 (line 706-751) ← 这里需要并行化
    for tax in tqdm(self.tax_ids, desc="Retrieving STRING data"):
        organism_string_ints = [...]
        self.string_ints.extend(organism_string_ints)
```

**当前问题**:
- ✓ 步骤1和2很快，不需要优化
- ✗ 步骤3是瓶颈：12,535次循环，每次~15秒

## 🎯 并行化目标

**只并行化步骤3**: 物种处理循环

**不改变**:
- 公共接口（方法签名）
- 数据结构（self.string_ints）
- 错误处理逻辑
- 向后兼容性

## 💡 实施方案选项

### 方案A：添加新参数 `n_workers` (推荐) ⭐⭐⭐⭐⭐

**思路**: 向 `download_string_data()` 添加可选参数

```python
def download_string_data(self, cache=False, n_workers=None):
    """
    Args:
        cache: 使用缓存
        n_workers: 并发进程数
            - None: 自动检测（推荐）
            - 1: 串行处理（向后兼容）
            - >1: 并行处理
    """
```

**优点**:
- ✓ 向后兼容（默认行为可配置）
- ✓ 灵活性高（用户可控制）
- ✓ 易于测试（串行/并行都能测）
- ✓ 不影响现有代码

**缺点**:
- ✗ 需要修改函数签名
- ✗ 需要文档更新

**是否破坏兼容性**: 
- 不会，因为是可选参数
- 默认值可以是 `None`（自动）或 `1`（串行）

### 方案B：添加新方法 `download_string_data_parallel()`

**思路**: 保留原方法，添加并行版本

```python
# 保持原有方法不变
def download_string_data(self, cache=False):
    # 原有串行实现
    pass

# 新增并行方法
def download_string_data_parallel(self, cache=False, n_workers=None):
    # 并行实现
    pass
```

**优点**:
- ✓ 完全不影响现有代码
- ✓ 容易对比测试
- ✓ 渐进式迁移

**缺点**:
- ✗ 代码重复
- ✗ 维护两个版本
- ✗ 用户需要知道选哪个

### 方案C：自动检测，智能切换

**思路**: 根据物种数量自动选择串行/并行

```python
def download_string_data(self, cache=False, parallel_threshold=100):
    if len(self.tax_ids) >= parallel_threshold:
        self._download_string_parallel(...)
    else:
        self._download_string_serial(...)
```

**优点**:
- ✓ 用户无感知
- ✓ 自动优化

**缺点**:
- ✗ 用户失去控制
- ✗ 调试困难
- ✗ 行为不可预测

## 🔧 核心实现细节讨论

### 1. 参数传递策略

**问题**: `string_to_uniprot` 是大字典，如何传递给子进程？

#### 选项1.1: 直接传递（简单）
```python
def process_organism(tax_id, string_to_uniprot_dict, cache):
    # 每个进程接收完整字典的副本
    pass
```

**优点**: 简单直接
**缺点**: 内存开销（N个进程 = N个副本）

#### 选项1.2: 共享内存（优化）
```python
from multiprocessing import shared_memory

# 主进程创建共享内存
shm = shared_memory.ShareableList(...)
# 子进程访问共享内存
```

**优点**: 内存高效
**缺点**: 
- 复杂度高
- Python 3.8+ 才支持
- dict 不能直接共享（需要序列化）

#### 选项1.3: 部分传递（智能）
```python
# 只传递该物种相关的映射
relevant_mapping = filter_for_organism(string_to_uniprot, tax_id)
```

**优点**: 内存占用小
**缺点**: 需要预先知道哪些蛋白质属于哪个物种

**推荐**: 选项1.1（简单直接），除非遇到内存问题

### 2. 结果收集策略

**问题**: 如何收集和合并所有进程的结果？

#### 选项2.1: 内存中累积（当前方式）
```python
all_results = []
for future in as_completed(futures):
    organism_ints = future.result()
    all_results.extend(organism_ints)
self.string_ints = all_results
```

**优点**: 简单，与当前代码一致
**缺点**: 可能内存占用大

#### 选项2.2: 流式写入临时文件
```python
# 每个进程写入自己的文件
process_1.pkl, process_2.pkl, ...
# 最后合并
```

**优点**: 内存占用小
**缺点**: 
- I/O 开销
- 需要清理临时文件
- 增加复杂度

**推荐**: 选项2.1，除非数据量真的很大

### 3. 错误处理策略

**当前代码已有良好的错误处理** (line 730-751):
```python
except EOFError as e:
    logger.warning(...)
except zlib.error as e:
    logger.warning(...)
except (IndexError, ValueError) as e:
    logger.warning(...)
```

**并行化后的挑战**:
- 子进程的异常需要传回主进程
- 需要记录哪些物种失败了

#### 解决方案: 统一错误处理
```python
def process_organism(tax_id, ...):
    try:
        # 处理逻辑
        return (tax_id, organism_ints, None)  # 成功
    except Exception as e:
        return (tax_id, [], str(e))  # 失败，返回错误信息

# 主进程
for future in as_completed(futures):
    tax_id, organism_ints, error = future.result()
    if error:
        logger.warning(f"Failed {tax_id}: {error}")
        failed_organisms.append(tax_id)
    else:
        all_results.extend(organism_ints)
```

### 4. 进度显示

**当前**: 使用 `tqdm` 显示进度条

**并行化后的问题**: 
- `ProcessPoolExecutor` 的结果是无序的
- 进度条需要在主进程更新

**解决方案**: 使用 `as_completed` + `tqdm`
```python
with tqdm(total=len(futures), desc="Processing organisms") as pbar:
    for future in as_completed(futures):
        # 处理结果
        pbar.update(1)
        pbar.set_postfix(
            success=success_count,
            failed=failed_count,
            total_ints=len(all_results)
        )
```

## 📝 具体实现建议

### 推荐实现路径

**阶段1: 最小改动实现** (优先)
```python
def download_string_data(self, cache=False, n_workers=1):  # 默认串行
    """
    n_workers=1: 串行（向后兼容）
    n_workers>1: 并行
    n_workers=None: 自动检测
    """
    if n_workers == 1:
        # 保持原有串行代码不变
        self._download_string_serial(cache)
    else:
        # 新的并行代码
        self._download_string_parallel(cache, n_workers)
```

**优点**:
- ✓ 最小改动
- ✓ 完全向后兼容
- ✓ 容易测试和回滚

**阶段2: 优化实现** (可选)
- 共享内存优化
- 流式写入
- 更智能的负载均衡

### 关键代码位置

需要修改的地方：
1. **`__init__`** (line 92-166): 添加 `n_workers` 参数
2. **`download_string_data`** (line 651-759): 添加并行逻辑
3. 新增辅助函数:
   - `_process_single_organism()`: 静态方法，处理单个物种
   - `_download_string_serial()`: 原有串行逻辑
   - `_download_string_parallel()`: 新并行逻辑
   - `_get_optimal_workers()`: 自动检测最优进程数

### 函数设计草图

```python
class PPI:
    def __init__(self, ..., n_workers=None):
        # ...
        self.n_workers = n_workers  # None表示自动
    
    def download_string_data(self, cache=False):
        """保持向后兼容"""
        # 决定使用串行还是并行
        n_workers = self._determine_workers()
        
        if n_workers == 1:
            self._download_string_serial(cache)
        else:
            self._download_string_parallel(cache, n_workers)
    
    def _determine_workers(self):
        """自动检测或使用用户指定的值"""
        if self.n_workers is not None:
            return self.n_workers
        
        # 自动检测
        if len(self.tax_ids) < 10:  # 太少不值得并行
            return 1
        
        return self._get_optimal_workers()
    
    @staticmethod
    def _process_single_organism(tax_id, string_to_uniprot, cache):
        """处理单个物种 - 必须是静态方法或顶层函数"""
        # 这个函数会在子进程中运行
        pass
    
    def _download_string_serial(self, cache):
        """原有串行实现 - 保持不变"""
        # 当前的 for 循环代码
        pass
    
    def _download_string_parallel(self, cache, n_workers):
        """新的并行实现"""
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        # 准备数据
        string_to_uniprot_dict = dict(self.string_to_uniprot)
        
        # 并行处理
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {...}
            # 收集结果...
```

## 🤔 需要讨论的问题

### Q1: 默认行为
**选项A**: `n_workers=None` (自动检测)
- 优点: 自动优化
- 缺点: 可能出乎意料

**选项B**: `n_workers=1` (串行)
- 优点: 完全向后兼容
- 缺点: 用户需要手动启用

**你的倾向**: ?

### Q2: 自动检测策略
当 `n_workers=None` 时，如何决定进程数？

```python
def _get_optimal_workers(self):
    import multiprocessing
    import psutil
    
    cpu_count = multiprocessing.cpu_count()
    available_memory = psutil.virtual_memory().available / (1024**3)
    
    # 方案A: 保守 (避免影响其他任务)
    return min(cpu_count // 2, int(available_memory / 2))
    
    # 方案B: 激进 (最大化性能)
    return min(cpu_count - 2, int(available_memory / 1))
    
    # 方案C: 可配置
    return min(
        self.max_workers or cpu_count,
        int(available_memory / self.memory_per_worker)
    )
```

**你的倾向**: ?

### Q3: 是否需要在 `__init__` 中添加参数？

**选项A**: 添加到 `__init__`
```python
def __init__(self, ..., n_workers=None, max_workers=None):
    self.n_workers = n_workers
    self.max_workers = max_workers
```
- 优点: 一次设置，多次使用
- 缺点: 可能用户只下载一次

**选项B**: 只在 `download_string_data` 中
```python
def download_string_data(self, cache=False, n_workers=None):
    # n_workers 是局部参数
```
- 优点: 更灵活
- 缺点: 每次调用都要指定

**你的倾向**: ?

### Q4: 测试模式的处理

当前 `test_mode=True` 时只处理1个物种 (line 666)

**问题**: 测试模式下还需要并行吗？

**选项A**: 测试模式强制串行
```python
if self.test_mode:
    n_workers = 1
```

**选项B**: 允许测试并行代码
```python
# 即使 test_mode，也允许并行（但只有1个物种）
```

**你的倾向**: ?

### Q5: 失败物种的处理

当前代码会跳过失败的物种并继续

**问题**: 并行化后是否需要收集失败的物种？

**选项A**: 只记录日志（当前行为）
**选项B**: 返回失败列表，供用户重试
**选项C**: 添加自动重试机制

**你的倾向**: ?

## 📊 实施优先级

### 第一阶段 (必须)
1. ✓ 基本并行化实现
2. ✓ 向后兼容性
3. ✓ 错误处理
4. ✓ 进度显示

### 第二阶段 (应该)
5. ⭐ 自动检测最优进程数
6. ⭐ 完整的错误收集和报告
7. ⭐ 单元测试

### 第三阶段 (可选)
8. ○ 共享内存优化
9. ○ 流式写入优化
10. ○ 更细粒度的性能监控

## 💭 我的建议

基于分析，我建议：

1. **实施方案**: 方案A（添加 `n_workers` 参数）
   - 向后兼容，灵活性高

2. **参数位置**: `download_string_data(cache, n_workers=None)`
   - 不改 `__init__`，保持简单

3. **默认行为**: `n_workers=None` → 自动检测
   - 但提供 `n_workers=1` 选项回退到串行

4. **错误处理**: 收集失败列表 + 日志
   - 方便用户了解哪些物种失败了

5. **实施顺序**: 
   - 先实现基本版本
   - 验证正确性和性能
   - 再考虑优化

## 📋 下一步

在我开始编码之前，请告诉我：

1. ✅ 你同意方案A（添加 `n_workers` 参数）吗？
2. ✅ 默认 `n_workers=None`（自动）还是 `n_workers=1`（串行）？
3. ✅ 是否需要在 `__init__` 中添加参数，还是只在方法中？
4. ✅ 测试模式下要并行吗？
5. ✅ 还有其他特殊需求或顾虑吗？

**讨论完这些，我再开始写代码！** 🚀
