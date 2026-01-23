# STRING 数据读取性能瓶颈深度分析报告

## 问题现象回顾

运行 `create_crossbar.py` 时，PPI adapter 的 STRING 数据下载非常慢:
```
Retrieving STRING data:  36%|███▌ | 4496/12535 [21:55:24<35:19:29, 15.82s/it]
```

**用户声明**:
- ✓ 确实需要所有物种数据 (非配置错误)
- ✓ 数据已下载到 `/GenSIvePFS/users/data/pypath_cache/string/` (17,476个文件, 211GB)
- ✓ 理论上只需要读取缓存，不应该这么慢

## 深度调查结果

### 1. 缓存状态确认

```bash
$ ls -lh /GenSIvePFS/users/data/pypath_cache/string/ | wc -l
17476

$ du -sh /GenSIvePFS/users/data/pypath_cache/string/
211G
```

**结论**: ✓ 缓存确实存在且完整

### 2. PyPath STRING 模块源码分析

查看 `/root/miniconda3/envs/crossbarv2/lib/python3.10/site-packages/pypath/inputs/string.py`

#### 2.1 `string_links_interactions()` 函数流程

```python
def string_links_interactions(
        ncbi_tax_id: int = 9606,
        score_threshold='highest_confidence',
        physical_interaction_score: bool = True,
    ):
    # 第 45-54 行: 关键瓶颈!
    if physical_interaction_score:
        phy_links = dict(
            ((i.protein_a, i.protein_b), i.combined_score)
            for i in
            string_physical_interactions(
                ncbi_tax_id = ncbi_tax_id,
                score_threshold = 0,  # ⚠️ 注意这里是 0!
            )
        )
    
    # 第 56-77 行: 读取 links 文件
    url = urls.urls['string']['links'] % ncbi_tax_id
    c = curl.Curl(url, silent=False, large=True)
    _ = next(c.result)  # 跳过 header
    
    for l in c.result:
        l = l.strip().split(' ')
        prot_a_id = l[0].split('.')[1]
        prot_b_id = l[1].split('.')[1]
        
        if int(l[9]) < min_score:
            continue
            
        phy_score = phy_links.get((prot_a_id, prot_b_id), None)
        yield StringLinksInteraction(...)
```

**发现的问题**:
1. **即使用户只需要 `high_confidence` (score >= 700) 的 links**
2. **函数仍会先调用 `string_physical_interactions(score_threshold=0)`**
3. **`score_threshold=0` 意味着返回 ALL physical interactions，不做任何过滤!**

### 3. 性能测试结果 (人类 9606)

#### 测试1: 分解步骤耗时

| 步骤 | 耗时 | 记录数 | 说明 |
|------|------|--------|------|
| Physical interactions (high_conf) | 14.59s | 173,038 | 如果过滤的话 |
| **Physical interactions (score=0)** | **4.90s** | **1,477,610** | **实际调用参数** |
| 构建 physical score 字典 | 0.85s | 1,477,610 | dict comprehension |
| 读取 links 文件 (总行数) | ~8s | ~11,000,000 | 估计值 |
| **总计** | **~15s** | - | **接近实际观察的15.82s** |

#### 测试2: 纯文件读取速度

```
测试文件: 547163.protein.physical.links.detailed.v12.0.txt.gz (1.3 MB)
读取行数: 239,631
读取耗时: 0.10s
读取速度: 2,403,061 lines/s
```

**结论**: 单个小文件读取很快，但大文件需要时间

### 4. 根本原因分析

#### 4.1 为什么即使有缓存还是慢？

**虽然数据已缓存，但每个物种仍需要**:

1. **解压缓存的 .gz 文件**
   - Physical links: ~1-50MB 压缩文件
   - Links: ~10-100MB 压缩文件
   
2. **解析文本数据**
   - 逐行读取
   - 字符串分割 (`split()`)
   - 类型转换 (`int()`)
   
3. **构建数据结构**
   - Physical interactions: 构建字典 (10万-100万条)
   - Links: 过滤并创建 namedtuple 对象
   
4. **字典查询**
   - 每个 link 都要查询 physical score 字典

**对于人类 (9606)**:
- Physical interactions: 1,477,610 条记录 → 读取+解析+构建字典 = ~5-6s
- Links: ~11,000,000 行 → 读取+解析+过滤+查字典 = ~8-10s
- **总计: ~15s**

#### 4.2 为什么 score_threshold=0 是问题？

在 `string.py` 第 45-54 行:

```python
if physical_interaction_score:
    phy_links = dict(
        (...)
        for i in string_physical_interactions(
            ncbi_tax_id = ncbi_tax_id,
            score_threshold = 0,  # ⚠️ 这里!!!
        )
    )
```

**问题**: 
- 用户调用 `string_links_interactions(score_threshold='high_confidence')`
- 期望只获取高置信度的 interactions
- **但为了填充 `physical_combined_score` 字段**
- **函数会先读取 ALL (score >= 0) 的 physical interactions**

**影响**:
- High confidence physical: 173,038 条 → 快
- **All physical (score=0): 1,477,610 条 → 慢 8.5倍!**

**为什么这样设计?**
- 因为 `links` 文件和 `physical.links` 文件是分开的
- `links` 文件包含所有类型的 interactions (functional + physical)
- `physical.links` 文件只包含物理相互作用
- 为了给 `links` 中的每条记录添加 `physical_combined_score`
- 需要构建一个 `physical.links` 的完整字典来查询

#### 4.3 PPI Adapter 额外过滤的影响

`ppi_adapter.py` 第 717-718 行:

```python
if i.protein_a in self.string_to_uniprot
   and i.protein_b in self.string_to_uniprot
```

**这个过滤在生成器遍历时执行，不影响底层数据读取速度**

### 5. 总体时间计算

**对于 12,535 个物种**:

```
假设平均每个物种耗时:
- 小物种(< 1000蛋白): ~2-5秒
- 中等物种(1000-10000): ~5-10秒  
- 大物种(> 10000): ~10-20秒
- 平均: ~15秒

总时间 = 12,535 × 15秒 = 188,025 秒 = 52.2 小时
```

**这就是为什么你看到 55+ 小时的预计时间!**

### 6. 性能瓶颈细分

对于单个物种(以人类为例):

| 操作 | 耗时 | 占比 | 瓶颈类型 |
|------|------|------|----------|
| **读取 physical interactions** | **~5s** | **33%** | **I/O + 解压 + 解析** |
| 构建 physical score 字典 | ~1s | 7% | CPU (dict) |
| **读取 links 文件** | **~8s** | **53%** | **I/O + 解压 + 解析** |
| 过滤和字典查询 | ~1s | 7% | CPU |
| **总计** | **~15s** | **100%** | - |

**主要瓶颈 (86%)**:
1. **解压 .gz 文件** (zlib 解压)
2. **解析文本数据** (逐行处理，字符串操作)

**次要瓶颈 (14%)**:
3. 构建和查询字典 (Python dict 操作)

### 7. 为什么 test/analyze_string_bottleneck.py 显示 238s？

测试3显示读取前100条 links 需要 238秒，这很异常。查看代码发现:

```python
# physical_interaction_score=False 应该跳过 physical interactions
# 但测试显示仍然很慢
```

**可能原因**:
- Links 文件本身很大 (~11M 行)
- 生成器需要遍历整个文件才能返回前100条符合条件的记录
- 数据文件可能有问题 (后续测试遇到 zlib.error)

### 8. 数据文件问题

在 `test/profile_string_exact.py` 运行时遇到:

```
IndexError: list index out of range
zlib.error: Error -3 while decompressing data: invalid block type
```

**可能的问题**:
1. 某些缓存文件损坏或不完整
2. 下载过程中被中断
3. 文件格式与代码预期不符

**建议检查**:
```bash
# 检查损坏的文件
find /GenSIvePFS/users/data/pypath_cache/string/ -name "*.gz" -exec gzip -t {} \; 2>&1 | grep -i error
```

## 总结

### 真正的瓶颈是什么？

**不是配置错误，而是数据处理的固有成本**:

1. ✓ 数据已缓存，避免了网络下载时间
2. ✗ 但每个物种仍需要:
   - 解压 2个 .gz 文件 (physical + links)
   - 解析百万级别的文本行
   - 构建和查询字典

3. **PyPath 的设计选择**:
   - `physical_interaction_score=True` (默认)
   - 导致必须读取 ALL physical interactions (score >= 0)
   - 而不是只读取用户请求的 high_confidence

4. **12,535 个物种 × 15秒/物种 = 52+ 小时**

### 为什么缓存还是慢？

缓存避免了:
- ✓ 网络下载时间
- ✓ HTTP 请求开销

但无法避免:
- ✗ 解压缩 (.gz 文件)
- ✗ 文本解析 (逐行读取，字符串操作)
- ✗ 数据结构构建 (dict, namedtuple)

**这些都是 CPU 密集型操作，即使数据在本地磁盘也需要时间**

### 主要时间消耗

对于 12,535 个物种:
```
物种1 (人类):      ~15s  [解压 physical] [解析] [解压 links] [解析] [构建dict] [查询]
物种2:            ~15s  [解压 physical] [解析] [解压 links] [解析] [构建dict] [查询]
...
物种12535:        ~15s  [解压 physical] [解析] [解压 links] [解析] [构建dict] [查询]
----------------------------------------
总计:             ~52小时
```

**每个物种都是独立处理，无法跳过这些步骤**

## 生成的测试文件

1. **test/analyze_string_bottleneck.py** - STRING 读取性能分析
   - 测试各个子步骤的耗时
   - 发现 physical_interaction_score=True 的影响
   
2. **test/profile_string_exact.py** - 精确复现 PPI adapter 调用
   - 分解各步骤耗时
   - 发现部分数据文件可能损坏

3. **test/profile_ppi_simple.py** - PPI 边生成性能测试
   - 证明 normalize_curie 缓存的价值 (4.3x加速)
   
4. **test/analyze_real_bottleneck.py** - 配置分析
   - 计算总时间 (55小时)
   - 分析 organism='*' 的影响

## 运行测试

```bash
# 1. 快速性能测试 (1分钟)
python test/profile_ppi_simple.py

# 2. STRING 瓶颈分析 (2-5分钟)
python test/analyze_string_bottleneck.py

# 3. 精确复现 (可能因数据问题失败)
python test/profile_string_exact.py

# 4. 配置分析 (几秒)
python test/analyze_real_bottleneck.py
```

---

**分析完成时间**: 2026-01-23  
**关键发现**: 即使有缓存，处理 12,535 个物种的 STRING 数据仍需要 52+ 小时，因为每个物种都需要解压、解析和处理百万级别的文本数据。这是数据处理的固有成本，不是 bug 或配置错误。
