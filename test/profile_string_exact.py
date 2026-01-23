#!/usr/bin/env python3
"""
精确复现 ppi_adapter.py 中的 STRING 调用
找出为什么每个物种需要 15.82 秒
"""

import sys
import time
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import string

print("\n" + "="*80)
print("精确复现 PPI Adapter 的 STRING 调用")
print("="*80)

# 这是 ppi_adapter.py 第712-715行的实际调用
TEST_ORGANISM = 9606

print(f"\n测试物种: {TEST_ORGANISM} (人类)")
print("调用参数: ncbi_tax_id={}, score_threshold='high_confidence'".format(TEST_ORGANISM))
print("-" * 80)

# ========================================================================
# 步骤1: 先测试完整生成器转列表的耗时
# ========================================================================
print("\n【步骤1】调用 string.string_links_interactions 并转为列表")
print("  这是 ppi_adapter.py line 710-718 实际做的事")

start = time.time()
try:
    # 完全复现 ppi_adapter.py 的调用
    organism_string_ints = [
        i
        for i in string.string_links_interactions(
            ncbi_tax_id=int(TEST_ORGANISM),
            score_threshold="high_confidence",
        )
    ]
    elapsed = time.time() - start
    print(f"  ✓ 完成: {len(organism_string_ints):,} 条记录")
    print(f"  ⏱  总耗时: {elapsed:.2f}s")
    print(f"  ⏱  平均速度: {len(organism_string_ints)/elapsed:.0f} records/s")
except Exception as e:
    elapsed = time.time() - start
    print(f"  ✗ 失败: {e}")
    import traceback
    traceback.print_exc()
    print(f"  ⏱  耗时: {elapsed:.2f}s")

# ========================================================================
# 步骤2: 分解各个子步骤的耗时
# ========================================================================
print("\n【步骤2】分解子步骤耗时")

# 2.1: 读取 physical interactions (这会在 string_links_interactions 内部调用)
print("\n  2.1 读取 physical interactions (score=0)")
start = time.time()
physical_ints = list(string.string_physical_interactions(
    ncbi_tax_id=TEST_ORGANISM,
    score_threshold=0  # 这是 string_links_interactions 内部使用的参数
))
elapsed_physical = time.time() - start
print(f"      ✓ {len(physical_ints):,} 条记录")
print(f"      ⏱  {elapsed_physical:.2f}s")

# 2.2: 构建 physical score 字典
print("\n  2.2 构建 physical score 字典")
start = time.time()
phy_links = dict(
    ((i.protein_a, i.protein_b), i.combined_score)
    for i in physical_ints
)
elapsed_dict = time.time() - start
print(f"      ✓ {len(phy_links):,} 个键值对")
print(f"      ⏱  {elapsed_dict:.2f}s")

# 2.3: 只统计 links 文件有多少行(不过滤)
print("\n  2.3 统计 links 文件行数")
from pypath.resources import urls
from pypath.share import curl
url = urls.urls['string']['links'] % TEST_ORGANISM
start = time.time()
c = curl.Curl(url, silent=False, large=True)
_ = next(c.result)  # 跳过header
line_count = sum(1 for _ in c.result)
elapsed_count = time.time() - start
print(f"      ✓ {line_count:,} 行")
print(f"      ⏱  {elapsed_count:.2f}s")
print(f"      ⏱  速度: {line_count/elapsed_count:.0f} lines/s")

# ========================================================================
# 步骤3: 测试不同的过滤条件
# ========================================================================
print("\n【步骤3】测试过滤的影响")

# 3.1: 过滤 high_confidence (score >= 700)
print("\n  3.1 过滤 high_confidence (score >= 700)")
url = urls.urls['string']['links'] % TEST_ORGANISM
start = time.time()
c = curl.Curl(url, silent=False, large=True)
_ = next(c.result)
filtered_count = 0
for l in c.result:
    l = l.strip().split(' ')
    if int(l[9]) >= 700:  # high_confidence threshold
        filtered_count += 1
elapsed_filter = time.time() - start
print(f"      ✓ {filtered_count:,} 条符合条件")
print(f"      ⏱  {elapsed_filter:.2f}s")
print(f"      ⏱  速度: {line_count/elapsed_filter:.0f} lines/s")

# 3.2: 完整解析(包括查询字典)
print("\n  3.2 完整解析(解析+过滤+查字典)")
url = urls.urls['string']['links'] % TEST_ORGANISM
start = time.time()
c = curl.Curl(url, silent=False, large=True)
_ = next(c.result)
parsed_count = 0
for l in c.result:
    l = l.strip().split(' ')
    prot_a_id = l[0].split('.')[1]
    prot_b_id = l[1].split('.')[1]
    
    if int(l[9]) >= 700:
        # 查询字典
        phy_score = phy_links.get((prot_a_id, prot_b_id), None)
        parsed_count += 1

elapsed_parse = time.time() - start
print(f"      ✓ {parsed_count:,} 条记录")
print(f"      ⏱  {elapsed_parse:.2f}s")
print(f"      ⏱  速度: {line_count/elapsed_parse:.0f} lines/s")

# ========================================================================
# 总结
# ========================================================================
print("\n" + "="*80)
print("性能分解总结 (人类)")
print("="*80)

print(f"\n各步骤耗时:")
print(f"  1. Physical interactions (score=0):  {elapsed_physical:>8.2f}s  ({len(physical_ints):>10,} records)")
print(f"  2. 构建 physical score 字典:         {elapsed_dict:>8.2f}s  ({len(phy_links):>10,} entries)")
print(f"  3. 读取 links 文件:                 {elapsed_count:>8.2f}s  ({line_count:>10,} lines)")
print(f"  4. 过滤 high_confidence:            {elapsed_filter:>8.2f}s  ({filtered_count:>10,} filtered)")
print(f"  5. 完整解析(含字典查询):             {elapsed_parse:>8.2f}s  ({parsed_count:>10,} records)")

total_estimated = elapsed_physical + elapsed_dict + elapsed_parse
print(f"\n估计总耗时: {total_estimated:.2f}s")
print(f"实际测试值: {elapsed:.2f}s" if 'elapsed' in locals() else "")

print(f"\n🔍 主要瓶颈分析:")
bottlenecks = [
    ("Physical interactions 读取", elapsed_physical),
    ("Links 文件读取", elapsed_count if elapsed_count > elapsed_filter else elapsed_filter),
    ("字典查询和解析", elapsed_parse - elapsed_count if elapsed_parse > elapsed_count else 0),
]
bottlenecks.sort(key=lambda x: x[1], reverse=True)

for i, (name, time_cost) in enumerate(bottlenecks, 1):
    pct = (time_cost / total_estimated * 100) if total_estimated > 0 else 0
    print(f"  {i}. {name:<30} {time_cost:>8.2f}s ({pct:>5.1f}%)")

print("\n" + "="*80)
print("结论")
print("="*80)
print(f"""
对于人类(9606):
- Physical interactions: {len(physical_ints):,} 条记录
- Links total: {line_count:,} 行
- High confidence links: {filtered_count:,} 条

每个物种的处理时间 ≈ {total_estimated:.1f}s，主要耗时在:
1. 读取和解析 physical interactions ({elapsed_physical:.1f}s)
2. 读取和解析 links 文件 ({elapsed_count:.1f}s)

即使数据已缓存，每个物种仍需要:
- 解压 .gz 文件
- 解析文本数据
- 构建数据结构

对于 12,535 个物种:
- 总时间 ≈ 12,535 × {total_estimated:.1f}s = {12535 * total_estimated / 3600:.1f} 小时
""")

print("\n" + "="*80)
