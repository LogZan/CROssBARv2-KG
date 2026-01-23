#!/usr/bin/env python3
"""
深度测试：分析 STRING 数据读取的真正瓶颈
测试 string.string_links_interactions() 的各个步骤耗时
"""

import sys
import time
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import string
from pypath.share import curl

print("\n" + "="*80)
print("STRING 数据读取性能深度分析")
print("="*80)

# 测试单个物种
TEST_ORGANISM = 9606  # 人类

print(f"\n测试物种: {TEST_ORGANISM} (人类)")
print("-" * 80)

# =============================================================================
# 测试 1: 只读取 physical interactions (在 string_links_interactions 中被调用)
# =============================================================================
print("\n【测试 1】读取 physical interactions (high_confidence)")
start = time.time()
try:
    physical_ints = list(string.string_physical_interactions(
        ncbi_tax_id=TEST_ORGANISM,
        score_threshold='high_confidence'
    ))
    elapsed_physical = time.time() - start
    print(f"  ✓ 完成: {len(physical_ints):,} 条记录")
    print(f"  ⏱  耗时: {elapsed_physical:.2f}s")
except Exception as e:
    elapsed_physical = time.time() - start
    print(f"  ✗ 失败: {e}")
    print(f"  ⏱  耗时: {elapsed_physical:.2f}s")

# =============================================================================
# 测试 2: 读取 physical interactions (score_threshold=0, 这是实际调用的)
# =============================================================================
print("\n【测试 2】读取 physical interactions (score_threshold=0)")
print("  这是 string_links_interactions 内部实际调用的参数!")
start = time.time()
try:
    physical_ints_all = list(string.string_physical_interactions(
        ncbi_tax_id=TEST_ORGANISM,
        score_threshold=0  # 获取所有分数的记录!
    ))
    elapsed_physical_all = time.time() - start
    print(f"  ✓ 完成: {len(physical_ints_all):,} 条记录")
    print(f"  ⏱  耗时: {elapsed_physical_all:.2f}s")
except Exception as e:
    elapsed_physical_all = time.time() - start
    print(f"  ✗ 失败: {e}")
    print(f"  ⏱  耗时: {elapsed_physical_all:.2f}s")

# =============================================================================
# 测试 3: 只读取 links (不包含 physical score)
# =============================================================================
print("\n【测试 3】读取 links (high_confidence, physical_interaction_score=False)")
start = time.time()
try:
    # 模拟关闭 physical_interaction_score 的情况
    links_count = 0
    for interaction in string.string_links_interactions(
        ncbi_tax_id=TEST_ORGANISM,
        score_threshold='high_confidence',
        physical_interaction_score=False  # 不获取 physical score
    ):
        links_count += 1
        if links_count >= 100:  # 只统计前100个看速度
            break
    
    elapsed_links_sample = time.time() - start
    print(f"  ✓ 完成前 {links_count:,} 条记录")
    print(f"  ⏱  耗时: {elapsed_links_sample:.2f}s")
except Exception as e:
    elapsed_links_sample = time.time() - start
    print(f"  ✗ 失败: {e}")
    print(f"  ⏱  耗时: {elapsed_links_sample:.2f}s")

# =============================================================================
# 测试 4: 完整调用 string_links_interactions (实际使用的方式)
# =============================================================================
print("\n【测试 4】完整调用 string_links_interactions (实际场景)")
print("  参数: score_threshold='high_confidence', physical_interaction_score=True")
start = time.time()
try:
    full_ints = list(string.string_links_interactions(
        ncbi_tax_id=TEST_ORGANISM,
        score_threshold='high_confidence',
        physical_interaction_score=True
    ))
    elapsed_full = time.time() - start
    print(f"  ✓ 完成: {len(full_ints):,} 条记录")
    print(f"  ⏱  耗时: {elapsed_full:.2f}s")
except Exception as e:
    elapsed_full = time.time() - start
    print(f"  ✗ 失败: {e}")
    print(f"  ⏱  耗时: {elapsed_full:.2f}s")

# =============================================================================
# 分析结果
# =============================================================================
print("\n" + "="*80)
print("性能分析结果")
print("="*80)

print(f"\n各步骤耗时对比:")
print(f"  Physical interactions (high_conf):  {elapsed_physical:.2f}s")
print(f"  Physical interactions (score=0):    {elapsed_physical_all:.2f}s  ⚠️")
print(f"  Links only (前100条):                {elapsed_links_sample:.2f}s")
print(f"  完整流程:                             {elapsed_full:.2f}s")

print(f"\n🔍 关键发现:")

if elapsed_physical_all > 10:
    print(f"\n⚠️  【重大发现】physical interactions (score=0) 耗时 {elapsed_physical_all:.2f}s!")
    print(f"    这是因为 string_links_interactions 内部调用:")
    print(f"    string_physical_interactions(ncbi_tax_id, score_threshold=0)")
    print(f"    score_threshold=0 意味着返回 ALL 记录，不做过滤!")
    print(f"\n    即使用户只需要 high_confidence 的 links,")
    print(f"    但为了获取 physical_combined_score,")
    print(f"    函数会先读取 ALL physical interactions!")

if elapsed_full > 15:
    print(f"\n⚠️  【瓶颈确认】完整流程耗时 {elapsed_full:.2f}s")
    print(f"    接近实际观察到的 15.82s/it")
    print(f"\n    瓶颈分解:")
    print(f"    1. 读取 physical interactions: ~{elapsed_physical_all:.1f}s")
    print(f"    2. 构建 physical score 字典")
    print(f"    3. 读取 links 并查询字典: 剩余时间")

print("\n" + "="*80)
print("深入原因分析")
print("="*80)

print("""
问题根源在 pypath/inputs/string.py 的 string_links_interactions 函数:

第45-54行:
    if physical_interaction_score:
        phy_links = dict(
            ((i.protein_a, i.protein_b), i.combined_score)
            for i in
            string_physical_interactions(
                ncbi_tax_id = ncbi_tax_id,
                score_threshold = 0,  # ⚠️ 这里是 0!
            )
        )

即使你只需要 high_confidence 的 links，
但为了填充 physical_combined_score 字段，
函数会先调用 string_physical_interactions(score_threshold=0)
这会读取并处理该物种的 ALL physical interactions!

对于人类:
- physical.links.detailed 文件可能有几百万条记录
- 需要解压、解析、构建字典
- 即使有缓存，解压和解析仍需要时间

这就是为什么每个物种需要 15+ 秒的真正原因!
""")

print("\n" + "="*80)
print("验证缓存读取速度")
print("="*80)

# 测试纯文件读取速度
import os
cache_dir = "/GenSIvePFS/users/data/pypath_cache/string"
files = os.listdir(cache_dir)
print(f"\n缓存目录中共有 {len(files):,} 个文件")

# 找一个 physical.links 文件测试读取速度
physical_files = [f for f in files if 'physical.links' in f and '9606' in f]
if physical_files:
    test_file = os.path.join(cache_dir, physical_files[0])
    print(f"\n测试文件: {physical_files[0]}")
    print(f"文件大小: {os.path.getsize(test_file) / 1024 / 1024:.1f} MB")
    
    # 测试解压和读取
    import gzip
    start = time.time()
    line_count = 0
    with gzip.open(test_file, 'rt') as f:
        for line in f:
            line_count += 1
    elapsed_read = time.time() - start
    
    print(f"读取行数: {line_count:,}")
    print(f"读取耗时: {elapsed_read:.2f}s")
    print(f"读取速度: {line_count / elapsed_read:.0f} lines/s")

print("\n" + "="*80)
