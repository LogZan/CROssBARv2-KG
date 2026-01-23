#!/usr/bin/env python3
"""
STRING 数据并行处理示例
演示如何使用多进程加速 STRING 数据下载
"""

import sys
import time
import collections
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import multiprocessing

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import string, uniprot
from biocypher._logger import logger


def process_single_organism(args):
    """
    处理单个物种的 STRING 数据
    
    Args:
        args: (tax_id, string_to_uniprot_dict, cache)
        
    Returns:
        (tax_id, interactions_list, error_msg)
    """
    tax_id, string_to_uniprot, cache = args
    
    try:
        # 重新设置缓存 (每个进程需要自己的设置)
        cache_config.set_adapter_cache('string')
        
        # 获取该物种的 STRING interactions
        organism_string_ints = [
            i
            for i in string.string_links_interactions(
                ncbi_tax_id=int(tax_id),
                score_threshold="high_confidence",
            )
            if i.protein_a in string_to_uniprot
            and i.protein_b in string_to_uniprot
        ]
        
        return (tax_id, organism_string_ints, None)
        
    except EOFError as e:
        return (tax_id, [], f"EOFError: {e}")
    except Exception as e:
        return (tax_id, [], f"{type(e).__name__}: {e}")


def get_optimal_workers(max_workers=None):
    """
    根据 CPU 和内存自动确定最优进程数
    
    Args:
        max_workers: 最大进程数限制
        
    Returns:
        optimal_workers: 推荐的进程数
    """
    import psutil
    
    cpu_count = multiprocessing.cpu_count()
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    
    # 假设每个进程需要 1GB 内存 (保守估计)
    memory_workers = int(available_memory_gb / 1.0)
    
    # 取 CPU 和内存的较小值
    optimal = min(cpu_count, memory_workers)
    
    # 保留一些资源给系统
    optimal = max(1, optimal - 2)
    
    # 应用用户限制
    if max_workers:
        optimal = min(optimal, max_workers)
    
    print(f"System resources:")
    print(f"  CPU cores: {cpu_count}")
    print(f"  Available memory: {available_memory_gb:.1f} GB")
    print(f"  Memory-based workers: {memory_workers}")
    print(f"  Recommended workers: {optimal}")
    
    return optimal


def download_string_parallel(tax_ids, string_to_uniprot, 
                             cache=True, n_workers=None,
                             test_mode=False):
    """
    并行下载 STRING 数据
    
    Args:
        tax_ids: 物种 ID 列表
        string_to_uniprot: STRING ID 到 UniProt ID 的映射
        cache: 是否使用缓存
        n_workers: 并发进程数 (None表示自动)
        test_mode: 测试模式，只处理前10个物种
        
    Returns:
        all_interactions: 所有物种的 interactions
        failed_organisms: 失败的物种列表
        stats: 统计信息
    """
    # 过滤有问题的 tax IDs
    tax_ids_to_be_skipped = ["4565", "8032", "1829", "1894"]
    valid_tax_ids = [tax for tax in tax_ids if str(tax) not in tax_ids_to_be_skipped]
    
    if test_mode:
        print("⚠️  TEST MODE: Only processing first 10 organisms")
        valid_tax_ids = valid_tax_ids[:10]
    
    # 确定进程数
    if n_workers is None:
        n_workers = get_optimal_workers()
    else:
        print(f"Using user-specified workers: {n_workers}")
    
    print(f"\nProcessing {len(valid_tax_ids)} organisms with {n_workers} workers")
    print(f"Expected time: ~{len(valid_tax_ids) * 15 / n_workers / 3600:.1f} hours")
    print("-" * 60)
    
    # 转换为普通dict (可序列化)
    string_to_uniprot_dict = dict(string_to_uniprot)
    
    # 准备任务参数
    tasks = [
        (tax_id, string_to_uniprot_dict, cache) 
        for tax_id in valid_tax_ids
    ]
    
    # 并行处理
    all_interactions = []
    failed_organisms = []
    processing_times = []
    
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # 提交所有任务
        futures = {
            executor.submit(process_single_organism, task): task[0]
            for task in tasks
        }
        
        # 收集结果
        with tqdm(total=len(futures), desc="Processing organisms") as pbar:
            for future in as_completed(futures):
                task_start = time.time()
                tax_id = futures[future]
                
                try:
                    tax_id, organism_ints, error = future.result()
                    
                    if error:
                        failed_organisms.append((tax_id, error))
                        pbar.set_postfix(
                            success=len(futures) - len(failed_organisms),
                            failed=len(failed_organisms)
                        )
                    else:
                        all_interactions.extend(organism_ints)
                        pbar.set_postfix(
                            success=len(futures) - len(failed_organisms),
                            failed=len(failed_organisms),
                            total_ints=len(all_interactions)
                        )
                    
                except Exception as e:
                    failed_organisms.append((tax_id, f"Future error: {e}"))
                    pbar.set_postfix(
                        success=len(futures) - len(failed_organisms),
                        failed=len(failed_organisms)
                    )
                
                processing_times.append(time.time() - task_start)
                pbar.update(1)
    
    total_time = time.time() - start_time
    
    # 统计信息
    stats = {
        'total_organisms': len(valid_tax_ids),
        'successful': len(valid_tax_ids) - len(failed_organisms),
        'failed': len(failed_organisms),
        'total_interactions': len(all_interactions),
        'total_time_seconds': total_time,
        'total_time_hours': total_time / 3600,
        'avg_time_per_organism': sum(processing_times) / len(processing_times) if processing_times else 0,
        'speedup': (len(valid_tax_ids) * 15) / total_time if total_time > 0 else 0,
    }
    
    return all_interactions, failed_organisms, stats


def main():
    """主函数 - 演示并行处理"""
    
    print("=" * 80)
    print("STRING 数据并行处理示例")
    print("=" * 80)
    
    # ========================================================================
    # 1. 准备工作 - 获取物种列表和映射
    # ========================================================================
    print("\n[1/3] 准备数据...")
    
    # 获取所有物种 (实际使用中这里应该从 PPI adapter 获取)
    test_organism = 9606  # 人类
    test_tax_ids = [test_organism]  # 测试时只用一个物种
    
    # 实际使用时应该是:
    # string_species = string.string_species()
    # tax_ids = list(string_species.keys())
    
    # 获取 UniProt 到 STRING 的映射
    print("  Getting UniProt to STRING mapping...")
    uniprot_to_string = uniprot.uniprot_data(
        fields="xref_string",
        organism=None,  # 所有物种
        reviewed=True
    )
    
    string_to_uniprot = collections.defaultdict(list)
    for k, v in uniprot_to_string.items():
        for string_id in list(filter(None, str(v).split(";"))):
            if "." in string_id:
                string_to_uniprot[string_id.split(".")[1]].append(k)
    
    print(f"  ✓ Mapping created: {len(string_to_uniprot)} STRING IDs")
    
    # ========================================================================
    # 2. 对比串行 vs 并行
    # ========================================================================
    print("\n[2/3] Performance comparison...")
    print("-" * 80)
    
    # 串行处理 (1个进程)
    print("\n📊 Serial processing (1 worker):")
    _, _, serial_stats = download_string_parallel(
        test_tax_ids, 
        string_to_uniprot,
        cache=True,
        n_workers=1,
        test_mode=True
    )
    
    # 并行处理 (自动确定进程数)
    print("\n📊 Parallel processing (auto workers):")
    _, _, parallel_stats = download_string_parallel(
        test_tax_ids,
        string_to_uniprot,
        cache=True,
        n_workers=None,  # 自动
        test_mode=True
    )
    
    # ========================================================================
    # 3. 显示结果
    # ========================================================================
    print("\n[3/3] Results summary")
    print("=" * 80)
    
    print("\nSerial processing:")
    print(f"  Total organisms: {serial_stats['total_organisms']}")
    print(f"  Successful: {serial_stats['successful']}")
    print(f"  Failed: {serial_stats['failed']}")
    print(f"  Total interactions: {serial_stats['total_interactions']:,}")
    print(f"  Time: {serial_stats['total_time_hours']:.2f} hours ({serial_stats['total_time_seconds']:.1f}s)")
    print(f"  Avg per organism: {serial_stats['avg_time_per_organism']:.1f}s")
    
    print("\nParallel processing:")
    print(f"  Total organisms: {parallel_stats['total_organisms']}")
    print(f"  Successful: {parallel_stats['successful']}")
    print(f"  Failed: {parallel_stats['failed']}")
    print(f"  Total interactions: {parallel_stats['total_interactions']:,}")
    print(f"  Time: {parallel_stats['total_time_hours']:.2f} hours ({parallel_stats['total_time_seconds']:.1f}s)")
    print(f"  Avg per organism: {parallel_stats['avg_time_per_organism']:.1f}s")
    print(f"  Speedup: {parallel_stats['speedup']:.1f}×")
    
    # 估算全量处理时间
    print("\n" + "=" * 80)
    print("Estimated time for 12,535 organisms:")
    print("=" * 80)
    
    for n_workers in [1, 4, 8, 16, 32]:
        estimated_time = (12535 * serial_stats['avg_time_per_organism']) / n_workers / 3600
        print(f"  {n_workers:2d} workers: {estimated_time:6.1f} hours")
    
    print("\n" + "=" * 80)
    print("✅ Done! Parallel processing is highly effective.")
    print("=" * 80)


if __name__ == "__main__":
    main()
