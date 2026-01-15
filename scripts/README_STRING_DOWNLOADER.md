# STRING Data Downloader

专门用于下载 STRING 数据库的独立脚本，具有以下特性：

## 功能特性

✅ **断点续传**: 自动保存下载进度，可以随时中断和恢复
✅ **错误重试**: 自动重试失败的下载（默认最多3次，指数退避）
✅ **并行下载**: 使用多线程并行下载（默认8个线程）
✅ **进度跟踪**: 实时显示下载进度和成功/失败统计
✅ **缓存复用**: 使用与 pypath 相同的缓存机制
✅ **日志记录**: 详细的日志记录到文件和控制台

## 使用方法

### 基本用法（使用默认设置）

```bash
cd /GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG
python scripts/download_string_data.py
```

### 高级用法

**使用 16 个并行线程:**
```bash
python scripts/download_string_data.py --workers 16
```

**设置最大重试次数为 5:**
```bash
python scripts/download_string_data.py --max-retries 5
```

**从头开始下载（忽略之前的进度）:**
```bash
python scripts/download_string_data.py --no-resume
```

**组合使用:**
```bash
python scripts/download_string_data.py --workers 16 --max-retries 5
```

## 参数说明

- `--workers N`: 并行下载线程数（默认: 8）
- `--max-retries N`: 每个物种的最大重试次数（默认: 3）
- `--cache-dir PATH`: 缓存目录路径（默认: /GenSIvePFS/users/data/pypath_cache）
- `--no-resume`: 忽略之前的进度，从头开始
- `--retry-failed`: 重试之前失败的下载（默认: 启用）

## 进度文件

下载过程中会生成以下文件：

1. **进度文件**: `/GenSIvePFS/users/data/pypath_cache/string_download_progress.json`
   - 记录已成功下载的物种 ID
   - 可以用于断点续传

2. **失败列表**: `/GenSIvePFS/users/data/pypath_cache/string_download_failed.json`
   - 记录下载失败的物种 ID 和错误信息
   - 下次运行时会自动重试

3. **日志文件**: `/GenSIvePFS/users/data/pypath_cache/string_download.log`
   - 详细的下载日志

## 中断和恢复

**中断下载:**
- 按 `Ctrl+C` 可以随时中断
- 进度会自动保存

**恢复下载:**
- 再次运行相同的命令即可从上次中断的地方继续
- 脚本会自动跳过已下载的物种

## 查看进度

**查看已完成数量:**
```bash
python -c "import json; print(json.load(open('/GenSIvePFS/users/data/pypath_cache/string_download_progress.json'))['total_completed'])"
```

**查看失败数量:**
```bash
python -c "import json; print(len(json.load(open('/GenSIvePFS/users/data/pypath_cache/string_download_failed.json'))))"
```

## 示例输出

```
================================================================================
Starting STRING data download
================================================================================
Fetching list of all species...
Total species: 12533
Already completed: 250
Pending downloads: 12283
Using 8 parallel workers
Max retries per species: 3

Downloading STRING: 45%|████████▌         | 5500/12283 [1:23:45<2:15:30, 1.2s/species]
  success: 5450, failed: 50, rate: 65.2 species/min
```

## 注意事项

1. **磁盘空间**: 确保 `/GenSIvePFS/users/data/pypath_cache` 有足够空间（约 50-100GB）
2. **网络稳定**: 建议在网络稳定的环境下运行
3. **并行数**: 不要设置过高的 `--workers`，可能导致网络拥塞或被限流
4. **运行时间**: 完整下载所有物种可能需要数小时到数天（取决于网络速度）

## 故障排查

**如果下载一直失败:**
1. 检查网络连接
2. 查看日志文件了解具体错误
3. 尝试减少并行线程数: `--workers 4`
4. 增加重试次数: `--max-retries 5`

**如果进度文件损坏:**
```bash
rm /GenSIvePFS/users/data/pypath_cache/string_download_progress.json
rm /GenSIvePFS/users/data/pypath_cache/string_download_failed.json
```
然后重新运行脚本。
