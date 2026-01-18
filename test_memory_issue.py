import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent
sys.path.insert(0, str(project_root))

try:
    from pypath.share import settings
    if not hasattr(settings, 'context') and hasattr(settings, 'settings'):
        settings.context = settings.settings.context
except ImportError:
    pass

from biocypher import BioCypher

bc = BioCypher(
    biocypher_config_path=str(project_root / "config/biocypher_config.yaml"),
    schema_config_path=str(project_root / "config/schema_config.yaml"),
)

# Test mode
TEST_MODE = True
CACHE = True
export_as_csv = True
output_dir_path = "/GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG/biocypher-out"

print("=" * 80)
print("Testing InterPro adapter (the problematic one)...")
print("=" * 80)

try:
    from bccb.interpro_adapter import InterPro
    
    interpro_adapter = InterPro(
        organism=9606,
        test_mode=TEST_MODE
    )
    print("✓ InterPro adapter initialized")
    
    print("Attempting to download InterPro data (this may fail with OOM)...")
    interpro_adapter.download_interpro_data(cache=CACHE)
    print("✓ InterPro data downloaded successfully")
    
    if export_as_csv:
        interpro_adapter.export_as_csv(path=output_dir_path)
    
    bc.write_nodes(interpro_adapter.get_interpro_nodes())
    bc.write_edges(interpro_adapter.get_interpro_edges())
    print("✓ InterPro adapter completed successfully")
    
except MemoryError as e:
    print(f"✗ MEMORY ERROR: {e}")
    print("The InterPro adapter ran out of memory!")
    sys.exit(1)
except Exception as e:
    print(f"✗ ERROR: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("Memory test completed")
print("=" * 80)
