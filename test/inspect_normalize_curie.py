#!/usr/bin/env python3
"""
Inspect normalize_curie function to understand what it does.
"""

from bioregistry import normalize_curie
import inspect

print("="*60)
print("normalize_curie Function Inspection")
print("="*60)

# Print docstring
print("\nDocstring:")
print(normalize_curie.__doc__)

# Print signature
print("\nSignature:")
print(inspect.signature(normalize_curie))

# Test with examples
print("\n" + "="*60)
print("Test Examples")
print("="*60)

test_cases = [
    "interpro:IPR000001",
    "uniprot:P12345",
    "go:0008150",
    "INTERPRO:IPR000001",  # uppercase prefix
    "InterPro:IPR000001",  # mixed case
    "invalid:test",
]

for test in test_cases:
    result = normalize_curie(test)
    print(f"Input:  {test:30s} -> Output: {result}")

print("\n" + "="*60)
print("Comparison: normalize_curie vs simple concatenation")
print("="*60)

# Compare with simple concatenation
prefixes = ["interpro", "uniprot", "go"]
identifiers = ["IPR000001", "P12345", "0008150"]

for prefix, identifier in zip(prefixes, identifiers):
    simple = prefix + ":" + identifier
    normalized = normalize_curie(prefix + ":" + identifier)
    match = "✓" if simple == normalized else "✗"
    print(f"{match} {prefix}:{identifier}")
    print(f"  Simple:     {simple}")
    print(f"  Normalized: {normalized}")
    if simple != normalized:
        print(f"  DIFFERENCE DETECTED!")
    print()
