import numpy as np
import sys

filename = sys.argv[1]
data = np.load(filename)

print(f"✓ Datei: {filename}")
print(f"✓ Shape: {data.shape}")
print(f"✓ Dtype: {data.dtype}")
print(f"✓ Min: {data.min():.4f}")
print(f"✓ Max: {data.max():.4f}")
print(f"✓ Mean: {data.mean():.4f}")
print(f"\n✓ Erste Werte:")
print(data.flat[:10])  # Erste 10 Werte
print(f"\n✓ Enthält NaN: {np.isnan(data).any()}")
print(f"✓ Enthält Inf: {np.isinf(data).any()}")

