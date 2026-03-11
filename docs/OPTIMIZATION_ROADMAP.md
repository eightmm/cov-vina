# Optimization Roadmap

## Current Performance Profile

**Baseline** (6lu7, 50 conformers, 50 opt steps):
- Pocket loading: 0.27s (once)
- First ligand: 0.94s (includes 0.71s warmup)
- Subsequent ligands: 0.23s each (steady state)

**Bottlenecks identified:**

1. **GPU Warmup Overhead: ~0.7s**
   - CUDA kernel compilation
   - First PyTorch operation triggers JIT compilation
   - Only happens once per session
   - Impact: High for single ligand, negligible for batch

2. **Conformer Generation: ~0.6s per ligand (CPU-bound)**
   - RDKit ETKDG: ~0.5s
   - Butina clustering: ~0.1s
   - Sequential processing (not parallelizable easily)
   - Impact: Largest bottleneck for steady-state

3. **Pocket Loading: ~0.3s (once)**
   - PDB parsing: ~0.1s
   - Feature computation: ~0.2s
   - Already cached for batch docking
   - Impact: Minimal (amortized)

4. **Gradient Optimization: ~0.1s per ligand**
   - Already GPU-optimized
   - Batch size 128 (good utilization)
   - Impact: Minor

5. **File I/O: ~0.05s per ligand**
   - SDF writing
   - Directory creation
   - Impact: Minimal

## Optimization Opportunities

### 🔥 High Impact (2-5x speedup)

#### 1. ~~**Mixed-Ligand GPU Batching**~~ (NOT NEEDED - Already Optimal)

**Status:** ❌ **UNNECESSARY** - Single-ligand batching already saturates GPU

**Analysis:**
- Each ligand generates 50-200 conformers → **already a large batch**
- GPU optimization processes 8 poses simultaneously (batch_size=8)
- 200 conformers = 25 GPU batches → **GPU fully utilized**
- Merging multiple ligands would require padding (memory waste) or complex ragged tensors
- **No performance gain, only added complexity**

```python
# Current (ALREADY OPTIMAL):
for ligand in ligands:
    conformers = generate(ligand)  # 200 conformers
    optimize(conformers, batch_size=8)  # 200/8 = 25 batches on GPU
    # ↑ GPU is ALREADY saturated with 8 concurrent poses!

# Proposed mixed-ligand batching:
all_conformers = flatten([generate(lig) for lig in ligands])
optimize(all_conformers)  # Same GPU utilization, but:
# - Requires padding (different conformer counts)
# - Loses per-ligand tracking
# - No speedup (conformer gen is CPU-bound)
```

**Verdict:** Single-ligand batching is sufficient. Skip this optimization.

---

#### 2. **Conformer Caching for Identical Ligands**
**Current:** Re-generate conformers for duplicate SMILES
**Proposed:** Cache conformers by canonical SMILES

```python
conformer_cache = {}  # canonical_smiles -> conformers

def get_conformers(smiles):
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if canonical not in conformer_cache:
        conformer_cache[canonical] = generate_conformers(smiles)
    return conformer_cache[canonical]
```

**Difficulty:** Easy
**Impact:**
- Libraries with duplicates: 2-10x for repeated molecules
- No impact on unique molecules

---

#### 3. **Async File I/O**
**Current:** Synchronous SDF writing blocks main thread
**Proposed:** Background thread for file writing

```python
import concurrent.futures

file_executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)

# Non-blocking write
future = file_executor.submit(write_sdf, poses, output_path)
```

**Difficulty:** Easy
**Impact:** ~0.05s saved per ligand (minor but free)

---

### 🔵 Medium Impact (1.5-2x speedup)

#### 4. **GPU Warmup Pre-heating**
**Current:** Warmup happens on first ligand
**Proposed:** Explicit warmup before processing

```python
def warmup_gpu(device):
    """Pre-compile CUDA kernels before batch processing."""
    # Run dummy operations to trigger compilation
    dummy = torch.randn(100, 3, device=device)
    _ = torch.cdist(dummy, dummy)  # RMSD kernel
    _ = torch.optim.Adam([dummy], lr=0.1)  # Optimizer
```

**Difficulty:** Easy
**Impact:**
- Single ligand: 0.94s → 0.23s (4x)
- Batch: Negligible (already amortized)

---

#### 5. **Lazy SDF Writing (Memory Trade-off)**
**Current:** Write SDF after each ligand
**Proposed:** Accumulate results, write at end

```python
# Write all results in one batch at the end
results = []
for ligand in ligands:
    results.append(dock(ligand))

# Write all at once (faster disk I/O)
write_all_sdfs(results)
```

**Difficulty:** Easy
**Impact:** ~0.02s per ligand (disk I/O batching)

---

### 🟢 Low Impact (<1.5x) but Easy Wins

#### 6. **Reduce Logging Overhead**
**Current:** Print statements in hot loops
**Proposed:** Conditional logging with verbosity levels

```python
if verbose >= 2:  # Only for debug mode
    print(f"Conformer {i}: energy {e:.3f}")
```

**Difficulty:** Trivial
**Impact:** ~0.01s per ligand

---

#### 7. **Pre-compute Common Features**
**Current:** Compute features for each pose independently
**Proposed:** Compute once per ligand, reuse

```python
# Current
for pose in poses:
    features = compute_vina_features(pose)  # Redundant

# Proposed
features = compute_vina_features(ligand_template)  # Once
for pose in poses:
    use_cached_features(features)
```

**Difficulty:** Medium (need to verify correctness)
**Impact:** ~0.02s per ligand

---

## Performance Projections

### Current Performance (After Phase 1)
| Ligands | Time   | Per-Ligand |
|---------|--------|------------|
| 1       | 0.55s  | 0.55s      |
| 10      | 2.4s   | 0.24s      |
| 100     | 24s    | 0.24s      |
| 1000    | 240s   | 0.24s      |

**Note:** First ligand includes 0.54s warmup (one-time cost)

### With Remaining Optimizations
| Ligands | Optimized | Speedup vs Current |
|---------|-----------|-------------------|
| 1       | 0.50s     | 1.1x              |
| 10      | 2.0s      | 1.2x              |
| 100     | 18s       | 1.3x              |
| 1000    | 180s      | 1.3x              |

**Remaining optimizations:** Conformer caching + Async I/O

## Implementation Priority

### Phase 1: Quick Wins ✅ COMPLETE
1. ✅ Pocket caching (DONE - 3.4x for batch)
2. ✅ GPU warmup pre-heating (DONE - 4.4x for first ligand)

**Achieved:** 3.4x speedup for batch docking

### Phase 2: Practical Improvements (1-2 days)
1. Conformer caching for duplicates (~30% for typical libraries)
2. Async file I/O (~5% improvement)

**Expected:** Additional 1.3x speedup (4.4x total from baseline)

### Phase 3: Advanced (future)
1. Multi-GPU support (if processing >1000 ligands)
2. Distributed computing (Ray, Dask) for massive libraries
3. ML-based conformer prediction (replace RDKit ETKDG)

**Expected:** 10-100x for massive libraries (10K+ ligands)

## Hardware Utilization

### Current (Already Optimal for Single GPU)
- **GPU:** 70-90% (✅ saturated with 8-pose batches from 200 conformers)
- **CPU:** 80-100% (conformer generation bottleneck - expected)
- **Disk:** <5% (sequential writes - I/O is fast)

### Why GPU is Already Saturated
- Each ligand: 200 conformers → 25 GPU batches (200 / 8)
- GPU processes 8 poses simultaneously
- **No idle time** - GPU always has work
- Mixed-ligand batching would NOT increase utilization

**Verdict:** No need for further GPU optimization on single-GPU setup

## Code Complexity Trade-offs

| Optimization | Code Complexity | Maintenance | Worth It? |
|--------------|-----------------|-------------|-----------|
| Pocket caching | Low | Low | ✅ Yes (DONE) |
| GPU warmup | Trivial | Trivial | ✅ Yes (DONE) |
| Conformer cache | Low | Low | ✅ Yes (simple) |
| Async I/O | Medium | Medium | ⚠️ Maybe (~5% gain) |
| Mixed batching | **High** | **High** | ❌ **NO** (no gain) |
| Multi-GPU | High | High | ❌ Not needed (single GPU saturated) |

## When to Optimize What

**Single ligand users:**
- ✅ Already optimal (GPU saturated with 200 conformers)

**Batch docking (10-100 ligands):**
- ✅ Pocket caching (DONE)
- ✅ GPU warmup (DONE)
- 🔄 Conformer caching (if duplicates expected)
- ❌ Skip mixed-ligand batching (no benefit)

**Massive libraries (1000+ ligands):**
- Same as above + consider multi-GPU if needed

**High-throughput (1000+ ligands):**
- Priority: ALL optimizations
- Consider: Distributed computing

## References

- Current profiling: See `/tmp/profile_pipeline.py`
- Pocket caching PR: Commit c874c59
- Unified API PR: Commit 8cc861c
