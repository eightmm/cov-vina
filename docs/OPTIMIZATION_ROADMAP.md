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

#### 1. **Mixed-Ligand GPU Batching**
**Current:** Process 1 ligand at a time sequentially
**Proposed:** Process 8 ligands simultaneously on GPU

```python
# Current
for ligand in ligands:
    conformers = generate(ligand)  # 0.6s
    optimize(conformers)           # 0.1s
# Total: N * 0.7s

# Proposed
for batch in chunks(ligands, batch_size=8):
    conformers = generate_batch(batch)  # 0.8s for 8 ligands
    optimize_batch(conformers)          # 0.2s for 8 ligands
# Total: (N/8) * 1.0s → 5.6x speedup for conformers
```

**Difficulty:** High (requires pipeline refactoring)
**Impact:**
- 10 ligands: 7s → 2s (3.5x)
- 100 ligands: 70s → 13s (5.4x)

**Implementation:**
1. Refactor `generate_conformers_and_cluster()` to accept List[Mol]
2. Batch RMSD computation across ligands
3. Batch optimization across all poses from all ligands
4. Track pose → ligand mapping

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

### Current Performance
| Ligands | Current |
|---------|---------|
| 1       | 1.2s    |
| 10      | 3.0s    |
| 100     | 24s     |
| 1000    | 240s    |

### With All Optimizations
| Ligands | Optimized | Speedup |
|---------|-----------|---------|
| 1       | 0.5s      | 2.4x    |
| 10      | 1.5s      | 2.0x    |
| 100     | 6s        | 4.0x    |
| 1000    | 45s       | 5.3x    |

## Implementation Priority

### Phase 1: Quick Wins (1 week)
1. ✅ Pocket caching (DONE - 3.4x for batch)
2. GPU warmup pre-heating
3. Conformer caching for duplicates
4. Async file I/O

**Expected:** 2x overall speedup

### Phase 2: Major Refactor (2-3 weeks)
1. Mixed-ligand GPU batching
2. Batch RMSD computation
3. Batch optimization

**Expected:** Additional 2-3x speedup (4-6x total)

### Phase 3: Advanced (future)
1. Multi-GPU support
2. Distributed computing (Ray, Dask)
3. ML-based conformer prediction
4. Disk-based pocket feature cache

**Expected:** 10-100x for very large libraries

## Hardware Utilization

### Current
- **GPU:** 10-30% (underutilized due to small batches)
- **CPU:** 80-100% (conformer generation bottleneck)
- **Disk:** <5% (sequential writes)

### Target (after optimization)
- **GPU:** 60-80% (mixed-ligand batching)
- **CPU:** 60-80% (parallel conformer generation)
- **Disk:** 10-20% (async writes)

## Code Complexity Trade-offs

| Optimization | Code Complexity | Maintenance | Worth It? |
|--------------|-----------------|-------------|-----------|
| Pocket caching | Low | Low | ✅ Yes (DONE) |
| GPU warmup | Trivial | Trivial | ✅ Yes |
| Conformer cache | Low | Low | ✅ Yes |
| Async I/O | Medium | Medium | ✅ Yes |
| Mixed batching | **High** | **High** | ⚠️ Depends on use case |
| Multi-GPU | High | High | ❌ Not yet |

## When to Optimize What

**Single ligand users:**
- Priority: GPU warmup pre-heating
- Skip: Mixed-ligand batching

**Batch docking (10-100 ligands):**
- Priority: Conformer caching, async I/O
- Consider: Mixed-ligand batching

**High-throughput (1000+ ligands):**
- Priority: ALL optimizations
- Consider: Distributed computing

## References

- Current profiling: See `/tmp/profile_pipeline.py`
- Pocket caching PR: Commit c874c59
- Unified API PR: Commit 8cc861c
