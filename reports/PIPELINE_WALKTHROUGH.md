# CovVina Pipeline 작동 방식

**실제 실행 예시**: SARS-CoV-2 Mpro (CYS145)에 acrylamide 워헤드 docking

```bash
uv run python scripts/run_covalent_pipeline.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o output/ \
  --optimize \
  --num_confs 1000 \
  --steps 200
```

---

## 파이프라인 단계별 실행 흐름

### 📥 **입력**

```
Protein PDB: examples/6lu7/6lu7_pocket.pdb
Query SMILES: C=CC(=O)N[C@@H](C)C(=O)O  (vinyl-alanine)
Reactive Residue: CYS145
```

---

### **Step 1: Pocket 추출 & Anchor Point 계산**

**파일**: `src/cov_vina/io/pocket.py:load_pocket_bundle()`

```python
# 1.1. PDB 파일 읽기
protein_mol = Chem.MolFromPDBFile("6lu7_pocket.pdb")

# 1.2. Reactive residue 찾기 (CYS145)
anchors = find_reactive_residues(protein_mol, residue_spec="CYS145")
# → anchor.residue_name = "CYS"
# → anchor.residue_num = 145
# → anchor.atom_name = "SG"  (sulfur atom)
# → anchor.coord = [-13.985, 15.763, 66.922]  (3D 좌표)

# 1.3. CB (Cβ) 찾기
# → anchor.cb_coord = [-13.2, 16.1, 67.5]  (Cβ 좌표)
# → anchor.bond_vector = normalize(SG - CB) = [0.7, -0.5, -0.5]

# 1.4. Pocket 추출 (12Å cutoff)
pocket_atoms = [atom for atom in protein_mol if distance(atom, anchor.coord) < 12.0]
# → 101 atoms selected
```

**출력**:
```
Anchor: CYS145:A atom SG (bond length 1.82 Å)
Extracting pocket within 12.0Å of CYS145:A...
  Pocket: 101 atoms (full protein: 101 atoms)
```

---

### **Step 2: Ligand 처리**

**파일**: `src/cov_vina/io/input.py:process_query_ligand()`

```python
# 2.1. SMILES → RDKit Mol
query_mol = Chem.MolFromSmiles("C=CC(=O)N[C@@H](C)C(=O)O")
canonical_smiles = Chem.MolToSmiles(query_mol)
# → "C=CC(=O)N[C@@H](C)C(=O)O"  (10 heavy atoms)
```

**출력**:
```
Query ligand: C=CC(=O)N[C@@H](C)C(=O)O
```

---

### **Step 3: Warhead 탐지**

**파일**: `src/cov_vina/molecular/anchor.py:detect_warheads()`

```python
# 3.1. SMARTS 패턴 매칭 (18개 warhead type)
query_no_h = Chem.RemoveHs(query_mol)  # H 제거

for smarts, map_num, name in _WARHEAD_REGISTRY:
    pattern = Chem.MolFromSmarts("[CH2:1]=[CH]C(=O)[N,n]")  # acrylamide
    matches = query_no_h.GetSubstructMatches(pattern)
    # → [(0, 1, 2, 3, 4), ...]  # Atom indices

# 3.2. Reactive atom 추출
reactive_atom_idx = 0  # CH2 carbon (beta carbon)
warhead = WarheadHit(
    warhead_type="acrylamide",
    reactive_atom_idx=0,
    matched_atoms=(0, 1, 2, 3, 4)
)

# 3.3. 호환성 체크
is_compatible = check_warhead_residue_compatibility("acrylamide", "CYS")
# → (True, "acrylamide with CYS: well-established combination")
```

**출력**:
```
Warhead: acrylamide (reactive atom idx 0)
```

---

### **Step 4: Adduct Template 생성** ⭐

**파일**: `src/cov_vina/molecular/adduct.py:create_adduct_template()`

**핵심**: Conformer 생성 **전에** CB-S를 리간드에 추가!

```python
# 4.1. Leaving group 제거 (acrylamide는 Michael addition이라 없음)
leaving_group_pattern = LEAVING_GROUPS.get("acrylamide", [])
# → []  (no leaving group)

# 4.2. Michael addition 처리
# C=C-C(=O) 의 C=C bond를 C-C single bond로 변경
editable_mol.RemoveBond(0, 1)  # Remove C=C double bond
editable_mol.AddBond(0, 1, Chem.BondType.SINGLE)  # Add C-C single bond

# 4.3. Protein anchor atoms 추가
cb_idx = editable_mol.AddAtom(Chem.Atom(6))   # CB (carbon)  → index 10
s_idx = editable_mol.AddAtom(Chem.Atom(16))   # S (sulfur)   → index 11

# 4.4. Bond 형성
editable_mol.AddBond(10, 11, SINGLE)  # CB-S bond
editable_mol.AddBond(11, 0, SINGLE)   # S-C(reactive) bond

# Before: 10 atoms (ligand only)
# After:  12 atoms (ligand + CB + S)
```

**출력**:
```
Creating adduct template (removing leaving group, adding CB-S)...
  Original atoms: 10
  Adduct atoms: 12
  Added CB atom at index: 10
  Added SG atom at index: 11
  Reactive atom index: 0 → 0
```

**구조 변화**:
```
Before:  H2C=CH-C(=O)-NH-...
After:   CB-S-CH2-CH2-C(=O)-NH-...
         ↑   ↑
       protein anchor
```

---

### **Step 5: CB-S CoordMap 생성** ⭐

**파일**: `src/cov_vina/molecular/anchor.py:create_covalent_coordmap()`

**핵심**: CB와 S만 고정, reactive atom은 자유!

```python
coord_map = {
    10: Point3D(-13.2, 16.1, 67.5),    # CB 위치 (protein Cβ)
    11: Point3D(-13.985, 15.763, 66.922),  # S 위치 (protein Sγ)
    # NOTE: Atom 0 (reactive C)는 coordmap에 없음!
}
```

**출력**:
```
Generating conformers with CB-S fixed at anchor position...
  Fixed atoms in coordmap: [10, 11]
```

**왜 중요한가?**
- CB-S는 공간상 anchor (움직일 수 없음)
- Reactive atom은 자유 → ligand가 360° 회전 가능
- **Conformer diversity 최대화!**

---

### **Step 6: Conformer 생성 + Butina Clustering**

**파일**: `src/cov_vina/molecular/conformer.py:generate_conformers_and_cluster()`

```python
# 6.1. RDKit EmbedMultipleConfs (CB-S 고정)
AllChem.EmbedMultipleConfs(
    adduct_mol,
    numConfs=1000,
    coordMap=coord_map,  # CB와 S 위치 고정
    randomSeed=42
)
# → 1000 conformers 생성 (CB-S는 고정, 나머지는 다양한 pose)

# 6.2. GPU-accelerated RMSD 계산
coords = torch.tensor([conf.GetPositions() for conf in conformers])  # [1000, 12, 3]
rmsd_matrix = compute_pairwise_rmsd_gpu(coords)  # [1000, 1000]

# 6.3. Butina clustering
clusters = Butina.ClusterData(
    rmsd_matrix.cpu().numpy(),
    nPts=1000,
    distThresh=1.0  # 1Å RMSD threshold
)
# → 7 clusters found

representative_cids = [cluster[0] for cluster in clusters]
# → [40, 15, 18, 24, 21, 47, 1]  (cluster centroids)
```

**출력**:
```
Generating 1000 conformers...
Applying rigid constraints for 2 atoms during generation...
Calculating PyTorch batched RMSD matrix for 1000 conformers...
Clustering conformers with RMSD threshold 1.0Å...
Selected 7 representative conformers from 7 clusters.
Generated 7 representative conformers (Butina clustering with CB-SG)
```

**결과**:
- 1000 conformers → 7 representatives (RMSD > 1Å씩 떨어진 diverse poses)
- CB-SG 포함된 **전체 구조**로 RMSD 계산 (물리적으로 정확)

---

### **Step 7: MMFF94 Relaxation (Optional)**

**파일**: `src/cov_vina/molecular/relax.py:relax_pose_with_fixed_core()`

```python
# 7.1. CB-S 고정한 채로 MMFF force field로 relaxation
anchor_query_indices = {10, 11}  # CB and S

for cid in representative_cids:
    AllChem.MMFFOptimizeMoleculeConfs(
        adduct_mol,
        confId=cid,
        fixedAtoms=anchor_query_indices,
        maxIters=500
    )
```

**출력**:
```
Relaxing conformers via MMFF94 with CB-SG fixed...
  Conformer 40: applied: MMFF
  Conformer 15: applied: MMFF
  ...
```

**목적**: Steric clash 제거, geometry 정리

---

### **Step 8: Vina Scoring with Exclusion Mask**

**파일**: `src/cov_vina/scoring/vina_scoring.py`

```python
# 8.1. Exclusion mask 생성
ligand_exclude_indices = get_covalent_exclusion_indices(
    adduct_mol, warhead, n_hop_exclude=2
)
# → [0, 1, 10, 11]  (reactive C, neighbor, CB, S)

protein_exclude_residues = get_protein_exclusion_residues(anchor)
# → [CYS145]

# 8.2. Intermolecular exclusion mask
exclusion_mask = torch.zeros(12, 101, dtype=torch.bool)  # [ligand_atoms, pocket_atoms]
exclusion_mask[ligand_exclude_indices, :] = True  # CB-S region 제외
exclusion_mask[:, cys145_atoms] = True  # CYS145 제외

# 8.3. Vina scoring
scores = compute_vina_score(
    ligand_coords, pocket_coords,
    ligand_features, pocket_features,
    exclusion_mask=exclusion_mask
)
# → [0.125, -0.342, 0.089, ..., -0.217]  (7 poses)
```

**출력**:
```
Scoring with 'vina' weights (excluding covalent bond region)...
  Excluding 547/1212 atom pairs (45.1%)
```

**왜 exclusion?**
- CB-S-C bond 자체는 scoring 제외 (이미 covalent이라 intermolecular 아님)
- CYS145 자체도 제외 (double-counting 방지)

---

### **Step 9: Gradient-Based Torsion Optimization** ⭐

**파일**: `src/cov_vina/optimization/torsion.py:optimize_torsions_vina()`

**핵심**: CB 고정, S-C bond 포함 **모든 rotatable bonds** 최적화

```python
# 9.1. Torsional DOF 추출
ref_indices = [10]  # CB만 고정 (S는 회전 가능!)
torsion_dofs = extract_torsion_dofs(adduct_mol, ref_indices)
# → 5 rotatable bonds found (including S-C bond!)

# 9.2. Parameterization
coords = torsion_forward_kinematics(torsion_angles, torsion_dofs, ref_coords)
# coords = f(θ1, θ2, θ3, θ4, θ5)  (differentiable)

# 9.3. Gradient descent
optimizer = torch.optim.Adam([torsion_angles], lr=0.05)

for step in range(200):
    optimizer.zero_grad()

    # Forward
    coords = forward_kinematics(torsion_angles)
    score = compute_vina_score(coords, pocket_coords, ...)

    # Backward
    score.backward()
    optimizer.step()

    # Early stopping
    if converged:
        break
```

**출력**:
```
--- Gradient-Based Torsion Optimization (7 poses) ---
Optimizing 7 poses in 1 batches (batch_size=128)...
  Batch 1/1: Optimizing poses 0-6...
    Step 50: best=-0.325, avg=-0.187
    Step 100: best=-0.423, avg=-0.301
    Step 150: best=-0.468, avg=-0.389
    All 7 poses converged at step 173
✓ Optimization complete!
  Best pose: 6 score -0.475 kcal/mol (delta -0.550)
  Avg improvement: -0.585 kcal/mol
```

**결과**:
- Initial score: ~0.0 kcal/mol
- Optimized score: **-0.475 kcal/mol** (개선 -0.55)
- 173 steps에 수렴

---

### **Step 10: 결과 저장**

**파일**: `src/cov_vina/pipeline.py` (end)

```python
# 10.1. Score로 정렬
sorted_indices = torch.argsort(scores)  # [6, 3, 1, 4, 0, 2, 5]

# 10.2. SDF 파일로 저장
writer = Chem.SDWriter("output/covalent_poses_all.sdf")
for rank, idx in enumerate(sorted_indices):
    mol_copy = Chem.Mol(adduct_mol, confId=idx)
    mol_copy.SetProp("Rank", str(rank+1))
    mol_copy.SetProp("VinScore", f"{scores[idx]:.3f}")
    mol_copy.SetProp("Warhead", "acrylamide")
    writer.write(mol_copy)
writer.close()
```

**출력**:
```
Saving all 7 poses sorted by energy...

Top 5 energies:
Rank 1: Conformer 1, Energy: -0.475 kcal/mol
Rank 2: Conformer 24, Energy: -0.417 kcal/mol
Rank 3: Conformer 15, Energy: -0.409 kcal/mol
Rank 4: Conformer 40, Energy: -0.407 kcal/mol
Rank 5: Conformer 21, Energy: -0.345 kcal/mol
✓ Saved 7 poses to output/covalent_poses_all.sdf

-> Covalent docking completed in 2.01s
```

---

## 📊 성능 특성

| Stage | Time | GPU? | Notes |
|-------|------|------|-------|
| Pocket extraction | 0.01s | ❌ | RDKit (CPU) |
| Warhead detection | <0.01s | ❌ | SMARTS matching |
| Adduct template | <0.01s | ❌ | Topology edit |
| Conformer generation | 0.8s | ❌ | RDKit ETKDG |
| RMSD + Butina | 0.2s | ✅ | PyTorch GPU |
| MMFF relaxation | 0.1s | ❌ | RDKit force field |
| Vina scoring | 0.01s | ✅ | PyTorch GPU |
| Gradient optimization | 0.9s | ✅ | PyTorch GPU |
| **Total** | **2.0s** | - | **RTX PRO 6000** |

---

## 🔑 핵심 차별점

### 1. **Adduct-First Approach**
- ❌ **Old**: Conformer 생성 → 나중에 CB-S 추가 → Coordinate transfer 필요
- ✅ **New**: CB-S 추가 → Conformer 생성 → Clean & Simple

### 2. **CB-S CoordMap Strategy**
- ❌ **Old**: Reactive atom 고정 → Low diversity (단일 방향)
- ✅ **New**: CB-S만 고정 → High diversity (360° 회전)

### 3. **Physical Butina Clustering**
- ❌ **Old**: Ligand만으로 RMSD → CB-S 무시
- ✅ **New**: Adduct 전체로 RMSD → CB-S 포함 (정확)

### 4. **Flexible S-C Bond Optimization**
- ❌ **Old**: Reactive atom 고정 → S-C bond 회전 불가
- ✅ **New**: CB만 고정 → S-C bond 회전 가능 (자유도 증가)

---

## 📝 요약

```
Input: SMILES + Protein PDB + Residue
  ↓
[Step 1-3] 준비: Pocket 추출, Warhead 탐지, 호환성 체크
  ↓
[Step 4-5] ⭐ Adduct-first: CB-S 추가 → CB-S coordmap 생성
  ↓
[Step 6-7] Sampling: 1000 conformers → Butina clustering → MMFF relax
  ↓
[Step 8-9] ⭐ Optimization: Vina scoring → Gradient descent (CB 고정)
  ↓
Output: Ranked poses (SDF) with scores
```

**핵심 혁신**:
- CB-S를 **처음부터** 리간드에 포함 (Adduct-first)
- CB-S만 고정하고 reactive atom은 자유 (Maximum diversity)
- GPU-accelerated gradient optimization (Fast & Accurate)

**결과**: 2초 만에 physically accurate한 covalent complex 예측! 🚀
