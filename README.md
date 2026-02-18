# Lifted-Product Quantum LDPC Code Construction and Simulation

MATLAB toolkit for constructing and simulating Lifted-Product Quantum Low-Density Parity-Check (LP-QLDPC) codes with min-sum (belief propagation) decoding. Includes protograph-based classical LDPC code construction, CSS code lifting via the algebra homomorphism, symplectic Gaussian elimination, and Monte Carlo simulation of decoding performance under depolarizing noise.

**Key Features:**
- Protograph-based code construction via circulant lifting  
- CSS stabilizer matrix generation from classical seed codes  
- Binary matrix row-reduction with column pivoting (Gaussian elimination mod 2)  
- Min-sum decoder (message-passing algorithm) for syndrome decoding  
- Logical operator extraction via symplectic LQ decomposition  
- Performance simulation: plots error probability vs. physical error rate

---

## Contents

| File | Purpose |
|------|---------|
| `LPcodeConstructDecode.m` | Main script: constructs [[1054,140,20]] LP-QLDPC code, simulates decoding |
| `syndrome_MSA_seq_vars_5.m` | Min-sum algorithm (belief propagation) decoder |
| `ArrayWithHandle.m` | Handle class wrapper for large sparse matrices (pass-by-reference) |
| `read_L.m` | Reads logical operator matrices from `.alist` format |
| `LP_Tanner155_lx.alist` | Precomputed logical X-operators for the [[1054,140,20]] code |

**Image assets:**
- `QCLDPC__3_5__PCmatrix.png` — visualizes the (3,5) protograph base matrix **B**
- `QCLDPC__3_5__GeneratorMatrix.png` — generator matrix G for the classical [155,64,20] QC-LDPC code
- `__1054_140_20__LPQLDPC_performancecurve_critical_p.png` — error threshold plot for the quantum code
- `_155_64_20_QCLDPC_performancecurve.png` — performance curve for the classical precursor code

---

## Quick Start

```matlab
% Run the main script to construct and simulate the [[1054,140,20]] LP-QLDPC code:
LPcodeConstructDecode

% This will:
% 1. Build a [155,64,20] QC-LDPC code from a (3,5) protograph (lift size l=31)
% 2. Lift the protograph into a [[1054,140,20]] CSS quantum code
% 3. Simulate min-sum decoding performance vs. Z-flip probability
% 4. Plot the decoding error probability curve
```

**Expected output:** A plot showing decoding error probability vs. physical Z-error rate, with threshold near p ≈ 0.06.

---

## Code Construction Pipeline

### Step 1: Classical QC-LDPC Code from Protograph

The (3,5) **base matrix** (protograph) is:

<img width="280" height="95" alt="image" src="https://github.com/user-attachments/assets/a1784544-2e8a-4cde-8864-93fa6c9d5623" />

Each entry `x^k` represents a circulant permutation matrix of size l × l (here l = 31). The full parity-check matrix **H** is obtained by replacing each entry with its corresponding circulant:

```matlab
l = 31;
x = circshift(eye(l), 1);   % Cyclic permutation
protoBpowers = [1 2 4 8 16; 5 10 20 9 18; 25 19 7 14 28];

% Lift protograph to get full H matrix
HB.data = sparse(cell2mat(arrayfun(@(int) x^int, protoBpowers, 'UniformOutput', false)));
```

The resulting classical code is a **[155, 64, 20]** QC-LDPC code (155 bits, 64 information bits, minimum distance 20). The `reducePCmat()` function performs Gaussian elimination to extract the generator matrix **G** in systematic form [P | I].

### Step 2: CSS Quantum Code via Lifted-Product

The **Lifted-Product (LP)** construction embeds the classical protograph into a quantum stabilizer code. Given base matrix **B**, define:

**B_X** = [B ⊗ I_n | I_m ⊗ (l−B)^T]  
**B_Z** = [I_n ⊗ B | (l−B)^T ⊗ I_m]

where m = 3 (rows), n = 5 (columns). Each block is lifted to size l × l. The CSS parity-check matrix is:

**H_CSS** = [**H_X**, 0; 0, **H_Z**]

encoded in symplectic form (stacking X-type and Z-type stabilizers). The resulting code is **[[1054, 140, 20]]**:
- 1054 physical qubits  
- 140 logical qubits  
- Minimum distance d = 20

```matlab
mb = height(protoBpowers);  nb = width(protoBpowers);
BX = [kron(protoBpowers, eye(nb)), kron(eye(mb), (l - protoBpowers)')];
BZ = [kron(eye(nb), protoBpowers), kron((l - protoBpowers)', eye(mb))];

% Lift to full CSS matrix (X and Z stabilizers stacked)
HB.data = sparse(cell2mat(arrayfun(@(int) (int==0)*zeros(l) + (int~=0)*x^int, ...
    [BX, zeros(nb*mb, nb^2+mb^2); zeros(nb*mb, nb^2+mb^2), BZ], 'UniformOutput', false)));
```

### Step 3: Logical Operators

Logical operators are extracted in two ways:

1. **Computed via symplectic Gaussian elimination** (on-the-fly, see commented-out section in script)
2. **Loaded from precomputed file** (`LP_Tanner155_lx.alist`) — the logical X-operators are provided by Narayanan Rengaswamy's research group and satisfy the commutation relations with Z-type stabilizers.

The script verifies that the loaded logical operators commute correctly with the stabilizer matrix:

```matlab
LX = read_L('LP_Tanner155_lx.alist');
fprintf('Logicals are inconsistent with Z-type stabilizers: %d\n', ...
    full(any(mod(LX * HB.data(446:end, 1055:end).', 2), 'all')));
```

---

## Function Reference

### `LPcodeConstructDecode.m` (Main Script)

**Sections:**

1. **Classical (3,5)-QC-LDPC code** — constructs the [155,64,20] code, verifies parity-check conditions, optionally simulates BSC decoding
2. **LP-QLDPC code generation** — lifts the protograph to [[1054,140,20]], computes CSS stabilizer matrix
3. **Logical operator extraction** — row-reduces augmented system to find logical X-operators (or loads from `.alist` file)
4. **Error performance simulation** — Monte Carlo trials of Z-error correction, plots error probability vs. physical error rate

**Key helper functions (defined inline):**

#### `reducePCmat()`

Binary Gaussian elimination with column pivoting to compute the row-reduced echelon form (RREF) of the global parity-check matrix `HB.data`. Returns the generator matrix in systematic form.

**Algorithm:**
1. Lower triangulation: for each pivot row r, eliminate 1s below the diagonal via row additions (mod 2)
   - If diagonal entry is 0, perform row swap (or column swap if entire lower column is 0)
2. Upper triangulation: eliminate 1s above each pivot
3. Extract kernel basis: columns corresponding to free variables form the generator matrix

**Outputs:**
- `genmat` — generator matrix G (sparse, systematic form)
- `rowOPSmatrix` — matrix encoding the sequence of row operations (useful for tracking transformations)
- `colOPSmatrix` — matrix encoding the column swaps (to restore original column ordering)

**Usage:**
```matlab
global HB HBheight HBwidth rowspan colspan;
HB = ArrayWithHandle();
HB.data = sparse(H_matrix);  % your parity-check matrix
HBheight = height(HB.data);  HBwidth = width(HB.data);
rowspan = 1:HBheight;  colspan = 1:HBwidth;

genmat = reducePCmat();  % performs in-place row-reduction on HB.data
fprintf('Number of information bits: %d\n', height(genmat.data));
```

**Important:** This function modifies `HB.data` in-place. To preserve the original matrix, save a copy before calling:
```matlab
HBcopy = HB.data;
genmat = reducePCmat();
HB.data = HBcopy;  % restore if needed
```

#### `simulateMinSum(p, b, nsamples)`

Monte Carlo simulation of min-sum decoding for the classical LDPC code under a biased binary symmetric channel (BSC).

**Inputs:**
- `p` — bit-flip probability (BSC crossover probability)
- `b` — channel bias (asymmetry: 0 → 1 occurs with prob p+b, 1 → 0 with prob p−b)
- `nsamples` — number of Monte Carlo trials

**Outputs:**
- `errorprob` — fraction of trials where decoding failed (returned codeword ≠ transmitted codeword)

**Workflow:**
1. Generate random information word, encode via G
2. Pass through BSC: flip each bit independently with prob p±b
3. Compute syndrome s = H·y (received word y)
4. Run min-sum decoder (see `syndrome_MSA_seq_vars_5`) to estimate error vector
5. Check if decoder output matches true error; count failures

**Example:**
```matlab
% Sweep over channel error rates
pvals = logspace(log10(0.05), log10(0.2), 10);
errorprobs = arrayfun(@(p) simulateMinSum(p, 0, 40), pvals);

% Plot performance curve
plot(pvals, errorprobs, '-o');
xlabel('BSC bit-flip probability (p)');
ylabel('Decoding error probability');
title('Error performance of [155,64,20] QC-LDPC code');
grid on;
```

#### `SimulateMinSum_Zchannel(p, nsamples)`

Monte Carlo simulation of min-sum decoding for the **quantum** LP-QLDPC code under a Z-type (phase-flip) depolarizing channel.

**Inputs:**
- `p` — single-qubit Z-error probability
- `nsamples` — number of Monte Carlo trials

**Outputs:**
- `errorprob` — fraction of trials where decoding resulted in a logical error (non-trivial logical Z operator)

**Workflow:**
1. Sample a physical Z-error pattern from i.i.d. Bernoulli(p) on all qubits
2. Compute syndrome s = H_X · e (X-type stabilizers measure Z-errors)
3. Run min-sum decoder to estimate error vector ê
4. Check if residual error e + ê is a logical operator: (e + ê) · LX^T ≠ 0 (mod 2)
5. Count logical failures

**Example:**
```matlab
% Sweep over physical error rates
pvals = logspace(log10(0.04), log10(0.09), 6);
errorprobs = arrayfun(@(p) SimulateMinSum_Zchannel(p, 45), pvals);

% Plot threshold curve
plot(pvals, errorprobs, '-o', 'MarkerFaceColor', 'blue');
xlabel('Z-flip probability');
ylabel('Decoding error probability');
title('Error performance of [[1054,140,20]] LP-QLDPC code against Z-errors');
grid on;
```

#### `BSCwithbias(x, p, b)`

Single-bit channel simulator for the biased BSC.

**Inputs:**
- `x` — input bit (0 or 1)
- `p` — base flip probability
- `b` — bias parameter

**Outputs:**
- `y` — output bit (flipped with prob p + (−1)^x · b)

**Logic:** With probability 1 − p − (−1)^x·b, output x unchanged; otherwise flip to 1−x.

---

### `syndrome_MSA_seq_vars_5.m` — Min-Sum Decoder

Iterative belief propagation (BP) decoder using the **min-sum approximation** to the sum-product algorithm. Decodes syndromes for classical LDPC codes or the X/Z sectors of CSS quantum codes.

**Inputs:**
- `H` — m × n parity-check matrix (sparse)
- `vn_neighbors` — cell array: `vn_neighbors{j}` = indices of check nodes adjacent to variable node j
- `cn_neighbors` — cell array: `cn_neighbors{i}` = indices of variable nodes adjacent to check node i
- `syndrome` — 1 × m binary syndrome vector
- `llr` — log-likelihood ratio log((1−p)/p) for the channel (prior belief)
- `max_iter` — maximum number of BP iterations

**Outputs:**
- `decoder_output` — 1 × n estimated error vector
- `Calculated_syndrome` — syndrome of decoder output (for verification)
- `decoding_success` — 1 if syndrome matches, 0 otherwise

**Algorithm:**

1. **Initialization:** Set all variable node beliefs to channel LLR; send to check nodes
2. **Check node update (min-sum rule):**
   - For check node i: compute extrinsic LLRs by taking the product of signs and minimum magnitude of incoming messages (excluding the outgoing edge)
   - Scale by syndrome parity: flip sign if syndrome bit = 1
   - Apply damping factor λ = 0.8 to stabilize convergence
3. **Variable node update:**
   - For variable node j: sum all incoming check node messages + channel LLR
   - Propagate updated belief to check nodes
4. **Hard decision:** threshold total LLR at 0 → estimate error bit
5. **Stopping criterion:** recompute syndrome; if it matches input, halt early (successful decode)

**Usage:**
```matlab
% Precompute neighbor lists
cn_neighbors = cellfun(@(row) find(row > 0), mat2cell(H.', n, ones(1,m)), 'UniformOutput', false);
vn_neighbors = cellfun(@(col) find(col > 0), mat2cell(H, m, ones(1,n)), 'UniformOutput', false);

% Decode
llr = log2((1 - p) / p);  % channel LLR for error probability p
[error_estimate, syndrome_check, success] = syndrome_MSA_seq_vars_5(H, vn_neighbors, cn_neighbors, syndrome, llr, 40);

if success
    disp('Decoding converged to valid codeword');
else
    disp('Decoding failed or reached max iterations');
end
```

**Performance note:** The nested loop structure (`for j = 1:n; for i = 1:length(col1s)`) is deliberately sequential to update messages immediately (flooding schedule variant). For large codes (n > 10^4), consider vectorization or C-mex acceleration.

---

### `read_L.m` — Logical Operator Import

Reads a logical operator matrix from the `.alist` format (David MacKay's sparse matrix representation, commonly used in coding theory).

**`.alist` Format:**
```
k  n                       % dimensions: k rows, n columns
row_wt_max  col_wt_max     % max row/column weights
<row weights>              % k integers (weight of each row)
<column weights>           % n integers (weight of each column)
<row support lists>        % k lines: 1-indexed column indices of 1s in each row
```

**Function signature:**
```matlab
L = read_L(filename)
```

**Output:**
- `L` — k × n sparse binary matrix (1s at specified positions)

**Example:**
```matlab
LX = read_L('LP_Tanner155_lx.alist');
disp(size(LX));  % should print [140 1054] for the [[1054,140,20]] code
```

---

### `ArrayWithHandle.m` — Sparse Matrix Wrapper

MATLAB handle class to wrap large sparse matrices. Unlike value classes, handle objects are **passed by reference**, which is essential for in-place modification of global matrices in helper functions like `reducePCmat()`.

**Usage:**
```matlab
global HB;
HB = ArrayWithHandle();
HB.data = sparse(my_matrix);  % assign matrix

% Pass to function that modifies HB.data in-place
reducePCmat();  % no explicit return value needed; HB.data is modified globally
```

**Why handle class?** MATLAB's default behavior copies arrays when passed to functions. For matrices with millions of entries (e.g., 1054 × 1054 CSS stabilizer matrices), copying is prohibitively expensive. Handle classes sidestep this by maintaining a single shared instance.

---

## Simulation Workflow & Reproducing Plots

### Reproducing the Classical [155,64,20] QC-LDPC Performance Curve

```matlab
% Section 1: QC-LDPC construction
l = 31;
x = circshift(eye(l), 1);
protoBpowers = [1 2 4 8 16; 5 10 20 9 18; 25 19 7 14 28];
global HB HBcopy HBheight HBwidth rowspan colspan genmat cn_neighbors vn_neighbors;
HB = ArrayWithHandle();
HB.data = sparse(cell2mat(arrayfun(@(int) x^int, protoBpowers, 'UniformOutput', false)));
HBcopy = HB.data;
HBheight = height(HB.data);  HBwidth = width(HB.data);
rowspan = 1:HBheight;  colspan = 1:HBwidth;
cn_neighbors = cellfun(@(row) colspan(row > 0), mat2cell(HB.data.', HBwidth, ones(1, HBheight)), 'UniformOutput', false);
vn_neighbors = cellfun(@(col) rowspan(col > 0), mat2cell(HB.data, HBheight, ones(1, HBwidth)), 'UniformOutput', false);
genmat = reducePCmat();

% Simulate decoding performance
pvals = logspace(log10(0.05), log10(0.2), 10);
errorprobs = arrayfun(@(p) simulateMinSum(p, 0, 40), pvals);

% Plot
figure;
plot(pvals, errorprobs, '-o', 'MarkerFaceColor', 'blue');
xlabel('BSC bit-flip probability (p)');
ylabel('Decoding error probability');
title('Performance curve of the [155,64,20] QC-LDPC code');
grid on;
```

**Expected result:** Matches `_155_64_20_QCLDPC_performancecurve.png` — threshold around p ≈ 0.11.

### Reproducing the Quantum [[1054,140,20]] LP-QLDPC Threshold Plot

Run the full script up through the quantum code construction sections, then:

```matlab
% Section: Error performance simulation (ensure HB.data = X-stabilizers only)
pvals = logspace(log10(0.04), log10(0.09), 6);
errorprobs = arrayfun(@(p) SimulateMinSum_Zchannel(p, 45), pvals);

% Plot
figure;
plot(pvals, errorprobs, '-o', 'MarkerFaceColor', 'blue');
xlabel('Z-flip probability');
ylabel('Decoding error probability');
title('Error performance of the [[1054,140,20]]-LP-QLDPC code against Z-errors');
grid on;
```

**Expected result:** Matches `!(__1054_140_20__LPQLDPC_performancecurve_critical_p.png)` — threshold near p ≈ 0.06, with <1% logical error rate at p = 0.06.

---

## Code Structure Visualization

The uploaded matrix images illustrate the sparse structure of the constructed codes:

| Image | Description |
|-------|-------------|
| `![QCLDPC__3_5__PCmatrix.png](QCLDPC__3_5__PCmatrix.png)` | **Parity-check matrix H** of the [155,64,20] QC-LDPC code. Block-circulant structure visible (3 × 5 blocks, each 31 × 31). Yellow = 1, blue = 0. |
| `![QCLDPC__3_5__PCmatrix_rowreduced.png](QCLDPC__3_5__PCmatrix_rowreduced.png)` | **Row-reduced form** of H after Gaussian elimination. Upper-right 64 × 64 block is identity (systematic form). |
| `![QCLDPC__3_5__GeneratorMatrix.png](QCLDPC__3_5__GeneratorMatrix.png)` | **Generator matrix G** extracted from RREF of H. First 91 columns are parity bits (matrix P), last 64 columns are identity (information bits). |

These visualizations confirm:
- The quasi-cyclic structure is preserved under lifting (protograph → full matrix)
- Row reduction successfully isolates the kernel (information subspace)
- The code has full rank (height of G = 64, matching the dimension)

---

## Dependencies

- **MATLAB R2020a or later** (uses `cellfun`, `cell2mat`, `sparse`, `eigs`, handle classes)
- **No toolboxes required** (all functions use base MATLAB)

**Performance notes:**
- Classical [155,64,20] code: simulation of 10 error rates × 40 trials ≈ 30 seconds
- Quantum [[1054,140,20]] code: simulation of 6 error rates × 45 trials ≈ 5 minutes (scales quadratically with n due to BP message passing)

For production-scale simulations (> 10^4 trials), consider:
- Parallelization via `parfor` (requires Parallel Computing Toolbox)
- Early stopping when confidence interval converges
- GPU-accelerated sparse matrix operations (requires custom kernels)

---

## Theoretical Background

### Lifted-Product Codes

The lifted-product construction (Panteleev & Kalachev, 2021) generates quantum LDPC codes with constant rate and linear distance by:

1. Starting with two classical LDPC codes (or a single protograph)
2. Forming the tensor product in a way that preserves the CSS structure
3. Lifting via a group algebra homomorphism (here, cyclic group Z_l)

**Key advantage:** Unlike surface codes (distance d ~ √n), LP-QLDPC codes achieve d ~ √n with nearly constant encoding rate k/n ≈ constant. This enables fault-tolerant quantum computation with asymptotically lower overhead.

**[[1054,140,20]] code parameters:**
- Rate k/n = 140/1054 ≈ 0.133
- Relative distance d/n = 20/1054 ≈ 0.019
- Stabilizer weight (row weight of H): ≈ 5 (sparse, enabling efficient syndrome extraction)

### Min-Sum Decoding Threshold

For a quantum code under independent X/Z errors with rate p, the **decoding threshold** is the supremum error rate p_th below which error correction succeeds with high probability as n → ∞. The [[1054,140,20]] code exhibits:

- **Min-sum threshold:** p_th ≈ 0.06 (empirically measured, see performance plot)
- **Comparison:** Surface codes have p_th ≈ 0.01 under similar depolarizing noise; LP-QLDPC codes offer a 6× improvement due to higher code distance and better connectivity

The threshold is *not* the minimum distance d/n because:
- Min-sum decoding (iterative BP) is suboptimal (maximum-likelihood decoding achieves higher threshold)
- Finite block length effects (n = 1054 is large but not asymptotic)

**Reference:** For optimal BP thresholds of LP-QLDPC codes, see:
- Panteleev & Kalachev, "Quantum LDPC Codes With Almost Linear Minimum Distance," *IEEE TIT* 2022
- Leverrier & Zémor, "Quantum Tanner Codes," *STOC* 2022

---

## Citation

If you use this code in your research, please cite the original LP-QLDPC construction:

```
@article{panteleev2021quantum,
  title={Quantum LDPC codes with almost linear minimum distance},
  author={Panteleev, Pavel and Kalachev, Gleb},
  journal={IEEE Transactions on Information Theory},
  volume={68},
  number={1},
  pages={213--229},
  year={2021},
  publisher={IEEE}
}
```

And the min-sum decoding analysis:

```
@article{roffe2020decoding,
  title={Decoding across the quantum LDPC code landscape},
  author={Roffe, Joschka and White, David R and Burton, Simon and Campbell, Earl},
  journal={Physical Review Research},
  volume={2},
  number={4},
  pages={043423},
  year={2020},
  publisher={APS}
}
```

---

## Acknowledgments

The [[1054,140,20]] LP-QLDPC code construction and logical operator files (`LP_Tanner155_lx.alist`) are based on the work of **Prof. Narayanan Rengaswamy** (University of Arizona, ECE). The min-sum decoder implementation (`syndrome_MSA_seq_vars_5.m`) is adapted from the sparse BP codebase developed in his Quantum Error Correction course (Fall 2025).
