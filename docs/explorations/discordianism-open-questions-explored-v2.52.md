# Discordianism: Open Questions Explored (v2.52): the joining arrow freezes the vertex block

*This round starts ChatGPT's Program 2, the joined-ends configuration Q₀ = 1⇉2⇇3 + 1→3. On the unjoined 1⇉2⇇3,
Theorem 3's modulus is the vertex block of a fold: the coefficients of vertices in P(e_v) (v2.27). Here I prove that on
Q₀ in characteristic 0 the vertex block takes exactly 21 values, all 0/1 matrices, and that all 21 occur. The
joining arrow does not absorb Theorem 3's direction into a symmetry orbit; it forbids it. The proof is short:
necessary conditions from the kernel of Δ and from three arrow Δ-conditions, eliminated over ℚ, plus 21 integer folds
checked in Lean.*

---

## 0. Results

The predictions are in logs/v2_52/PREDICTIONS.md; 47 and 48 were added after 45 and 46 had held.

| # | Statement | Status |
|---|---|---|
| 45a | On Q₀ the vertex blocks of all census folds form one finite set of integer matrices, the same at q = 3, 5, 7 | held: 21 blocks at each prime, identical sets |
| 45b | On 1⇉2⇇3 the number of distinct vertex blocks grows with q | held: 31, 39, 47, 63 at q = 3, 5, 7, 11, which is 4q + 19 |
| 46 | On Q₀ every vertex block has entries in {−1, 0, 1} | held, and more: the entries are in {0, 1} |
| 47 | P_K-closure plus the Δ-condition of c alone cut the vertex blocks to the 21 | **failed**: that system is 1-dimensional. Adding the Δ-conditions of a and b gives exactly the 21 points (below). |
| 48 | Without the arrow c the projection is positive-dimensional | held: 2-dimensional with P_K-closure alone, 1-dimensional with the a and b conditions added |

## 1. The theorem

**Setting.** Q₀ = 1⇉2⇇3 + 1→3, with arrows a₁, a₂: 1→2, b₁, b₂: 3→2, c: 1→3, and the paths cb₁, cb₂. Arrow deletion
acts as Δ(a_i) = e₁⊗e₂, Δ(b_j) = e₃⊗e₂, Δ(c) = e₁⊗e₃, Δ(cb_j) = e₁⊗b_j + c⊗e₂.

For a fold P, let its *vertex block* be V[i][j] = coefficient of e_i in P(e_j).

**Theorem (vertex rigidity on joined ends).** Let k be a field of characteristic 0. The vertex block of every fold of
kQ₀ is one of the following 21 matrices, and each of them is the vertex block of an integer fold, which is a fold over
every commutative ring. Below, rows are i and columns are j, so column j is the vertex part of P(e_j).

```
000 000 000 | 000 000 001 | 000 001 001 | 000 010 000 | 000 010 001 | 000 010 010 | 001 000 001
001 001 001 | 001 010 001 | 010 010 000 | 010 010 010 | 100 000 000 | 100 000 001 | 100 000 100
100 001 001 | 100 010 000 | 100 010 001 | 100 010 010 | 100 010 100 | 100 100 000 | 100 100 100
```

Each entry lists the three rows of V. For example, `100 010 100` means V = [[1,0,0],[0,1,0],[1,0,0]]: P(e₁) has vertex
part e₁ + e₃, P(e₂) has e₂, and P(e₃) has 0.

**Proof.**
1. K = ker Δ = span(e₁, e₂, e₃, a₁ − a₂, b₁ − b₂), over any field. ΔP = (P⊗P)Δ gives P(K) ⊆ K. Write P_K for the 5×5
   matrix of P on K. Its upper-left 3×3 block is V.
2. P_K satisfies:
   - P_K² = P_K;
   - for x, y ∈ K, P(Px·Py) = Px·Py and P((1−P)x·(1−P)y) = 0, because image and kernel are subalgebras and K is
     closed under products.
3. For each arrow x: s → t, Δ(P x) = (P⊗P)(e_s⊗e_t) = P(e_s) ⊗ P(e_t), where P(e_s) and P(e_t) are columns of P_K.
   Use this for x = a₂, b₂, c, with P(x) as 30 further unknowns.
4. Singular eliminates these 55 unknowns down to the 9 entries of V over ℚ (`analysis/v2_52_vertex_ideal.py`). The
   elimination ideal is zero-dimensional. Its associated primes (minAssGTZ) are 21 rational points, exactly the 21
   census blocks above. So V(E) over ℚ̄ is these 21 points, and every fold's V lies among them.
5. Realisation: for each block, a symmetric lift of a census fold over 𝔽₇ is a fold over ℤ, by an exact check
   (`analysis/v2_52_realise.py`). Lean then proves `vbk_isFold`, k = 1..21, by reflection, over every commutative ring
   (`Gen/pairs_in_short_vblocks.lean`, built in 1732 s). ∎

**What each condition contributes.**

| conditions | dimension of the vertex-block set |
|---|---|
| P_K only (this is also the unjoined 1⇉2⇇3: its K is the same algebra) | 2 |
| P_K + arrows a, b (the full arrow conditions of 1⇉2⇇3) | 1: Theorem 3's families |
| P_K + arrow c | 1 |
| P_K + arrows a, b, c | **0**: 21 points |

The mechanism, by hand: write P(e₁) = v + r₁ and P(e₃) = w + r₃, with vertex parts v, w and radical parts r₁, r₃ ∈ K.
Then Δ(P c) = P(e₁)⊗P(e₃) can only have the tensor terms that Δ produces: e₁⊗e₂, e₃⊗e₂, e₁⊗e₃, e₁⊗b_j and c⊗e₂. So:
- v wᵀ is supported on the positions (1,2), (3,2), (1,3);
- r₁ ⊗ r₃ = 0, so r₁ = 0 or r₃ = 0.

Theorem 3's families have v wᵀ with a nonzero (1,1) entry (for example V = [[−s,0,s],[0,0,0],[s,0,−s]]), so they are
excluded. v2.29 found the same in Lean as a failure "in the Δ condition of the arrow 1→3". The theorem shows that
nothing replaces them. On its own, though, the c-condition still leaves a 1-dimensional set; the a and b conditions are
needed too.

## 2. What this does and does not settle for Program 2

- **Settled.** Theorem 3's modulus is the one invariant known to grow on the unjoined quiver, and it is rigid on Q₀:
  21 values, independent of the field in characteristic 0. This is the precise sense in which the connecting arrow
  "kills" the Theorem 3 direction 1 + s(δ_b − δ_a). It forces the vertex block into a finite set, rather than having a
  symmetry absorb the direction.
- **Not settled.** Lemma (L) for Q₀ is the statement that every twist fibre is a finite union of census orbits. That
  needs the same rigidity for every other coordinate of P modulo the census group. The census is 366 local-torus
  orbits at q = 3, 5, 7, 11. A natural next step is to stratify by the 21 blocks and treat each stratum like V_A in
  v2.45–v2.50: the variety of folds with a given vertex block, its components, and the census-group orbits on them.

## 3. Next

1. Stratify the folds of Q₀ by the 21 vertex blocks: compute each stratum's variety (75 unknowns, Singular) and compare
   its components with the census orbits.
2. The same rigidity question for 1⇉2⇇3 + 1→3 + 2→4 (job 14's quiver): are its vertex blocks also finitely many?

## Files

- analysis/v2_52_vertex_blocks.py, v2_52_vertex_ideal.py, v2_52_realise.py
- lean/gen/gen_vblocks.py → lean/EatTheFibre/Gen/pairs_in_short_vblocks.lean
- logs/v2_52/: PREDICTIONS.md, vertex_blocks_{pairs_in_short,pairs_in}.txt, vertex_ideal_{short,short_abc,in,in_ab}.txt,
  realise_7.json, lean_vblocks.log
