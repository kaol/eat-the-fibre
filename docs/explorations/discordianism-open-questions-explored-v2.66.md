# Discordianism: Open Questions Explored (v2.66): the proof procedure on every small unobstructed quiver, and where the vertex blocks come from

*v2.54–v2.65 proved Conjecture 3ᶜ (finitely many census classes, with a bound uniform in q, in every characteristic)
for three quivers, one at a time. This round runs the whole procedure end to end on every small unobstructed quiver in
the census list, as a breadth and falsification test, and finds structure in its first step.
- **Conjecture 3ᶜ is proved for 12 more quivers** (§1): K₂, K₃ and the ten one-class quivers with 3 or 4 vertices.
  The thirteenth, 1⇛2⇇3 + 1→3, and the five-vertex 1⇉2→3→4→5 are still running (§1).
- **Duality (proved, char-free).** Path reversal turns folds of Q into folds of Q^op, compatibly with the census
  group. So every theorem here also holds for the opposite quivers (§2).
- **The collapse lemma (proved, char-free).** Collapsing every parallel class to one arrow is a retraction
  kQ ⇄ kQ^s of Δ-compatible algebras, and P ↦ ιPπ embeds Folds(kQ^s) into Folds(kQ), keeping the vertex block (§3).
- **Conjecture S (new, verified on 24 unobstructed quivers).** For unobstructed Q this embedding hits every vertex block: the
  vertex blocks of Q are exactly the folds of the simple quiver Q^s, one each. It fails, as it must, for the
  obstructed 1⇉2⇇3. A proof would replace the vertex step of the procedure, the only step that needed extra work at
  p = 2, by a statement about quivers without parallel arrows.*

---

## 0. Results

The predictions are in logs/v2_66/PREDICTIONS.md: 83–84 committed before any run, 85–86 and 87–90 added in two
later commits, each before its computation (one exception, see 85).

| # | Statement | Status |
|---|---|---|
| 83 | On each of the 13 quivers: zero-dimensional vertex ideal, every component passes, no exceptional prime, so 3ᶜ holds uniformly | **held on the 12 finished quivers**; 1⇛2⇇3 + 1→3 still running: no genuine failure so far, two pipeline steps needed a rerun (§1) |
| 84 | The vertex blocks are the census blocks at q = 3 | **held** on all 12 quivers with saved q = 3 representatives (10 from v2.60, pair_before and pair_mid computed here): equal as sets |
| 85 | Duality; pair_before has 71 classes at q = 3 | **held** (71). The duality is proved (§2). *Caveat:* the pair_before census had finished when the prediction was committed (same shell command), though I had not read it |
| 86 | Vertex blocks of pair_before, pair_mid = census blocks | **held** (23, 100) |
| 87 | The vertex ideal over ℤ does not depend on the class size m ≥ 2 | **held** on 11 pairs (K₂/K₃/K₄, and the m = 2, 3 versions of 7 three- and four-vertex quivers, plus Q₀/1⇛2⇇3 + 1→3): same leading terms of the strong Gröbner basis over ℤ, same vdim, same points |
| 88 | Control: with m = 1 the vertex ideal differs | **held** (A₂ vs K₂: vdim 6 vs 8; A₃ vs 1⇉2→3: 23 vs 27), **but the point sets are equal**, which led to 89 |
| 89 | Conjecture S on 7 unobstructed quivers | **held** on all 7, and afterwards on Q₁ and the 5-vertex path quiver too (§3) |
| 90 | S fails for 1⇉2⇇3 | **held**: its vertex variety is 1-dimensional, 1→2←3 has 27 |

## 1. The procedure, end to end (`analysis/v2_66_prove.py`)

For one quiver, the driver runs the scripts of v2.55–v2.59 in order:
1. **Integral structure** (`v2_58_integral.py`): Smith form of Δ, saturated K-basis, det [Kmat | e_C] = ±1, and the
   census group and weightings over ℤ.
2. **Vertex blocks over ℤ** (`v2_58_vertex_integral.py`): elimination, leading coefficients, radical certificate,
   point ideal. Where a leading coefficient is not ±1, the vertex ideal is decomposed directly over 𝔽_p for each prime
   p dividing it (`v2_66_vertex_badprimes.py`). This checks that every associated prime is a rational point and that
   the points are the reductions of the char 0 points, bijectively.
3. **The char 0 pipeline with certificates** on every vertex block (`v2_55_pipeline.py`; the integral case split if
   the default run fails), and the tally.
4. **Specialisation to every prime** (`v2_58_goodprimes.py`, p = 2 always checked directly).
5. **The Bézout bound** (`v2_58_bound.py` plus the stabiliser factor).

`analysis/v2_66_table.py` checks every condition from the summaries.

| quiver | vertex blocks | vertex step at p | components | specialisation claims | candidate primes (all good) | log₁₀ B | classes at q = 3 |
|---|---|---|---|---|---|---|---|
| K₂ = 1⇉2 | 6 | — | 14 | 76 | — | 7.5 | 20 |
| K₃ = 1⇛2 | 6 | — | 22 | 242 | — | 20.9 | 34 |
| 1⇉2→3 | 23 | — | 49 | 339 | — | 23.0 | 71 |
| 1→2⇉3 | 23 | — | 49 | 338 | — | 23.0 | 71 |
| 1⇉2, 1→3 | 27 | — | 67 | 402 | — | 11.8 | 101 |
| 1⇉2, 3→2 | 27 | — | 67 | 402 | — | 11.8 | 101 |
| 1⇛2→3 | 23 | — | 75 | 1073 | 2,3 | 53.6 | 120 |
| 1⇛2, 1→3 | 27 | — | 107 | 1287 | 2 | 29.4 | 183 |
| 1⇉2→3→4 | 100 | 2 | 239 | 1629 | — | 30.9 | 334 |
| 1→2⇉3→4 | 100 | 2 | 174 | 1194 | — | 41.4 | 257 |
| 1⇉2, 2→3, 2→4 | 148 | 2 | 263 | 2142 | — | 43.1 | 414 |
| 1⇉2, 1→3, 4→2 | 123 | 2 | 318 | 2071 | — | 15.3 | 485 |

Every row passes every check:
- a unimodular structure (Smith form 1, det B = ±1);
- a zero-dimensional vertex ideal whose radical is the ideal of its points (radical certificate);
- every component passing the v2.54 criterion, none unresolved;
- presolve pivots ±1 and no branch constants;
- no exceptional prime, and no unresolved cross-check.

"Vertex step at p" lists the primes dividing a non-unit leading coefficient. Only 2 occurs, in the four-vertex
quivers, and the direct decomposition over 𝔽₂ gives the same points. (The vertex ideals of quivers with parallel
arrows are not radical: vdim 116 against 100 points for 1⇉2→3→4, in characteristics 0 and 2 alike. Only the zero set
matters.) Point strata (blocks whose stratum is a single fold) carry no specialisation claim.

**Theorem.** For each quiver Q in the table, and for its opposite quiver (§2), the number of local-torus orbits of
Folds(𝔽_q Q) is at most B, for every prime power q.

Two remarks:
- The class counts at q = 3 come from the census, not from the proof. The proof bounds the count but does not
  compute it.
- Some quivers in the table also appear at more primes in earlier rounds, with constant counts (v2.60).

**Stop or continue.** This is the first time the procedure ran unattended over a list rather than being steered
quiver by quiver. Three gaps showed up, and each was closed generally rather than by hand:
1. **Non-unit leading coefficients at p = 2** (the four-vertex quivers). Closed by the direct char p check in step 2.
2. **No rational point on a component** (1⇛2⇇3 + 1→3, block 000,001,001). The v2.54 cut needs a point at which the
   orbit-dimension matrix has rank > r, only to choose which minor to add. On a component cut out by quadrics with no
   rational points over random values of an independent set, the pipeline now takes algebraic points numerically.
   The choice of minor is all they affect: the cut and its dimension check stay exact. The block then passes.
3. **Singular timeouts** (same quiver, block 100,010,001, a dimension-8 component). The pipeline now checks
   dim std(J) ≤ r before decomposing (`v2_55_pipeline.py`, same results and faster on small blocks). The driver
   treats a timeout as "try the case split", and the Singular timeout can be raised (SING_TMO).

The cost is dominated by the largest blocks: the zero block of 1⇛2⇇3 + 1→3 alone took 45 minutes, with 23
components up to dimension 8.

## 2. Duality (proved; `analysis/v2_66_duality.py`, `duality.txt`)

Let r: kQ → kQ^op reverse paths, r(a₁⋯a_k) = a_k⋯a₁ (with the arrows reversed). It is an anti-isomorphism of
algebras. Deleting a letter commutes with reversal up to swapping the tensor factors:

  Δ_op ∘ r = τ ∘ (r⊗r) ∘ Δ.

**Lemma D.** P ↦ rPr⁻¹ is a bijection Folds(kQ) → Folds(kQ^op) over every commutative ring. It conjugates the census
group and the class weightings, so local-torus classes correspond.

*Proof.*
- Idempotency is conjugation-invariant.
- The image and kernel of rPr⁻¹ are r(im P) and r(ker P), and an anti-isomorphism maps subalgebras to subalgebras.
- For the Δ-identity: Δ_op rPr⁻¹ = τ(r⊗r)ΔPr⁻¹ = τ(r⊗r)(P⊗P)Δr⁻¹. Since τ commutes with X⊗X for every X, this equals
  (rPr⁻¹⊗rPr⁻¹)τ(r⊗r)Δr⁻¹ = (rPr⁻¹⊗rPr⁻¹)Δ_op.
- Class substitutions and weightings act letter by letter, so they commute with r. ∎

**Checks.**
- The script checks both identities entry by entry over ℤ for pair_after → pair_before, pair_srcout → pair_tgtin,
  and the self-dual K₂, K₃, 1→2⇉3→4 and the double fork.
- It maps the q = 3 representatives of pair_after onto folds of pair_before, with 0 failures.
- It explains three equalities seen in the census, which had no explanation before:
  - pair_after and pair_before: 71 classes each;
  - pair_srcout and pair_tgtin: 101 each;
  - their vertex blocks, which agree after relabelling (checked).

So Q₀ = 1⇉2⇇3 + 1→3 proves 3ᶜ for 1⇇2⇉3 + 3→1 (diverging classes with joined ends) as well, and Q₁ and Q₂ for their
opposites.

## 3. Where the vertex blocks come from

### 3.1 The collapse lemma (proved; `analysis/v2_66_retraction.py`, `retraction.txt`)

Let Q^s be Q with every parallel class replaced by one arrow. Define:
- π: kQ → kQ^s, sending every arrow of a class to the class's arrow;
- ι: kQ^s → kQ, sending it back to a chosen arrow of the class.

Both are algebra maps, since they substitute letters. Both commute with Δ, because deleting one letter commutes with
substituting letters: Δ^s π = (π⊗π)Δ and Δι = (ι⊗ι)Δ^s. And πι = id.

**Lemma R.** For every fold P of kQ^s, ιPπ is a fold of kQ, over every commutative ring. P ↦ ιPπ is injective, and
it keeps the vertex block.

*Proof.*
- **Idempotent:** (ιPπ)² = ιP(πι)Pπ = ιPπ.
- **Image:** ι(im P) is a subalgebra, because ι is an algebra map.
- **Kernel:** ι is injective, so ker ιPπ = π⁻¹(ker P), the preimage of a subalgebra under an algebra map.
- **Δ-identity:** Δ ιPπ = (ι⊗ι)Δ^s Pπ = (ι⊗ι)(P⊗P)Δ^s π = (ι⊗ι)(P⊗P)(π⊗π)Δ = (ιPπ⊗ιPπ)Δ.
- **Injective:** π(ιPπ)ι = P.
- **Vertex block:** π and ι fix the vertex idempotents, and the vertex coordinates of ιPπ(e_i) = ι(P(e_i)) are those
  of P(e_i). ∎

**Checks.** The script checks the five hypotheses entry by entry over ℤ on 15 pairs (Q, Q^s), including the
obstructed 1⇉2⇇3. On all 15 it also lifts every fold of Q^s at q = 3 to a fold of Q: 0 failures, vertex blocks kept.

**Corollary.** For every quiver Q, VB(Q^s) ⊆ VB(Q), where VB(Q) is the set of vertex blocks of folds of Q.

### 3.2 Conjecture S

**Conjecture S.** If Q is unobstructed (in the sense of 3ᶜ), then VB(Q) = VB(Q^s), and every vertex block of Q^s is
the block of exactly one fold of Q^s. So the vertex-block strata of Folds(Q) are indexed by Folds(Q^s), and the
stratum of P^s contains the fold ιP^sπ.

**Evidence** (`multiplicity/`, `census_*_3.txt`, `*.blockcensus.txt`):

| Q | Q^s | |VB(Q)| = |VB(Q^s)| (char 0, as sets) | folds of Q^s at q = 3 | vdim of the vertex ideal (Q, Q^s) |
|---|---|---|---|---|
| K₂ (and K₃, K₄) | A₂ = 1→2 | 6 | 6 | 8, 6 |
| 1⇉2→3 (and 1⇛2→3, 1→2⇉3) | A₃ | 23 | 23 | 27, 23 |
| 1⇉2, 1→3 (and 1⇛2, 1→3) | 1→2, 1→3 | 27 | 27 | 30, 27 |
| 1⇉2, 3→2 | 1→2←3 | 27 | 27 | 30, 27 |
| 1⇉2⇇3 + 1→3 (and 1⇛2⇇3 + 1→3) | 1→2←3, 1→3 | 21 | 21 | 29, 21 |
| 1⇉2→3→4, 1→2⇉3→4 | A₄ | 100 | 100 | 116 / 106, 100 |
| 1⇉2, 2→3, 2→4 (and m = 3) | its simplification | 148 | 148 | 159, 148 |
| 1⇉2, 1→3, 4→2 (and m = 3) | its simplification | 123 | 123 | 127, 123 |
| Q₁ = Q₀ + 2→4 | its simplification | 90 | 90 | 105, 90 |
| Q₂ = Q₀ + 3→4 | its simplification | 99 | 99 | 119, 99 |
| 1⇉2→3→4→5 | A₅ | 563 | 563 | 627, — |
| **1⇉2⇇3 (obstructed)** | 1→2←3 | **1-dimensional** vs 27 | — | — |

Three things in this table:
- **Q^s is radical and Q is not.** Q^s's vertex ideal has vdim equal to its number of points, while Q's carries
  extra non-reduced structure from the class differences. That structure does not depend on m ≥ 2 (prediction 87).
- **Folds of Q^s are their vertex blocks.** At q = 3, Folds(𝔽₃Q^s) has exactly |VB| elements in every case
  computed (and the same counts at q = 5 for A₂, A₃, A₄ and 1→2←3 + 1→3). This is consistent with v2.16: every fold
  of a quiver without parallel arrows is canonical. The path quivers give 6, 23, 100, 563 for A₂ … A₅.
- **The obstructed case is where S fails.** Theorem 3 (v2.27) says the converging family of 1⇉2⇇3 has the vertex
  block as its invariant. In S's terms, the unjoined converging classes create vertex blocks that no fold of Q^s
  has. Joining the ends (Q₀) removes exactly those.

**Why it matters.**
- **For the procedure.** Step 2, the vertex elimination over ℤ, is the step that needed direct work at p = 2. Under
  S it becomes "VB(Q) ⊆ VB(Q^s)". Lemma R already gives ⊇.
- **For a structural induction (Program 3).** The strata of Folds(Q) get a canonical base point ιP^sπ. The question
  inside each stratum becomes a relative one over that point: the twist variety of v2.32, with an explicit centre.

### 3.3 What S reduces to (Lemma R′, proved; `analysis/v2_66_kerpi.py`, `kerpi.txt`)

**Lemma R′.** Let P be a fold of kQ with P(ker π) ⊆ ker π. Then P^s = πPι is a fold of kQ^s, with the vertex block of P.
Here ker π is the ideal generated by the class differences a − a′.

*Proof.*
- **A key identity.** For every x, x − ιπx ∈ ker π, so πP(x − ιπx) = 0. Hence πPιπ = πP.
- **Idempotent:** (πPι)² = πP(ιπ)Pι = πPPι = πPι.
- **Image:** π(im P), the image of a subalgebra under an algebra map.
- **Kernel:** it is ι⁻¹(P⁻¹(ker π)). Write x = (x − Px) + Px to see that P⁻¹(ker π) = ker P + ker π. That is a
  subalgebra plus an ideal, hence a subalgebra, and so is its preimage under ι.
- **Δ-identity:** Δ^s πPι = (π⊗π)(P⊗P)(ι⊗ι)Δ^s. Since (π⊗π)(P⊗P) = (πPι⊗πPι)(π⊗π) by the key identity, and
  πι = id, this is (P^s⊗P^s)Δ^s. ∎

Together with Lemma R: **VB(Q) = VB(Q^s) if and only if every vertex block of Q is the block of a fold that preserves
ker π.**
- Class substitutions with column sums 1 satisfy πσ = π, and the weightings scale whole classes, so both preserve
  ker π. So the property is a property of census classes.
- Not every fold has it, so S is an existence statement:

| quiver | folds or classes at q | preserve ker π | blocks | blocks with a ker π-preserving fold |
|---|---|---|---|---|
| K₂ | 20 classes, q = 3 and 5 | 14 | 6 | 6 |
| 1⇉2→3 | 219 folds, q = 3 | 99 | 23 | 23 |
| 1⇉2, 1→3 | 105 G-orbits, q = 3 | 71 | 27 | 27 |
| Q₀ | 407 / 513 G-orbits, q = 3 / 5 | 184 / 198 | 21 | 21 |
| **1⇉2⇇3** (obstructed) | 6376 / 59972 folds, q = 3 / 5 | 2128 / 12612 | 31 / 39 = 4q + 19 | **27** = VB(1→2←3) |

The obstructed row splits exactly: 27 blocks come from 1→2←3, and the other 4(q − 2) blocks carry only folds that
move ker π. The count 4(q − 2) is Theorem 3's converging family (449 + 4(q − 2) orbits, v2.27). So Theorem 3's
modulus is precisely the failure of ker π-stability at the vertex level.

**What a proof of S needs.**
- **The goal.** In each vertex-block stratum of an unobstructed Q, find a fold preserving ker π.
- **A natural route is degeneration.** The substitution σ_t: a_i ↦ a₁ + t(a_i − a₁) lies in the census group for
  t ≠ 0, and σ₀ = ιπ. It keeps the vertex block. The limit of σ_t P σ_t⁻¹ as t → 0 exists exactly when P is filtered
  by the powers of ker π, and it is then a ker π-preserving fold with the same block.
- **So S follows from:** every stratum of an unobstructed quiver contains a fold filtered by the powers of ker π.
- **Where the hypothesis must enter.** For 1⇉2⇇3 this fails on Theorem 3's blocks, so a proof has to use the joined
  ends there, as v2.52 did.

## 4. Next

1. **Prove S, or the inclusion VB(Q) ⊆ VB(Q^s), for unobstructed Q.** Start with one class (Conjecture 3ᵖ's
   setting), where the vertex ideal's extra structure is visible in K₂ (vdim 8 vs 6).
2. **Finish 1⇉2→3→4→5** (563 blocks; running: no failure in the first 424 blocks).
3. **Program 3**, as carried from v2.65: finiteness in Lemma B, the mixed kind, the λ = 1 bound. Conjecture S gives
   the leaf step a base point in every stratum.

## Files

- **analysis/:**
  - v2_66_prove.py (end-to-end driver, resumable with RESUME=1);
  - v2_66_vertex_badprimes.py, v2_66_table.py;
  - v2_66_duality.py, v2_66_retraction.py, v2_66_multiplicity.py, v2_66_kerpi.py.
- **Changed:**
  - v2_53_block_census.py (one-prime input);
  - v2_55_pipeline.py: numeric points to choose minors when no rational point is found, dim std(J) checked before
    minAssGTZ, SING_TMO;
  - pcensus.py: multiplicity variants and simplifications (K4, A2–A5, K3_*, S_*), and symmetry_plan for quivers
    without parallel arrows (no group generators).
- **logs/v2_66/:**
  - PREDICTIONS.md;
  - one folder per quiver: summary.json, integral.txt, vertex_integral.txt, vertex_char0.txt, vertex_char2.txt
    where needed, vertex_badprimes.json, blocks.txt, pipeline.txt, tally.txt, goodprimes.txt, bound.txt,
    certs.jsonl.xz;
  - *.blockcensus.txt, census_*_3.txt, reps/;
  - duality.txt, retraction.txt, kerpi.txt, multiplicity/.
