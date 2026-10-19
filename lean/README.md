# Lean 4 formalization (v2.29)

Lean 4 (v4.35.0) with Mathlib (v4.35.0, pinned in `lake-manifest.json`). Mathlib is used from its prebuilt cache; `.lake/`
(packages, oleans, ~7 GB) is git-ignored and never committed.

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean && ./build.sh        # lake exe cache get, then each module in turn (the kernel checks need ~5 GB each)
python3 gen/gen.py           # regenerate EatTheFibre/Gen/*.lean from the census engine (pcensus.build)
```

## What is proved

| file | content |
|---|---|
| `Basic.lean` | Coordinate algebras (structure constants `c`, arrow deletion `D`), `IsFold` (idempotent, image and kernel closed under products, `Δ(Pv) = P Δ(v) Pᵀ`). **Basis criterion** `isFold_iff_basis`: the conditions hold for all vectors iff they hold on basis vectors. That is the method of every Python check. **`IsFold.conj`**: conjugating by a Δ-compatible automorphism gives a fold. **`vertex_block`**: an automorphism fixing the vertices and preserving the radical does not change the vertex-to-vertex entries |
| `Poly.lean` | Computable integer polynomials; `eval` is a ring map into any commutative ring; `isZero p → eval t p = 0` |
| `Reflect.lean` | `check c D P : Bool` evaluates the four basis conditions on polynomial matrices. **`check_sound`**: `check = true` gives `IsFold` over every commutative ring at every parameter value. The check is split into rows, one `decide +kernel` each |
| `Gen/*.lean` | Generated from `pcensus.build`, i.e. the census engine's own basis, products and Δ. Families proved to be folds (`*_isFold`, with readable images `*_a2`, `*_kills` …), plus controls proved to fail (`*_check_fails`) |
| `Invariants.lean` | Distinct parameters give distinct orbits: `thm3_rigid_GT`, `thm4_rigid_GT` |
| `Audit.lean` | `#print axioms`: only `propext`, `Classical.choice`, `Quot.sound` (no `sorry`, no `native_decide`) |

| result (docs) | Lean |
|---|---|
| Theorem 1 / Theorem 2 normal form on 1⇉2⇉3 (v2.26, v2.27) | `Q_two_pairs.thm2_isFold` |
| Theorem 2 fails with the arrow 1→3 | `Q_two_pairs_short.thm2_check_fails` (fails at Δ of the arrow 1→3) |
| Theorem 3, families F1–F4 on 1⇉2⇇3 (v2.27) | `Q_pairs_in.thm3F{1,2,3,4}_isFold`, rigidity `thm3_rigid_GT` |
| Theorem 3 fails with the arrow 1→3 | `Q_pairs_in_short.thm3F{1,2,4}_check_fails` |
| Theorem 4 on 1⇉2⇉3 and on 1⇉2⇉3 + 1→3 (v2.28) | `Q_two_pairs.thm4_isFold`, `Q_two_pairs_short.thm4_isFold`, rigidity `Q_two_pairs.thm4_rigid_GT` |

Not formalized: the reduction from an arbitrary acyclic quiver to the subalgebra `W` in the proofs of Theorems 3 and 4.
That needs general path algebras. Census counts are not formalized either.
