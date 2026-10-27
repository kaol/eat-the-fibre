# v2.60 predictions (committed before the q = 7 run and the fibre checks)

Program 3 (leaf extension), Q₀ → Q₁ = Q₀ + ℓ: 2→4. kQ₁ = kQ₀ ⊕ J, J = paths through 4 (an ideal).
Observed at q = 3 and 5 (analysis/v2_60_leaf_census.py): of Q₁'s 1220 census classes, 617 preserve both J and kQ₀
("split"), 307 only kQ₀, 280 only J, 16 neither; the compression π P ι is a fold of Q₀ exactly in the first three
kinds; the 897 J-preserving classes lie over all 366 Q₀ classes, with fibre sizes {1: 47, 2: 150, 3: 129, 4: 37, 5: 3}.

66. The same decomposition and the same fibre-size distribution at q = 7.
67. Split folds have the form P = P̄ ⊕ P_J with P_J(e₄) = λe₄ (λ ∈ {0, 1}) and P_J(pℓ) = φ(p)ℓ + μ(p)e₄, where
    φ(p) = λ·P̄(p) + Σ μ(suffix)·P̄(prefix) over the deletions of one arrow of p (Δp = Σ prefix ⊗ suffix). Derived by
    hand from ΔP = (P⊗P)Δ; checked on every split census representative at q = 3, 5, 7.
68. Every Q₀ fold P̄ has the trivial lift P̄ ⊕ 0 (λ = 0, μ = 0), a fold in every characteristic; so Q₀'s classes inject
    into Q₁'s split classes (366 of the 617), and the other 251 split classes have λ = 1 or μ ≠ 0.

## Outcome

- **66 held**: q = 7 gives 617 / 307 / 280 / 16 and the same fibre distribution.
- **67 held**: (*) on every split representative (692, 890, 1152 at q = 3, 5, 7); proved by hand (doc §2).
- **68 held**: 366 trivial lifts (all of Q₀'s classes), and the other 251 split classes all have λ = 1 (228 with μ = 0,
  exactly over the Q₀ classes with P̄(B) ⊆ B; 23 with μ ≠ 0). λ = 0 forces the trivial lift (proved).
