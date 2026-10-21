# v2.34 pre-registration (committed before the rigidity computation is run)

At a fold P (over 𝔽_p, all linear algebra), with Der = derivations of kQ, Der_Δ = Δ-compatible derivations, E_w the
class-weight derivations and W_P = {w : P compatible with Δ_w}:
- T      = T_P Folds (Zariski tangent: kernel of the linearised fold equations, as in analysis/tangent.py)
- F      = [Der, P] ∩ T          (tangent to the behaviour fibre Aut(kQ)·P ∩ Folds)
- S      = [Der_Δ, P] + [E_{W_P}, P]   (tangent to the census-symmetry orbit; S ⊆ F is a sanity check)
- b(P)   = dim T − dim F   (behaviour excess: fold directions not produced by any automorphism)
- e(P)   = dim F − dim S   (fibre excess: automorphism directions staying among folds but not census symmetries)

Predictions, on the G-orbit representatives of the 21 censuses in logs/v2_32/reps.tar.gz:
1. Conjecture 3ᶜ, infinitesimal form: on every quiver without the obstruction (no consecutive classes of size ≥ 2, no
   converging/diverging pair with unjoined outer ends), e(P) = 0 at every representative.
2. On every quiver with the obstruction, e(P) ≥ 1 at some representative (the Theorem 2–4 directions).
3. b(P) = 0 at every representative of every census (behaviours infinitesimally rigid).
Failures of 3 can be non-reducedness of the fold scheme rather than moduli; they will be reported, not explained away.
