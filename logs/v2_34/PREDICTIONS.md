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

## Addendum (recorded after the first runs, before the joined-ends test)

Outcome so far: prediction 1 fails at first order (1⇉2, 1→3 has points with e > 0). Refined, post hoc, criterion: for
each behaviour fibre, the *second-order* excess, max over the fibre of (dim S + dimension of the excess directions
that extend to second order along the orbit) − max over the fibre of dim S. On 1⇉2, 1→3 (𝔽₁₁) it is 0 on all 91 fibres;
on 1⇉2⇇3 (𝔽₃) it is 1 on exactly the four Theorem 3 fibres and 0 on the other 315.

Prediction 4 (the joined-ends case): on 1⇉2⇇3 + 1→3 over 𝔽₃, using the classes of Inn ⋊ Sub (finer than behaviours,
since the substitutions a_i ↦ a_i + λ c b_j are left out, which can only make the test stricter), the second-order excess
is 0 on every class.

## Outcomes

1. Refuted: e > 0 at representatives of unobstructed quivers (1⇉2, 1→3 / 𝔽₁₁: 15 of 137).
2. Holds: e > 0 somewhere on each obstructed quiver tested (but it does not discriminate).
3. Refuted: b > 0 widely, on both sides of the obstruction.
4. Holds at the behaviour level (all 18 behaviours containing the 20 apparent excess classes have second-order excess 0);
   fails at the finer Inn ⋊ Sub class level.
