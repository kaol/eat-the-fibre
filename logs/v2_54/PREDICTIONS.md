# v2.54 predictions (committed before any computation)

Criterion (char 0): a connected algebraic group H acting on an irreducible variety Z has finitely many orbits if, for
every r < dim Z, the locus Z_r = {P ∈ Z : dim H·P ≤ r} = {rank [Lie H, P] ≤ r} has dimension ≤ r. Proof: on a
component Y of Z_r \ Z_{r−1} with dim Y = r, the orbits have dimension r, so each is open and dense in Y. Hence Y
holds a single orbit, and the rest of Y lies in Z_{r−1}. Induct.

An upper bound for dim Z_r suffices. N random integer combinations det(L·M(P)·R) of the (r+1)-minors of
M(P) = ([h_i, P])_i vanish on Z_r, so dim(Z ∩ V(combinations)) ≥ dim Z_r. Here L is (r+1)×m, R is d²×(r+1), and
N = dim Z + 2.

53. For every one of the 118 components Z of Folds(1⇉2⇇3 + 1→3) (v2.53) and every r < dim Z, the bound gives
    dim Z_r ≤ r for H_Z = ⟨Aut_Δ°, T_{W_Z}⟩. So every component has finitely many H_Z-orbits, and Folds(k̄Q₀) is a
    finite union of census-group orbits. Moderate confidence: on boundary loci the local torus can be larger than
    T_{W_Z}, and H_Z alone might not be enough.
54. The total number of H_Z-orbits is then bounded by Σ_Z Σ_r (number of components of Z_r of dimension r), and this
    bound is at least 366.
