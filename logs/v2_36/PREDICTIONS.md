# v2.36 pre-registration (committed before the runs on 1⇉2⇇3 and 1⇉2⇇3 + 1→3)

Inner-curve search (analysis/family_curves.py): for every liftable local-orbit representative P₀ and every direction
x = α δ_a + β δ_b with (α, β) among small integer pairs up to sign and scale, test whether
P_t = Inn(1 + t·x) P₀ Inn(1 + t·x)⁻¹ is a fold for all t (over ℚ[t]); if so, count the local-torus orbits its
𝔽_q-points meet. A curve meeting ≥ 3 local orbits is a family.

1. Control, 1⇉2⇇3 over 𝔽₅: at least two families are found (F2 and F4 are pure inner curves along δ_b − δ_a, v2.33),
   every family found lies in one of the two growing signature groups of v2.35, and none elsewhere.
2. Joined ends, 1⇉2⇇3 + 1→3 over 𝔽₅: no curve meets ≥ 3 local orbits.
