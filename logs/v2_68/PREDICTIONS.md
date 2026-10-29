# v2.68 predictions (committed before the runs)

Goal: uniformity in q (every characteristic) for block 100,001,001 of 1⇛2⇇3 + 1→3, the one block whose char 0
certificate is an orbit tree (v2.67) rather than v2.54 cuts. Method (analysis/v2_68_uniform.py), for every node Y
of the tree, over ℤ:
- C1 a strong Gröbner basis of I_Y (integral primitive generators); its leading coefficients;
- C2 the closure J_Y of I_Y + (f̃) under the integral vector fields, computed over ℤ (so its reduction mod p lies in
  the char p closure and vanishes on the points of Y_p with orbit dimension < dim Y, in any characteristic);
- C3 K_Y = the intersection over ℤ of the children's ideals, and ρ^k ∈ J_Y over ℤ for its generators ρ;
- C4 the local-torus weightings of Y compatible on I_Y over ℤ;
- C5 at the root: the stratum ideal contains a power of every generator of the intersection of the 6 roots.
For p not dividing a leading coefficient of any C1 basis (and of the C2–C5 bases where an exact ℤ-membership fails),
every node fibre Y_p is H-stable (flat closure of an H-stable variety under a smooth group scheme), and the orbit
tree specialises: finitely many classes over 𝔽̄_p, bounded independently of p by the degrees.

98. Every claim C1–C5 holds over ℤ, and the candidate bad primes form a finite set (prediction: a subset of the
    primes 2, 3, 5 that appear as leading coefficients in the other 20 blocks' certificates, v2.66).
99. At each candidate prime the claims hold when checked directly in characteristic p (the same computations over 𝔽_p).
