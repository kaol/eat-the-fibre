# v2.72 predictions (committed before the runs)

p = 2 needs a formulation first: vertex blocks have entries in {0, 1, −1}, and 1 ≡ −1 mod 2, so char 0 blocks
merge. Step one is the census itself at q = 2 (pcensus.py, local-torus classes), never run on these quivers, compared
with the char 0 geometric orbit counts (v2.67) and the q = 3 census.

109. On the four quivers without −1 entries in any vertex block (K₂, K₃, Q₀, 1⇛2⇇3 + 1→3), the census at q = 2
     equals the char 0 count (20, 34, 366, 752).
110. On the other 13 quivers the census at q = 2 differs from the char 0 count, and is smaller: folds from char 0
     blocks that agree mod 2 lie in one 𝔽_2-stratum.
111. The q = 2 count equals the number of geometric orbits in characteristic 2, computed by Lemma L at p = 2 on the
     merged blocks (char-2 trees, Lemma N), wherever that computation finishes.
