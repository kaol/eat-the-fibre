import sys, time
exec(open("q123s.py").read().split("\nt0 = time.time(); done")[0])
def fibre_prof(PK):
    P = np.zeros((d, d), dtype=np.int64); PKfull = Kmat @ PK % p
    for i in range(3): P[:, i] = PKfull[:, i]
    Pdelta = PKfull[:, 3]; La, Lb = lift(P, 3), lift(P, 5)
    if La is None or Lb is None: return None
    def set_a2(Pb): Pb[:, :, 4] = (Pb[:, :, 3] - Pdelta) % p
    c = [0, 0, 0, 0]; n = 0
    A = extend(P, 3, La, 4, set_a2); c[0] = len(A)
    for Pa in A:
        B = extend(Pa, 5, Lb, 5); c[1] += len(B)
        for Pab in B:
            L6, L7 = lift(Pab, 6), lift(Pab, 7)
            if L6 is None or L7 is None: continue
            C = extend(Pab, 6, L6, 6); c[2] += len(C)
            for P6 in C: n += len(extend(P6, 7, L7, 7))
    return c, n
big = []
for PK in Kfolds:
    t = time.time(); r = fibre_prof(PK); dt = time.time() - t
    if r and dt > 0.5 * (p / 3) ** 4: big.append((dt, r, PK.tolist()))
for b in sorted(big, key=lambda x: -x[0])[:6]: print(f"t={b[0]:.1f}s survivors per stage {b[1][0]} fibre {b[1][1]}  PK={b[2]}")
