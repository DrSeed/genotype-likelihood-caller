import os, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.makedirs("figures", exist_ok=True); os.makedirs("results", exist_ok=True)
rng = np.random.default_rng(0)
ERR = 0.01
def gl(alt, dp):
    # binomial genotype likelihoods for 0/0, 0/1, 1/1 given alt reads out of dp
    from math import comb
    ps = {"0/0": ERR, "0/1": 0.5, "1/1": 1 - ERR}
    return {g: comb(dp, alt) * p**alt * (1-p)**(dp-alt) for g, p in ps.items()}
def call(alt, dp):
    L = gl(alt, dp); return max(L, key=L.get)
# accuracy vs depth
depths = [2, 4, 6, 8, 10, 15, 20, 30, 50]
truth = ["0/0", "0/1", "1/1"]; pt = {"0/0": ERR, "0/1": 0.5, "1/1": 1-ERR}
acc = []
for dp in depths:
    ok = 0; N = 3000
    for _ in range(N):
        g = truth[rng.integers(3)]
        alt = rng.binomial(dp, pt[g]); ok += (call(alt, dp) == g)
    acc.append(ok / N)
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# example het site likelihoods across possible alt counts at depth 20
dp = 20; alts = np.arange(dp+1)
for g, c in [("0/0","#C44E52"),("0/1","#4C72B0"),("1/1","#55A868")]:
    ax[0].plot(alts, [gl(a, dp)[g] for a in alts], "o-", ms=3, color=c, label=g)
ax[0].set_xlabel("alt reads (of 20)"); ax[0].set_ylabel("likelihood"); ax[0].set_title("Genotype likelihoods at depth 20"); ax[0].legend(fontsize=8)
ax[1].plot(depths, acc, "o-", color="#4C72B0"); ax[1].set_xlabel("read depth"); ax[1].set_ylabel("calling accuracy"); ax[1].set_ylim(0.5,1.02); ax[1].set_title("Accuracy rises with depth")
fig.suptitle("Genotype-likelihood calling (demo data)"); fig.tight_layout(); fig.savefig("figures/demo.png", dpi=140)
open("results/summary.csv","w").write("depth,accuracy\n"+"\n".join(f"{d},{a:.3f}" for d,a in zip(depths,acc))+"\n")
print(f"acc @10x={acc[4]:.3f} @30x={acc[7]:.3f}"); print("ok")
