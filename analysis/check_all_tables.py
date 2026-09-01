"""Verify that every result table/number in the manuscript matches the raw
CSV/JSON outputs. Uses standard (round-half-up) rounding."""
import csv
import json
from decimal import Decimal, ROUND_HALF_UP

OK = True


def r3(x):
    return float(Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def report(name, flag):
    global OK
    OK = OK and flag
    print(("PASS: " if flag else "FAIL: ") + name)


md = open("manuscript/manuscript.md").read()

# ---------------- main simulation tables (02) ----------------
rows = {}
with open("output/02_sim_summary.csv") as f:
    for r in csv.DictReader(f):
        rows[(float(r["rG"]), float(r["h2"]))] = [
            r3(r["taskA_naive_mean"]), r3(r["taskA_mt_mean"]), r3(r["taskA_rn_mean"]),
            r3(r["taskB_single_mean"]), r3(r["taskB_mt_mean"])]


def table_after(header, keycols, valcols):
    seg = md.split(header, 1)[1]
    out = {}
    for ln in seg.splitlines():
        s = ln.strip()
        if not s.startswith("|"):
            if out:
                break
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if not cells or cells[0].startswith("-"):
            continue
        try:
            key = tuple(float(cells[i]) for i in range(keycols))
            out[key] = [float(c.replace("**", "")) for c in cells[keycols:keycols + valcols]]
        except (ValueError, IndexError):
            break
    return out


A = table_after("| rG | h\u00b2 | Naive mean | MT-GBLUP | RN-GBLUP |", 2, 3)
B = table_after("| rG | h\u00b2 | Single-env GBLUP | MT-GBLUP |", 2, 2)
bad = sum(1 for k, v in rows.items() if A.get(k) != v[:3] or B.get(k) != v[3:])
report(f"Main simulation tables (Task A + B) — {bad} mismatches", bad == 0)

# ---------------- structured (04) ----------------
s4 = {}
with open("output/04_structured_summary.csv") as f:
    for r in csv.DictReader(f):
        s4[(float(r["kappa"]), float(r["h2"]))] = [
            r3(r["taskA_naive_mean"]), r3(r["taskA_mt_mean"]), r3(r["taskA_rn_mean"])]
S4 = table_after("| \u03ba | h\u00b2 | Naive | MT-GBLUP | RN-GBLUP | MT (oracle) | RN (oracle) |", 2, 3)
bad = sum(1 for k, v in s4.items() if S4.get(k) != v)
report(f"Structured G\u00d7E table — {bad} mismatches", bad == 0)

# ---------------- unbalanced (05) ----------------
s5 = {}
with open("output/05_unbalanced_summary.csv") as f:
    for r in csv.DictReader(f):
        s5[(float(r["rg"]),)] = [
            r3(r["taskA_naive_mean"]), r3(r["taskA_mt_mean"]),
            r3(r["taskB_single_mean"]), r3(r["taskB_mt_mean"])]
S5 = table_after("| rG | Task A: naive | Task A: MT | Task B: single | Task B: MT |", 1, 4)
bad = sum(1 for k, v in s5.items() if S5.get(k) != v)
report(f"Unbalanced-data table — {bad} mismatches", bad == 0)

# ---------------- factor-analytic (07) ----------------
s7 = {}
with open("output/07_fa_summary.csv") as f:
    for r in csv.DictReader(f):
        s7[r["scenario"]] = [r3(r[f"{m}_mean"]) for m in
                             ["naive", "fa1", "fa2", "fa3", "rn", "mt"]]
# manuscript FA table has columns: naive FA1 FA2 FA3 RN MT (7 value cols after scenario)
S7 = table_after("| Scenario | varExp(1/2/3) | naive | FA1 | FA2 | FA3 | RN | MT |", 1, 6)
# note: the FA table's first column is a string, not float; parse by string key
S7s = {}
for ln in md.split("| Scenario | varExp(1/2/3) | naive | FA1 | FA2 | FA3 | RN | MT |", 1)[1].splitlines():
    s = ln.strip()
    if not s.startswith("|"):
        if S7s:
            break
        continue
    cells = [c.strip() for c in s.strip("|").split("|")]
    if not cells or cells[0].startswith("-") or cells[0] == "Scenario":
        continue
    try:
        S7s[cells[0]] = [float(c.replace("**", "")) for c in cells[2:8]]
    except (ValueError, IndexError):
        break
# map scenario keys (manuscript uses "CS, rG = 0.3" etc.) to CSV keys
keymap = {"CS, rG = 0.3": "rg0.3_h2_0.4", "CS, rG = 0.6": "rg0.6_h2_0.4",
          "RN, \u03ba = 0.3": "kappa0.3_h2_0.6", "RN, \u03ba = 0.8": "kappa0.8_h2_0.6"}
bad = sum(1 for k, csvk in keymap.items() if S7s.get(k) != s7.get(csvk))
report(f"Factor-analytic rank table — {bad} mismatches", bad == 0)

# ---------------- real data (03) ----------------
real = json.load(open("output/03_real_validation.json"))
checks = {
    "MT 0.379": abs(real["taskA"]["mt"] - 0.3792) < 0.01,
    "naive 0.313": abs(real["taskA"]["naive"] - 0.3132) < 0.01,
    "RN 0.262": abs(real["taskA"]["rn"] - 0.2623) < 0.01,
    "TaskB MT 0.443": abs(real["taskB"]["mt"] - 0.4429) < 0.01,
    "TaskB single 0.440": abs(real["taskB"]["single"] - 0.4405) < 0.01,
    "rG(E2,E4)=0.97": abs(real["typeB_correlation"][1][2] - 0.97) < 0.02,
    "rG(E4,E5)=0.94": abs(real["typeB_correlation"][2][3] - 0.94) < 0.02,
}
for name, flag in checks.items():
    report(f"Real-data number {name}", flag)

print("\nOVERALL:", "ALL CONSISTENT" if OK else "SOME MISMATCHES")
