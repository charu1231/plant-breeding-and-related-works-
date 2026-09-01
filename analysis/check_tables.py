import csv

csv_rows = {}
with open('output/02_sim_summary.csv') as f:
    for r in csv.DictReader(f):
        csv_rows[(float(r['rG']), float(r['h2']))] = [
            round(float(r['taskA_naive_mean']), 3),
            round(float(r['taskA_mt_mean']), 3),
            round(float(r['taskA_rn_mean']), 3),
            round(float(r['taskB_single_mean']), 3),
            round(float(r['taskB_mt_mean']), 3),
        ]

md = open('manuscript/manuscript.md').read()

def table_after(header):
    seg = md.split(header, 1)[1]
    out = {}
    for ln in seg.splitlines():
        s = ln.strip()
        if not s.startswith('|'):
            if out:
                break
            continue
        cells = [c.strip() for c in s.strip('|').split('|')]
        if not cells or cells[0].startswith('-') or cells[0] == 'rG':
            continue
        try:
            key = (float(cells[0]), float(cells[1]))
            out[key] = [float(c.replace('**', '')) for c in cells[2:]]
        except (ValueError, IndexError):
            break
    return out

A = table_after('| rG | h\u00b2 | Naive mean | MT-GBLUP | RN-GBLUP |')
B = table_after('| rG | h\u00b2 | Single-env GBLUP | MT-GBLUP |')

errs = []
for k, v in csv_rows.items():
    if A.get(k) != v[:3]:
        errs.append((k, 'A', A.get(k), v[:3]))
    if B.get(k) != v[3:]:
        errs.append((k, 'B', B.get(k), v[3:]))

print('cells checked:', len(csv_rows) * 5, '| mismatches:', len(errs))
for e in errs:
    print('  ', e)
print('RESULT:', 'PASS - manuscript tables exactly match CSV' if not errs else 'FAIL')
