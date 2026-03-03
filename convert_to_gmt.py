#/usr/bin/env python

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: python convert_to_gmt.py input.tsv output.gmt")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

go_terms = defaultdict(lambda: {"desc": "", "genes": set()})

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            continue

        #gene = parts[0]
        gene = parts[0].replace("|", "_").replace(".", "_") # replace '|' and '.' in bv-brc fig id with '_'.
        go_entries = parts[1].split(';')
        for entry in go_entries:
            if '|' not in entry:
                continue
            go_id, desc = entry.split('|', 1)
            go_terms[go_id]["desc"] = desc
            go_terms[go_id]["genes"].add(gene)

with open(output_file, 'w') as out:
    for go_id, data in sorted(go_terms.items()):
        line = [go_id, data["desc"]] + sorted(data["genes"])
        out.write('\t'.join(line) + '\n')

print(f"✅ Successfully wrote GMT to: {output_file}")

