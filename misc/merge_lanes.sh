#!/usr/bin/env bash
# Funde lanes de amostras com mais de um run, para o padrao samples/{amostra}_{1,2}.fastq.gz
# que o Snakefile espera. Amostras de run unico nao sao tocadas.
#
#   bash misc/merge_lanes.sh metadata/amostras_tierA_master.csv samples
#
# O nome final da amostra e o PRIMEIRO run accession em ordem alfabetica (coluna sample_fastq),
# exatamente como consta nos metadata/exp_design_*.csv.
set -euo pipefail

master="${1:?uso: bash misc/merge_lanes.sh <amostras_master.csv> [dir_samples]}"
dir="${2:-samples}"

python3 - "$master" "$dir" <<'PY'
import csv, os, subprocess, sys
master, d = sys.argv[1], sys.argv[2]
with open(master, newline='') as fh:
    for row in csv.DictReader(fh):
        runs = row["runs"].split(";")
        if len(runs) < 2:
            continue
        alvo = row["sample_fastq"]
        for mate in ("1", "2"):
            partes = [os.path.join(d, f"{r}_{mate}.fastq.gz") for r in runs]
            faltando = [p for p in partes if not os.path.exists(p)]
            if faltando:
                print(f"PULANDO {alvo} mate {mate}: ausentes {faltando}", file=sys.stderr)
                continue
            saida = os.path.join(d, f"{alvo}_{mate}.fastq.gz")
            tmp = saida + ".tmp"
            print(f"fundindo {len(partes)} lanes -> {saida}")
            with open(tmp, "wb") as out:
                subprocess.run(["cat"] + partes, stdout=out, check=True)
            os.replace(tmp, saida)
            # remove apenas as partes que NAO sao o alvo
            for p in partes:
                if os.path.abspath(p) != os.path.abspath(saida):
                    os.remove(p)
PY

echo "Pronto. Confira a contagem:"
ls "$dir"/*_1.fastq.gz | wc -l
