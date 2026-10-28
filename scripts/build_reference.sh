#!/usr/bin/env bash
#
# Constroi a referencia compartilhada UMA VEZ, sob um lock exclusivo.
#
#   bash scripts/build_reference.sh [args extras para o snakemake]
#
# Por que o flock fica AQUI e nao dentro de uma rule:
# o Snakemake apaga os outputs existentes de um job imediatamente antes de
# executa-lo. Dois runners que planejaram enquanto a referencia ainda nao
# existia fariam o segundo apagar o indice que o primeiro acabou de concluir.
# Um lock dentro do shell: da rule nao impede isso, porque a remocao acontece
# antes do shell rodar. O lock precisa envolver o processo snakemake inteiro.
# O lock interno do Snakemake (.snakemake/locks) e por workdir e nao protege
# um diretorio compartilhado entre workdirs/maquinas diferentes.

set -euo pipefail

cd "$(dirname "$0")/.."

# Mesma logica de camadas que a analise usa: config.yaml e, se existir,
# config.local.yaml por cima (ou o arquivo apontado por SNAKEMAKE_CONFIG).
CONFIGS=("config.yaml")
if [[ -n "${SNAKEMAKE_CONFIG:-}" ]]; then
    CONFIGS+=("$SNAKEMAKE_CONFIG")
elif [[ -f config.local.yaml ]]; then
    CONFIGS+=("config.local.yaml")
fi

read -r REF_DIR GENOME_ID < <(python3 -I - "${CONFIGS[@]}" <<'PY'
import sys, yaml
ref = {}
for path in sys.argv[1:]:
    with open(path) as fh:
        ref.update((yaml.safe_load(fh) or {}).get("reference", {}))
print(ref["dir"], ref["genome_id"])
PY
)

TARGET="${REF_DIR}/${GENOME_ID}"

if [[ ! -d "$REF_DIR" ]]; then
    echo "ERRO: $REF_DIR nao existe." >&2
    echo "Crie com:  sudo mkdir -p $REF_DIR && sudo chown \"\$USER\":\"\$USER\" $REF_DIR" >&2
    exit 1
fi

mkdir -p "$TARGET"

exec 9>"${TARGET}/.build.lock"
if ! flock -n 9; then
    echo "Outro build ja esta em andamento em ${TARGET}. Abortando." >&2
    exit 1
fi

echo "Lock obtido em ${TARGET}/.build.lock"
echo "Construindo referencia em ${TARGET}"

snakemake -s reference.smk \
    --configfile "${CONFIGS[@]}" \
    --cores "${SNAKEMAKE_CORES:-20}" \
    --resources mem_mb="${SNAKEMAKE_MEM_MB:-64000}" \
    reference_all "$@"

echo "Pronto. Sentinela e MANIFEST:"
ls -la "$TARGET" | grep -E "\.done_|MANIFEST" || true
