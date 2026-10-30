#!/usr/bin/env bash
#
# Baixa corridas paired-end do SRA listadas na primeira coluna de um CSV.
#
#   bash misc/download_sra.sh <lista.csv> [dir_saida] [dir_temp]
#
# O dir_temp precisa ficar em disco: /tmp neste cluster e tmpfs em RAM
# (ver plano de rede, secao 13.7) e o fasterq-dump estoura a memoria la.
#
# Variaveis de ambiente opcionais:
#   RETRIES=3        tentativas por run antes de desistir
#   RETRY_WAIT=60    segundos entre tentativas
#   FASTERQ_MEM=32G  teto de RAM do fasterq-dump (antes: a RAM TOTAL do no)
#   VERIFY=0|1|2     0 = so checa par nao-vazio (padrao)
#                    1 = + integridade gzip (pigz -t)
#                    2 = + confere a coluna read_count do CSV, se existir
#
# Um run que falha NAO aborta a fila: e tentado RETRIES vezes, registrado em
# <dir_saida>/../download_falhas.txt e o script segue. O exit code e 1 se
# algum run ficou para tras, 0 se todos vieram.

set -uo pipefail          # sem -e: o laco precisa sobreviver a um run ruim

csv="${1:?uso: bash misc/download_sra.sh <lista.csv> [dir_saida] [dir_temp]}"
outdir="${2:-samples}"
tmpdir="${3:-/mnt/scratch-shared/$USER/sra_tmp}"

retries="${RETRIES:-3}"
retry_wait="${RETRY_WAIT:-60}"
verify="${VERIFY:-0}"

# acima de ~16 threads o fasterq-dump nao ganha mais nada
ncores=$(nproc)
threads=$(( ncores > 16 ? 16 : ncores ))

# Teto explicito. Passar a RAM total (free -m) e perigoso num no compartilhado:
# o fasterq-dump tenta usa-la e concorre com os jobs do Snakemake.
memory="${FASTERQ_MEM:-32G}"

# --progress so faz sentido em terminal; redirecionado, enche o log de \b
progress_flag=""
[ -t 1 ] && progress_flag="--progress"

mkdir -p "$outdir" "$tmpdir"
falhas="$(dirname "$outdir")/download_falhas.txt"
: > "$falhas"

# read_count esperado por run, se o CSV trouxer essa coluna
declare -A expected_reads=()
if [ "$verify" -ge 2 ]; then
    col=$(awk -F',' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="read_count") {print i; exit}}' "$csv")
    if [ -n "${col:-}" ]; then
        while IFS=, read -r -a f; do
            [ -n "${f[0]:-}" ] && expected_reads["${f[0]}"]="${f[$((col-1))]}"
        done < <(awk -F',' 'NR>1' "$csv")
    else
        echo "AVISO: VERIFY=2 pedido mas o CSV nao tem coluna read_count" >&2
    fi
fi

baixar_um() {
    local id="$1"
    prefetch --output-directory "$outdir" "$id" || return 1
    # shellcheck disable=SC2086
    fasterq-dump $progress_flag \
        --threads "$threads" \
        --mem "$memory" \
        --temp "$tmpdir" \
        --outdir "$outdir" \
        "$outdir/$id" || return 1

    # single-end entra aqui como um unico arquivo; e erro, nao silencio
    if [ ! -s "$outdir/${id}_1.fastq" ] || [ ! -s "$outdir/${id}_2.fastq" ]; then
        echo "   $id: fasterq-dump nao produziu o par _1/_2 (single-end?)" >&2
        return 1
    fi
    pigz -f -p "$threads" "$outdir/${id}_1.fastq" "$outdir/${id}_2.fastq" || return 1
    return 0
}

verificar() {
    local id="$1"
    [ -s "$outdir/${id}_1.fastq.gz" ] && [ -s "$outdir/${id}_2.fastq.gz" ] || return 1
    if [ "$verify" -ge 1 ]; then
        pigz -t "$outdir/${id}_1.fastq.gz" && pigz -t "$outdir/${id}_2.fastq.gz" || {
            echo "   $id: gzip corrompido" >&2; return 1; }
    fi
    if [ "$verify" -ge 2 ] && [ -n "${expected_reads[$id]:-}" ]; then
        local obs
        obs=$(( $(pigz -dc -p "$threads" "$outdir/${id}_1.fastq.gz" | wc -l) / 4 ))
        if [ "$obs" != "${expected_reads[$id]}" ]; then
            echo "   $id: read_count divergente (esperado ${expected_reads[$id]}, obtido $obs)" >&2
            return 1
        fi
    fi
    return 0
}

samples=$(awk -F"," 'NR>1 {print $1}' "$csv")
total=$(printf '%s\n' ${samples} | grep -c . || true)
n=0; ok=0; pulados=0; ruins=0

for id in ${samples}; do
    n=$(( n + 1 ))
    if verificar "$id" 2>/dev/null; then
        echo "==> ${id} [${n}/${total}] ja presente, pulando"
        pulados=$(( pulados + 1 )); continue
    fi
    sucesso=0
    for (( try=1; try<=retries; try++ )); do
        echo "==> ${id} [${n}/${total}] tentativa ${try}/${retries}"
        if baixar_um "$id" && verificar "$id"; then sucesso=1; break; fi
        rm -rf "$outdir/$id" "$outdir/${id}_1.fastq" "$outdir/${id}_2.fastq"
        [ "$try" -lt "$retries" ] && { echo "    falhou; aguardando ${retry_wait}s"; sleep "$retry_wait"; }
    done
    rm -rf "$outdir/$id"
    if [ "$sucesso" -eq 1 ]; then
        ok=$(( ok + 1 ))
    else
        ruins=$(( ruins + 1 )); echo "$id" >> "$falhas"
        echo "    DESISTINDO de $id apos ${retries} tentativas"
    fi
done

echo "Pronto. FASTQ em ${outdir}/"
echo "Resumo: ${ok} baixados, ${pulados} ja presentes, ${ruins} falharam de ${total}."
if [ "$ruins" -gt 0 ]; then
    echo "Runs que falharam (relance com esta lista): $falhas" >&2
    exit 1
fi
