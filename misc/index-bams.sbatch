#!/usr/bin/env bash
#
# Gera o .bai de BAMs ordenados por coordenada (um BAM por tarefa de um job array).
# Roda no SLURM, desanexado da sessao (CLAUDE.md da tese, secao 4.7).
#
# Este arquivo e AUTOCONTIDO de proposito: o sbatch copia o script para o no de
# execucao, mas o /home do admin e local a cada no (nao e NFS). Nao chame
# nenhum outro arquivo de /home daqui.
#
#   LIST=/mnt/scratch-shared/ronald/bai_piloto/lists/lote.txt   # um caminho de BAM por linha
#   sbatch --array=0-80%8 --export=ALL,LIST=$LIST,PHASE=final misc/index-bams.sbatch
#
# Variaveis (todas opcionais, exceto LIST):
#   LIST        arquivo com os BAMs; a tarefa N processa a linha N+1
#   PHASE       rotulo gravado no log de vazao (default: run)
#   LOGDIR      onde gravar o log por tarefa (default: /mnt/scratch-shared/ronald/bai_piloto)
#   BAI_THREADS threads extras do samtools index (default: 4)
#   EVICT       1 = descarta do page cache do no o BAM antes de ler (medicao de vazao
#               a frio). Default 0.
#   WAIT_FILE   se definido, espera este arquivo existir (ate 10 min) antes de comecar.
#               Usado para o monitor limpar o cache do servidor NFS antes da medicao.
#
# Propriedades:
#   - idempotente: um BAM cujo .bai ja existe, e mais novo que o BAM e VALIDA e pulado
#   - o temporario mantem a extensao final (.tmp.bai) e o formato vai explicito (-b):
#     o htslib escolhe formato pela extensao (CLAUDE.md 4.1.1)
#   - confere a ferramenta antes de interpretar qualquer codigo de saida (4.2.1)
#   - o resultado e VERIFICADO (idxstats + contagem por contig), nao so o exit code
#   - escreve somente <bam>.bai ao lado do BAM (e o log por tarefa em LOGDIR)

#SBATCH --job-name=bai
#SBATCH --partition=gpu,cpu,bigmem
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%A_%a.out

set -uo pipefail

# O locale dos nos NAO e uniforme: no r550-1 o awk imprimiu "1793576119,139" (virgula decimal)
# e quebrou o parse do log de vazao. LC_ALL=C fixa ponto decimal e ordenacao em qualquer no.
export LC_ALL=C
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH
LOGDIR="${LOGDIR:-/mnt/scratch-shared/ronald/bai_piloto}"
PHASE="${PHASE:-run}"
THREADS="${BAI_THREADS:-4}"
EVICT="${EVICT:-0}"
HOST="$(hostname -s)"

# --- 1. a ferramenta existe? (distinguir "nao rodou" de "falhou") ---------------------
for t in samtools python3 stat awk; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERRO: '$t' ausente no PATH em $HOST; nada foi verificado" >&2; exit 127; }
done
[ -n "${LIST:-}" ] && [ -s "$LIST" ] || { echo "ERRO: LIST vazio ou ausente: '${LIST:-}'" >&2; exit 2; }
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] || { echo "ERRO: rode como job array (--array)" >&2; exit 2; }

BAM="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")"
[ -n "$BAM" ] || { echo "ERRO: LIST nao tem linha $((SLURM_ARRAY_TASK_ID + 1))" >&2; exit 2; }
[ -f "$BAM" ] || { echo "ERRO: BAM nao existe: $BAM" >&2; exit 2; }
case "$BAM" in *.bam) ;; *) echo "ERRO: nao e .bam: $BAM" >&2; exit 2;; esac

BAI="${BAM}.bai"
TMP="${BAM%.bam}.tmp.bai"          # termina em .bai: o formato nao muda pela extensao
mkdir -p "$LOGDIR/tasks"
REC="$LOGDIR/tasks/${PHASE}_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID}.tsv"

cleanup() { rm -f "$TMP"; }
trap cleanup EXIT

# --- 2. validacao de um indice: nao basta existir --------------------------------------
# Confere: idxstats abre, tem um registro por @SQ (+ a linha '*'), e a contagem por
# contig bate com 'samtools view -c' em tres contigs (inicio, meio e FIM do arquivo:
# o ultimo contig e NC_012920.1). Um .bai truncado ou deslocado falha aqui.
valid_bai() {   # $1=bam  $2=bai
    local bam="$1" bai="$2" nsq nidx c a b
    [ -s "$bai" ] || return 1
    [ "$(stat -c %s "$bai")" -ge 500000 ] || { echo "  indice pequeno demais: $(stat -c %s "$bai") B" >&2; return 1; }
    nsq=$(samtools view -H "$bam" | grep -c '^@SQ') || return 1
    nidx=$(samtools idxstats -X "$bam" "$bai" | wc -l) || return 1
    [ "$nidx" -eq $((nsq + 1)) ] || { echo "  idxstats: $nidx linhas, esperado $((nsq + 1))" >&2; return 1; }
    for c in NC_000021.9 NC_000022.11 NC_012920.1; do
        a=$(samtools view -c -X "$bam" "$bai" "$c") || return 1
        b=$(samtools idxstats -X "$bam" "$bai" | awk -v c="$c" '$1==c{print $3+$4}')
        [ -n "$b" ] && [ "$a" = "$b" ] || { echo "  contig $c: view=$a idxstats=$b" >&2; return 1; }
    done
    return 0
}

evict() {   # descarta o arquivo do page cache deste no (somente paginas limpas)
    python3 - "$1" <<'PY'
import os, sys
fd = os.open(sys.argv[1], os.O_RDONLY)
os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
os.close(fd)
PY
}

now() { date +%s.%N; }
SIZE=$(stat -c %s "$BAM")
ACTION="indexed"; STATUS="ok"

if [ -n "${WAIT_FILE:-}" ]; then
    for _ in $(seq 1 120); do [ -e "$WAIT_FILE" ] && break; sleep 5; done
    [ -e "$WAIT_FILE" ] || echo "AVISO: $WAIT_FILE nao apareceu em 10 min; seguindo sem" >&2
fi

echo "[$(date -Is)] $HOST task=${SLURM_ARRAY_TASK_ID} bam=$BAM size=$SIZE threads=$THREADS phase=$PHASE"

# --- 3. idempotencia -------------------------------------------------------------------
T0=$(now)
if [ -f "$BAI" ] && [ "$BAI" -nt "$BAM" ] && valid_bai "$BAM" "$BAI"; then
    ACTION="skipped"
    echo "  .bai ja existe e e valido: pulado"
else
    # BAM ainda sendo escrito? quickcheck exige o marcador EOF.
    samtools quickcheck -v "$BAM" || { echo "ERRO: quickcheck falhou (BAM incompleto/corrompido): $BAM" >&2; STATUS="quickcheck_falhou"; }
    if [ "$STATUS" = "ok" ] && [ -n "$(find "$BAM" -mmin -5 2>/dev/null)" ]; then
        echo "ERRO: BAM modificado ha menos de 5 min (ainda em escrita?): $BAM" >&2; STATUS="bam_recente"
    fi
    if [ "$STATUS" = "ok" ]; then
        [ "$EVICT" = "1" ] && evict "$BAM"
        T0=$(now)                                   # o relogio comeca na leitura
        rm -f "$TMP"
        if samtools index -b -@ "$THREADS" "$BAM" "$TMP"; then
            if valid_bai "$BAM" "$TMP"; then
                mv -f "$TMP" "$BAI" || STATUS="mv_falhou"
            else
                echo "ERRO: indice gerado nao passou na verificacao: $BAM" >&2; STATUS="invalido"
            fi
        else
            echo "ERRO: samtools index saiu com $?" >&2; STATUS="index_falhou"
        fi
    fi
fi
T1=$(now)

# --- 4. registro de vazao (um arquivo por tarefa: append concorrente em NFS nao e atomico)
awk -v ph="$PHASE" -v jid="${SLURM_ARRAY_JOB_ID:-0}" -v tid="${SLURM_ARRAY_TASK_ID}" -v host="$HOST" \
    -v bam="$BAM" -v size="$SIZE" -v t0="$T0" -v t1="$T1" -v act="$ACTION" -v st="$STATUS" -v th="$THREADS" \
    'BEGIN{d=t1-t0; printf "%s\t%s\t%s\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.1f\t%s\t%s\t%s\n", ph,jid,tid,host,bam,size,t0,t1,d,(d>0?size/d/1e6:0),act,st,th}' > "$REC"

[ "$STATUS" = "ok" ] || exit 1
echo "[$(date -Is)] concluido ($ACTION)"
exit 0
