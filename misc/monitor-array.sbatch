#!/usr/bin/env bash
#
# Monitor do array de armazenamento durante uma medicao de vazao.
# Roda no r660xs (que e o servidor NFS: o volume /mnt/shared e o /dev/mapper/mpatha local dele),
# amostra a cada INTERVAL s os contadores do dispositivo e da NIC de dados, e para quando
# existir $LOGDIR/STOP (ou ao fim do --time).
#
#   sbatch -w r660xs misc/monitor-array.sbatch
#   ... (rodar a medicao) ...
#   touch /mnt/scratch-shared/ronald/bai_piloto/STOP
#
# Antes de amostrar, descarta do page cache do r660xs os BAMs listados em EVICT_LIST
# (um por linha): sem isso, um BAM lido recentemente sai da RAM do servidor e a medicao
# mede a memoria, nao o array. Quando termina o descarte, cria $LOGDIR/evict_done.
#
# Saida: $LOGDIR/monitor/monitor_<jobid>.tsv com colunas
#   epoch  rd_sectors  wr_sectors  in_flight  io_ticks_ms  nic_rx_bytes  nic_tx_bytes  loadavg1
# (setores de 512 B; io_ticks_ms/intervalo = %util do dispositivo)

#SBATCH --job-name=bai-monitor
#SBATCH --partition=bigmem
#SBATCH --nodelist=r660xs
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%j.out

set -uo pipefail
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH
LOGDIR="${LOGDIR:-/mnt/scratch-shared/ronald/bai_piloto}"
INTERVAL="${INTERVAL:-5}"
OUT="$LOGDIR/monitor/monitor_${SLURM_JOB_ID:-manual}.tsv"
mkdir -p "$LOGDIR/monitor"
rm -f "$LOGDIR/STOP" "$LOGDIR/evict_done"

for t in python3 awk ip readlink; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERRO: '$t' ausente em $(hostname -s)" >&2; exit 127; }
done
[ "$(hostname -s)" = "r660xs" ] || { echo "ERRO: este monitor so faz sentido no r660xs (servidor do array); estou em $(hostname -s)" >&2; exit 2; }

DEV="$(basename "$(readlink -f /dev/mapper/mpatha)")"             # dm-2
NIC="$(ip -o -4 addr show | awk '$4 ~ /^192\.168\.100\./{print $2; exit}')"
[ -n "$DEV" ] && [ -n "$NIC" ] || { echo "ERRO: nao achei dispositivo ($DEV) ou NIC ($NIC)" >&2; exit 2; }
echo "monitorando dispositivo=$DEV nic=$NIC intervalo=${INTERVAL}s -> $OUT"

if [ -n "${EVICT_LIST:-}" ] && [ -s "$EVICT_LIST" ]; then
    python3 - "$EVICT_LIST" <<'PY'
import os, sys
n = 0
for line in open(sys.argv[1]):
    p = line.strip()
    if not p: continue
    fd = os.open(p, os.O_RDONLY)
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    os.close(fd); n += 1
print("descartados do page cache do servidor:", n, "arquivos")
PY
fi
touch "$LOGDIR/evict_done"

printf 'epoch\trd_sectors\twr_sectors\tin_flight\tio_ticks_ms\tnic_rx_bytes\tnic_tx_bytes\tloadavg1\n' > "$OUT"
while [ ! -e "$LOGDIR/STOP" ]; do
    d=$(awk -v dev="$DEV" '$3==dev{print $6"\t"$10"\t"$12"\t"$13}' /proc/diskstats)
    n=$(awk -v nic="$NIC:" '$1==nic{print $2"\t"$10}' /proc/net/dev)
    l=$(cut -d' ' -f1 /proc/loadavg)
    printf '%s\t%s\t%s\t%s\n' "$(date +%s)" "$d" "$n" "$l" >> "$OUT"
    sleep "$INTERVAL"
done
echo "STOP detectado; fim."
