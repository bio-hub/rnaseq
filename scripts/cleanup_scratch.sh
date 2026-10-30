#!/usr/bin/env bash
#
# Apaga o work_dir no scratch — SO depois de confirmar que os entregaveis
# chegaram ao final_dir.
#
#     bash scripts/cleanup_scratch.sh <work_dir> <final_dir> <keep_bams>
#
# Chamado pelo handler onsuccess do Snakefile quando
# config["output"]["cleanup_scratch"] e true.
#
# Esta e uma remocao irreversivel de centenas de GB. A regra e: na duvida,
# NAO apaga. Qualquer verificacao que nao feche preserva o work_dir e sai 0,
# porque falhar aqui nao deve marcar como falha um workflow que deu certo.

set -uo pipefail

WORK="${1:?uso: cleanup_scratch.sh <work_dir> <final_dir> <keep_bams>}"
FINAL="${2:?falta final_dir}"
KEEP="${3:-coord}"

abortar() {
    echo "[cleanup] ABORTADO: $1"
    echo "[cleanup] work_dir PRESERVADO: $WORK"
    exit 0
}

echo
echo "=============================================="
echo "[cleanup] work_dir : $WORK"
echo "[cleanup] final_dir: $FINAL"
echo "=============================================="

#------------------------------------------------------------------------------
# Guardas contra apagar a coisa errada
#------------------------------------------------------------------------------
[ -n "$WORK" ] && [ "$WORK" != "/" ]        || abortar "work_dir vazio ou raiz"
[ -d "$WORK" ]                               || abortar "work_dir nao existe (ja limpo?)"
[ "$WORK" != "$FINAL" ]                      || abortar "work_dir e final_dir sao o mesmo"
case "$WORK" in
    /mnt/scratch-shared/*) ;;
    *) abortar "work_dir fora de /mnt/scratch-shared — recusando apagar" ;;
esac
# final_dir nao pode estar DENTRO do work_dir, senao apagariamos o resultado
case "$FINAL" in
    "$WORK"/*) abortar "final_dir esta dentro do work_dir" ;;
esac

#------------------------------------------------------------------------------
# A sentinela do export e a evidencia de que a copia terminou e foi verificada
#------------------------------------------------------------------------------
[ -f "$FINAL/.export_done" ] || abortar "sentinela $FINAL/.export_done ausente"

#------------------------------------------------------------------------------
# Conferencia independente: nao confia so na sentinela
#------------------------------------------------------------------------------
n_orig=$(find "$WORK/counts" -name "*.genes.results" -type f 2>/dev/null | wc -l)
n_dest=$(find "$FINAL/counts" -name "*.genes.results" -type f 2>/dev/null | wc -l)
echo "[cleanup] contagens RSEM: origem=$n_orig destino=$n_dest"
[ "$n_dest" -gt 0 ]          || abortar "nenhuma contagem no destino"
[ "$n_dest" -eq "$n_orig" ]  || abortar "contagem de arquivos diverge"

if [ "$KEEP" != "none" ]; then
    n_bam_o=$(find "$WORK/mapped" -name "*_rg.bam" -type f 2>/dev/null | wc -l)
    n_bam_d=$(find "$FINAL/mapped" -name "*_Aligned.sortedByCoord.out.bam" -type f 2>/dev/null | wc -l)
    echo "[cleanup] BAMs: origem=$n_bam_o destino=$n_bam_d"
    [ "$n_bam_d" -eq "$n_bam_o" ] || abortar "contagem de BAMs diverge"
fi

vazios=$(find "$FINAL" -type f -size 0 ! -name ".export_done" 2>/dev/null | wc -l)
echo "[cleanup] arquivos vazios no destino: $vazios"
[ "$vazios" -eq 0 ] || abortar "ha arquivo de zero bytes no destino"

[ -f "$FINAL/logs/multiqc_report.html" ] || abortar "multiqc_report.html ausente no destino"

#------------------------------------------------------------------------------
# Tudo conferido: apaga
#------------------------------------------------------------------------------
tam=$(du -sh "$WORK" 2>/dev/null | cut -f1)
echo "[cleanup] todas as verificacoes passaram; removendo $tam de $WORK"
rm -rf "$WORK"

if [ -d "$WORK" ]; then
    echo "[cleanup] AVISO: $WORK ainda existe apos rm -rf"
else
    echo "[cleanup] liberado: $tam"
fi
df -h /mnt/scratch-shared 2>/dev/null | tail -1 | sed 's/^/[cleanup] /'
