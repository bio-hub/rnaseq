#!/usr/bin/env bash
#
# Recomprime BAMs que foram escritos sem compressao.
#
#   bash misc/recompress-bams.sh <dir_bams> [check|run] [paralelismo] [limiar_GB]
#
# Por que isto existe: o samtools escolhe o formato de saida pela EXTENSAO do
# arquivo. Uma ferramenta que escreve em "saida.tmp" ou "saida.part" em vez de
# "saida.bam" produz BAM VALIDO e SEM COMPRESSAO — tipicamente ~13x maior.
#
# O arquivo passa em tudo: samtools quickcheck aprova, o header esta correto,
# as reads estao todas la. So o tamanho denuncia, e normalmente isso aparece
# quando o disco acaba.
#
# Seguranca: escreve .recomp.bam, verifica contagem de reads E header contra o
# original, e so entao substitui por rename atomico. Nunca apaga antes de
# confirmar que o novo arquivo tem o mesmo conteudo.

set -uo pipefail

DIR="${1:?uso: recompress-bams.sh <dir_bams> [check|run] [paralelismo] [limiar_GB]}"
MODO="${2:-check}"
PAR="${3:-3}"
LIMIAR_GB="${4:-10}"

case "$MODO" in check|run) ;; *) echo "modo invalido: $MODO"; exit 2 ;; esac
[ -d "$DIR" ] || { echo "diretorio nao existe: $DIR"; exit 2; }
command -v samtools >/dev/null || { echo "samtools nao esta no PATH"; exit 2; }

LIMIAR=$(( LIMIAR_GB * 1000000000 ))

echo "=============================================="
echo "recompressao de BAM"
echo "  dir     : $DIR   modo: $MODO   limiar: ${LIMIAR_GB} GB"
echo "  inicio  : $(date -Is)"
echo "=============================================="

# Mesma guarda do reparo de @RG: nao mexer em arquivo sendo escrito.
if [ "$MODO" = "run" ]; then
    recentes=$(find "$DIR" -maxdepth 1 -name "*.bam" -mmin -5 2>/dev/null | wc -l)
    ativos=$(pgrep -f "STAR .*--readFilesIn" 2>/dev/null | wc -l)
    if [ "$recentes" -gt 0 ] || [ "$ativos" -gt 0 ]; then
        echo "ABORTANDO: $recentes BAM(s) modificados ha <5min, $ativos STAR ativo(s)."
        exit 1
    fi
fi

mapfile -t ALVOS < <(find "$DIR" -maxdepth 1 -name "*.bam" -size +"${LIMIAR_GB}"G -type f | sort)
echo
echo "[1/2] BAMs acima de ${LIMIAR_GB} GB: ${#ALVOS[@]}"
[ "${#ALVOS[@]}" -gt 0 ] || { echo "  nada a fazer"; exit 0; }
for b in "${ALVOS[@]}"; do printf "  %-52s %s\n" "$(basename "$b")" "$(du -h "$b" | cut -f1)"; done

if [ "$MODO" = "check" ]; then
    echo
    echo "Para recomprimir:  bash $0 $DIR run $PAR $LIMIAR_GB"
    exit 0
fi

recomp_um() {
    local b="$1"
    local nome tmp n_orig n_novo
    nome=$(basename "$b")
    tmp="${b%.bam}.recomp.bam"     # extensao .bam: e ela que define a compressao

    n_orig=$(samtools view -c -@ 2 "$b" 2>/dev/null) || { echo "FALHA  $nome (contagem origem)"; return 1; }

    rm -f "$tmp"
    if ! samtools view -b -O bam -@ 3 -o "$tmp" "$b" 2>/dev/null; then
        echo "FALHA  $nome (view)"; rm -f "$tmp"; return 1
    fi
    samtools quickcheck "$tmp" 2>/dev/null || { echo "FALHA  $nome (quickcheck)"; rm -f "$tmp"; return 1; }

    # O conteudo precisa ser identico: mesmo numero de reads e mesmo @RG.
    n_novo=$(samtools view -c -@ 2 "$tmp" 2>/dev/null)
    if [ "$n_orig" != "$n_novo" ]; then
        echo "FALHA  $nome (reads $n_orig -> $n_novo)"; rm -f "$tmp"; return 1
    fi
    if ! diff -q <(samtools view -H "$b" | grep '^@RG') <(samtools view -H "$tmp" | grep '^@RG') >/dev/null 2>&1; then
        echo "FALHA  $nome (@RG divergiu)"; rm -f "$tmp"; return 1
    fi

    local antes depois
    antes=$(stat -c %s "$b"); depois=$(stat -c %s "$tmp")
    mv -f "$tmp" "$b" || { echo "FALHA  $nome (mv)"; rm -f "$tmp"; return 1; }
    samtools index -@ 2 "$b" 2>/dev/null || echo "AVISO  $nome (index)"

    echo "OK     $nome  $(( antes/1024/1024 ))MB -> $(( depois/1024/1024 ))MB"
    return 0
}
export -f recomp_um

echo
echo "[2/2] recomprimindo, $PAR em paralelo"
echo
printf '%s\n' "${ALVOS[@]}" | xargs -P "$PAR" -I{} bash -c 'recomp_um "$@"' _ {}

echo
echo "=============================================="
echo "conferencia final"
rest=$(find "$DIR" -maxdepth 1 -name "*.bam" -size +"${LIMIAR_GB}"G -type f | wc -l)
orf=$(find "$DIR" -maxdepth 1 -name "*.recomp.bam" | wc -l)
echo "  ainda acima de ${LIMIAR_GB} GB : $rest"
echo "  .recomp.bam orfaos          : $orf"
echo "  total do diretorio          : $(du -sh "$DIR" | cut -f1)"
sync
echo "fim: $(date -Is)"
[ "$rest" -eq 0 ] && [ "$orf" -eq 0 ]
