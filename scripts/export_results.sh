#!/usr/bin/env bash
#
# Copia os entregaveis do work_dir (scratch) para o final_dir.
#
#     bash scripts/export_results.sh <work_dir> <final_dir> <keep_bams>
#
# keep_bams: none | coord | all
#
# Chamado pela rule export_results. Copia, verifica por CONTEUDO e so entao
# retorna 0 — a sentinela .export_done depende disso.

set -uo pipefail

WORK="${1:?uso: export_results.sh <work_dir> <final_dir> <keep_bams>}"
FINAL="${2:?falta final_dir}"
KEEP="${3:-coord}"

falhou=0

echo "=============================================="
echo "export: $WORK  ->  $FINAL"
echo "keep_bams=$KEEP"
echo "inicio: $(date -Is)"
echo "=============================================="

mkdir -p "$FINAL"/{counts,logs}

copiar() {
    local origem="$1" destino="$2" rotulo="$3"
    local n
    n=$(find "$origem" -maxdepth 1 -name "$4" -type f 2>/dev/null | wc -l)
    if [ "$n" -eq 0 ]; then
        echo "  AVISO: nenhum arquivo '$4' em $origem ($rotulo)"
        return 0
    fi
    mkdir -p "$destino"
    # -u: nao recopia o que ja esta igual, deixa o export idempotente
    find "$origem" -maxdepth 1 -name "$4" -type f -exec cp -u -t "$destino" {} +
    echo "  $rotulo: $n arquivo(s)"
}

echo
echo "[1/4] contagens"
copiar "$WORK/counts" "$FINAL/counts" "contagens STAR"  "*_ReadsPerGene.out.tab"
copiar "$WORK/counts" "$FINAL/counts" "contagens RSEM"  "*.genes.results"
copiar "$WORK/counts" "$FINAL/counts" "isoformas RSEM"  "*.isoforms.results"

# _SJ.out.tab: juncoes de splicing detectadas pelo STAR.
# Nao serve ao DEG. Serve ao STAR 2-pass e a qualquer reanalise que queira
# juncoes novas. Sao ~12 MB por amostra (854 MB em 73) -- barato demais para
# justificar perder junto com o scratch. Exportado por PRECAUCAO: a decisao
# sobre 2-pass fica reversivel.
copiar "$WORK/counts" "$FINAL/counts" "juncoes STAR"    "*_SJ.out.tab"

# {amostra}.stat/ e DIRETORIO, entao `copiar` (que so pega -type f) nao serve.
# Guarda o modelo do RSEM: .cnt, .model, .theta. Util para auditar a
# quantificacao depois, ~5 MB por amostra.
n_stat=0
for d in "$WORK"/counts/*.stat; do
    [ -d "$d" ] || continue
    if cp -ru "$d" "$FINAL/counts/" 2>/dev/null; then
        n_stat=$((n_stat+1))
    else
        echo "  FALHA ao copiar $d"; falhou=1
    fi
done
echo "  modelo RSEM (.stat): $n_stat diretorio(s)"

echo
echo "[2/4] relatorios e logs"
copiar "$WORK/logs" "$FINAL/logs" "multiqc"        "multiqc_report.html"
copiar "$WORK/logs" "$FINAL/logs" "STAR final"     "*_Log.final.out"
copiar "$WORK/logs" "$FINAL/logs" "fastp json"     "*.json"
copiar "$WORK/logs" "$FINAL/logs" "RSEM cnt"       "*.cnt"
if [ -d "$WORK/logs/multiqc_report_data" ]; then
    cp -ru "$WORK/logs/multiqc_report_data" "$FINAL/logs/" && echo "  multiqc_data: ok"
fi

echo
echo "[3/4] BAMs (keep_bams=$KEEP)"
case "$KEEP" in
    none)
        echo "  nenhum BAM exportado"
        ;;
    coord|all)
        # *_rg.bam: o BAM coordenado com @RG utilizavel pelo variant calling.
        # Renomeia para o nome convencional no destino.
        mkdir -p "$FINAL/mapped"
        n=0
        for b in "$WORK"/mapped/*_rg.bam; do
            [ -e "$b" ] || continue
            s=$(basename "$b" _rg.bam)
            d="$FINAL/mapped/${s}_Aligned.sortedByCoord.out.bam"
            if [ ! -f "$d" ] || [ "$b" -nt "$d" ]; then
                cp -u "$b" "$d" || { echo "  FALHA ao copiar $b"; falhou=1; continue; }
            fi
            n=$((n+1))
        done
        echo "  BAM coordenado: $n"
        if [ "$KEEP" = "all" ]; then
            copiar "$WORK/mapped" "$FINAL/mapped" "BAM transcriptoma" "*_Aligned.toTranscriptome.out.bam"
        fi
        ;;
    *)
        echo "  keep_bams invalido: $KEEP"; falhou=1
        ;;
esac

echo
echo "[4/4] DEG, se houver"
if [ -f "$WORK/rdata/DEG.RData" ]; then
    mkdir -p "$FINAL/rdata"
    cp -u "$WORK/rdata/DEG.RData" "$FINAL/rdata/" && echo "  DEG.RData: ok"
else
    echo "  (DEG desligado ou ausente)"
fi

#------------------------------------------------------------------------------
# Verificacao por CONTEUDO.
# Depois do episodio dos 1867 arquivos de zero bytes (queda de energia,
# 2026-10-29), "o arquivo existe" nao e evidencia de nada.
#------------------------------------------------------------------------------
echo
echo "verificacao do destino"
vazios=$(find "$FINAL" -type f -size 0 ! -name ".export_done" 2>/dev/null | wc -l)
if [ "$vazios" -gt 0 ]; then
    echo "  $vazios arquivo(s) de ZERO bytes em $FINAL:"
    find "$FINAL" -type f -size 0 ! -name ".export_done" | head -10 | sed 's/^/    /'
    falhou=1
else
    echo "  nenhum arquivo vazio"
fi

n_counts=$(find "$FINAL/counts" -name "*.genes.results" -type f 2>/dev/null | wc -l)
n_orig=$(find "$WORK/counts" -name "*.genes.results" -type f 2>/dev/null | wc -l)
echo "  contagens RSEM: $n_counts no destino, $n_orig na origem"
[ "$n_counts" -eq "$n_orig" ] || { echo "  DIVERGE"; falhou=1; }

n_sj=$(find "$FINAL/counts" -name "*_SJ.out.tab" -type f 2>/dev/null | wc -l)
n_sj_orig=$(find "$WORK/counts" -name "*_SJ.out.tab" -type f 2>/dev/null | wc -l)
echo "  juncoes STAR:   $n_sj no destino, $n_sj_orig na origem"
[ "$n_sj" -eq "$n_sj_orig" ] || { echo "  DIVERGE"; falhou=1; }

if [ "$KEEP" != "none" ]; then
    n_bam=$(find "$FINAL/mapped" -name "*_Aligned.sortedByCoord.out.bam" -type f 2>/dev/null | wc -l)
    echo "  BAMs no destino: $n_bam"

    # O quickcheck so vale se o samtools existir. Sem esta guarda, um samtools
    # fora do PATH sai com 127 e o `||` abaixo o chamava de BAM corrompido --
    # 5 falsos positivos observados em 2026-11-01, rodando o script fora do
    # ambiente conda. O inverso e pior: zero BAMs amostrados dava "ok".
    if ! command -v samtools >/dev/null 2>&1; then
        echo "  AVISO: samtools fora do PATH -- quickcheck NAO executado."
        echo "         Rode dentro do ambiente conda (biohub_rnaseq) para validar os BAMs."
        falhou=1
    elif [ "$n_bam" -eq 0 ]; then
        echo "  AVISO: nenhum BAM no destino para verificar"
        falhou=1
    else
        # amostra: BAM truncado passa em 'existe' e falha no quickcheck
        ruins=0
        n_check=0
        while read -r b; do
            n_check=$((n_check+1))
            samtools quickcheck "$b" || { echo "    CORROMPIDO: $b"; ruins=$((ruins+1)); }
        done < <(find "$FINAL/mapped" -name "*.bam" -type f | shuf -n 5)
        if [ "$ruins" -eq 0 ] && [ "$n_check" -gt 0 ]; then
            echo "  quickcheck amostral: ok ($n_check de $n_bam)"
        else
            falhou=1
        fi
    fi
fi

sync

echo
if [ "$falhou" -eq 0 ]; then
    echo "EXPORT OK"
    du -sh "$FINAL" 2>/dev/null | sed 's/^/  /'
    echo "fim: $(date -Is)"
    exit 0
else
    echo "EXPORT FALHOU — o work_dir NAO deve ser apagado"
    exit 1
fi
