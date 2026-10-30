#!/usr/bin/env bash
#
# Repara o @RG de BAMs produzidos por versoes do workflow anteriores ao commit
# que corrigiu --outSAMattrRGline.
#
#     bash misc/repair-bam-readgroups.sh <dir_bams> [modo] [paralelismo]
#
#     modo: check (padrao) | repair
#
# Exemplos:
#     bash misc/repair-bam-readgroups.sh ~/rnaseq/output/mapped check
#     bash misc/repair-bam-readgroups.sh ~/rnaseq/output/mapped repair 4
#
#------------------------------------------------------------------------------
# O problema
#------------------------------------------------------------------------------
# O STAR foi chamado com --outSAMattrRGline ID:X\tSM:X\tLB:...  O `\t` virou o
# literal 't' (o bash converte \t fora de aspas), colapsando os quatro campos
# num unico:
#
#     @RG  ID:SRR25556754tSM:SRR25556754tLB:TRUSEQtPL:ILLUMINA
#
# As contagens (RSEM, ReadsPerGene) NAO sao afetadas — elas ignoram read group.
# Mas MarkDuplicates e BQSR do interomics_snakemake exigem SM: e PL: como
# campos proprios. Sem eles o pipeline de variant calling nao roda.
#
# Reescrever o header custa uma passada de I/O por BAM; realinhar custaria
# horas por amostra. Por isso reparar, nao remapear.
#
#------------------------------------------------------------------------------
# Seguranca
#------------------------------------------------------------------------------
# - modo `check` (padrao) nao escreve nada; so relata
# - escreve em .tmp e so entao faz mv: o BAM original nunca fica num estado
#   intermediario. Uma interrupcao deixa um .tmp orfao, nao um BAM truncado
# - idempotente: BAM com @RG ja valido e pulado
# - verifica o resultado (header + quickcheck) ANTES de substituir o original

set -uo pipefail

DIR="${1:?uso: repair-bam-readgroups.sh <dir_bams> [check|repair] [paralelismo]}"
MODO="${2:-check}"
PAR="${3:-4}"

case "$MODO" in check|repair) ;; *) echo "modo invalido: $MODO (use check ou repair)"; exit 2 ;; esac
[ -d "$DIR" ] || { echo "diretorio nao existe: $DIR"; exit 2; }
command -v samtools >/dev/null || { echo "samtools nao encontrado no PATH"; exit 2; }

PADRAO="*_Aligned.sortedByCoord.out.bam"

echo "=============================================="
echo "reparo de @RG"
echo "  dir         : $DIR"
echo "  modo        : $MODO"
echo "  paralelismo : $PAR"
echo "  inicio      : $(date -Is)"
echo "=============================================="

#------------------------------------------------------------------------------
# GUARDA 1: nao reparar com o workflow rodando.
#
# Um BAM que o STAR esta escrevendo ainda nao tem @RG e aparece como "precisa
# de reparo". Reescreve-lo corromperia o alinhamento em curso. Verificado em
# 2026-10-30: 6 BAMs de zero bytes, cujos nomes batiam com os STAR ativos.
#------------------------------------------------------------------------------
if [ "$MODO" = "repair" ]; then
    bloqueio=""

    # (a) BAM escrito nos ultimos 5 min NESTE diretorio. E o sinal mais direto
    #     de escrita em curso, e nao depende de adivinhar nomes de processo.
    recentes=$(find "$DIR" -maxdepth 1 -name "$PADRAO" -mmin -5 2>/dev/null | wc -l)
    [ "$recentes" -eq 0 ] || bloqueio="$recentes BAM(s) modificados nos ultimos 5 min"

    # (b) STAR alinhando. O build de indice (reference.smk) tambem usa STAR mas
    #     nao escreve BAM, entao e excluido para nao bloquear o reparo a toa.
    star=$(pgrep -f "STAR .*--readFilesIn" 2>/dev/null | wc -l)
    [ "$star" -eq 0 ] || bloqueio="${bloqueio:+$bloqueio; }$star processo(s) STAR alinhando"

    if [ -n "$bloqueio" ]; then
        echo
        echo "ABORTANDO: $bloqueio"
        echo
        echo "Rode o reparo SO depois que o workflow terminar. Um BAM que o STAR"
        echo "esta escrevendo ainda nao tem @RG e apareceria como 'precisa de"
        echo "reparo' — reescreve-lo destruiria o alinhamento em curso."
        echo "Verificado em 2026-10-30: 6 BAMs de zero bytes cujos nomes batiam"
        echo "com os STAR ativos."
        echo
        pgrep -a -f "STAR .*--readFilesIn" 2>/dev/null | cut -c1-120 | head -3 | sed 's/^/  /'
        exit 1
    fi
fi

#------------------------------------------------------------------------------
# Levantamento
#------------------------------------------------------------------------------
mapfile -t BAMS < <(find "$DIR" -maxdepth 1 -name "$PADRAO" -type f | sort)
total=${#BAMS[@]}
[ "$total" -gt 0 ] || { echo "nenhum BAM '$PADRAO' em $DIR"; exit 0; }

echo
echo "[1/3] inspecionando $total BAM(s)"

quebrados=(); validos=(); semrg=(); incompletos=()
for b in "${BAMS[@]}"; do
    # GUARDA 2: BAM sem marcador de EOF esta incompleto (em escrita) ou
    # corrompido. Em nenhum dos casos pode ser reparado — reparar um BAM em
    # escrita o destroi, e um corrompido precisa ser refeito, nao remendado.
    if ! samtools quickcheck "$b" 2>/dev/null; then
        incompletos+=("$b")
        continue
    fi
    rg=$(samtools view -H "$b" 2>/dev/null | grep -m1 '^@RG' || true)
    if [ -z "$rg" ]; then
        semrg+=("$b")
    elif printf '%s' "$rg" | grep -q $'\t''SM:'; then
        validos+=("$b")
    else
        quebrados+=("$b")
    fi
done

echo "  @RG valido (tem SM: proprio) : ${#validos[@]}"
echo "  @RG malformado               : ${#quebrados[@]}"
echo "  sem @RG, mas BAM completo    : ${#semrg[@]}"
echo "  INCOMPLETO/CORROMPIDO        : ${#incompletos[@]}  (pulados)"

if [ ${#incompletos[@]} -gt 0 ]; then
    echo
    echo "  estes nao serao tocados — em escrita agora ou precisam ser refeitos:"
    for b in "${incompletos[@]}"; do
        printf '    %-48s %s  mtime=%smin\n' \
          "$(basename "$b")" "$(du -h "$b" | cut -f1)" \
          "$(( ($(date +%s) - $(stat -c %Y "$b")) / 60 ))"
    done
fi

if [ ${#quebrados[@]} -gt 0 ]; then
    echo
    echo "  exemplo do header quebrado:"
    samtools view -H "${quebrados[0]}" | grep -m1 '^@RG' | sed 's/^/    /'
fi

a_reparar=( "${quebrados[@]}" "${semrg[@]}" )
if [ ${#a_reparar[@]} -eq 0 ]; then
    echo
    echo "Nada a fazer: todos os $total BAMs tem @RG valido."
    exit 0
fi

if [ "$MODO" = "check" ]; then
    echo
    echo "=============================================="
    echo "${#a_reparar[@]} de $total BAM(s) precisam de reparo."
    echo "Para aplicar:"
    echo "    bash $0 $DIR repair $PAR"
    echo "=============================================="
    exit 0
fi

#------------------------------------------------------------------------------
# Espaco: o .tmp coexiste com o original ate o mv
#------------------------------------------------------------------------------
maior=$(for b in "${a_reparar[@]}"; do stat -c %s "$b"; done | sort -n | tail -1)
livre=$(df -B1 --output=avail "$DIR" | tail -1)
margem=$(( maior * PAR * 2 ))
echo
echo "[2/3] espaco: maior BAM $(numfmt --to=iec "$maior"), livre $(numfmt --to=iec "$livre")"
if [ "$livre" -lt "$margem" ]; then
    echo "  ABORTANDO: precisa de ~$(numfmt --to=iec "$margem") para $PAR em paralelo"
    echo "  reduza o paralelismo ou libere espaco"
    exit 1
fi
echo "  ok"

#------------------------------------------------------------------------------
# Reparo
#------------------------------------------------------------------------------
echo
echo "[3/3] reparando ${#a_reparar[@]} BAM(s), $PAR em paralelo"
echo

reparar_um() {
    local b="$1"
    local s tmp rg_novo
    s=$(basename "$b" _Aligned.sortedByCoord.out.bam)
    tmp="${b}.rgfix.tmp"

    rm -f "$tmp"
    # -w sobrescreve a linha @RG existente (que esta quebrada)
    if ! samtools addreplacerg \
            -r "ID:${s}" -r "SM:${s}" -r "LB:TRUSEQ" -r "PL:ILLUMINA" \
            -w -o "$tmp" -@ 2 "$b" 2>/dev/null; then
        echo "FALHA  $s  (addreplacerg)"
        rm -f "$tmp"
        return 1
    fi

    # Verifica ANTES de tocar no original
    if ! samtools quickcheck "$tmp" 2>/dev/null; then
        echo "FALHA  $s  (quickcheck do temporario)"
        rm -f "$tmp"
        return 1
    fi
    rg_novo=$(samtools view -H "$tmp" 2>/dev/null | grep -m1 '^@RG' || true)
    if ! printf '%s' "$rg_novo" | grep -q $'\t''SM:'; then
        echo "FALHA  $s  (@RG continua sem SM: proprio)"
        rm -f "$tmp"
        return 1
    fi

    # mv no mesmo filesystem e rename atomico
    mv -f "$tmp" "$b" || { echo "FALHA  $s  (mv)"; rm -f "$tmp"; return 1; }

    # .bai: o STAR nao gera, e o variant calling precisa
    samtools index -@ 2 "$b" 2>/dev/null || echo "AVISO  $s  (index falhou)"

    echo "OK     $s"
    return 0
}
export -f reparar_um

printf '%s\n' "${a_reparar[@]}" \
  | xargs -P "$PAR" -I{} bash -c 'reparar_um "$@"' _ {} \
  | tee /tmp/.rgfix_saida_$$ || true

ok=$(grep -c '^OK' /tmp/.rgfix_saida_$$ 2>/dev/null || echo 0)
falhas=$(grep -c '^FALHA' /tmp/.rgfix_saida_$$ 2>/dev/null || echo 0)
rm -f /tmp/.rgfix_saida_$$

#------------------------------------------------------------------------------
# Conferencia final, independente do que o laco reportou
#------------------------------------------------------------------------------
echo
echo "=============================================="
echo "reparados: $ok   falhas: $falhas"
echo
echo "conferencia final de todos os $total BAM(s):"
ainda=0; indexados=0
for b in "${BAMS[@]}"; do
    rg=$(samtools view -H "$b" 2>/dev/null | grep -m1 '^@RG' || true)
    printf '%s' "$rg" | grep -q $'\t''SM:' || { echo "  AINDA QUEBRADO: $(basename "$b")"; ainda=$((ainda+1)); }
    [ -f "${b}.bai" ] && indexados=$((indexados+1))
done
echo "  com @RG valido : $(( total - ainda ))/$total"
echo "  com .bai       : $indexados/$total"

orfaos=$(find "$DIR" -maxdepth 1 -name "*.rgfix.tmp" | wc -l)
[ "$orfaos" -eq 0 ] || echo "  AVISO: $orfaos arquivo(s) .rgfix.tmp orfaos — remova apos conferir"

sync
echo "fim: $(date -Is)"
echo "=============================================="
[ "$ainda" -eq 0 ] && [ "$falhas" -eq 0 ]
