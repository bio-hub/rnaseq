#!/usr/bin/env bash
#
# Gera o .dict (e o .fai) da referencia NCBI GRCh38.p14 usada no mapeamento STAR.
# O diretorio da referencia (/mnt/reference/rnaseq/GRCh38_p14/) e somente-leitura de proposito
# e NAO e tocado. O GATK exige FASTA + .fai + .dict no MESMO diretorio, entao montamos um
# diretorio duravel ao lado dos dados da tese:
#
#   /mnt/shared/projetos/tese_natalia/referencia_ncbi/
#       GRCh38_p14.fa       -> symlink para o FASTA original (nada e copiado: 3,3 GB)
#       GRCh38_p14.fa.fai
#       GRCh38_p14.dict
#
# Por que ali: e volume duravel e compartilhado (mpatha, visivel em todos os nos), fica junto
# dos BAMs que dependem dele, e nao e o scratchpad da sessao (tmpfs, some no reboot) nem
# /tmp dos nos. O symlink garante que o FASTA e UM so: se a referencia mudar, o .fai e o .dict
# deixam de bater e a validacao abaixo falha em vez de passar em silencio.
#
#   sbatch misc/make-ncbi-dict.sbatch
#
# Idempotente: se .fai e .dict ja existem e passam na validacao, nao faz nada.

#SBATCH --job-name=ncbi-dict
#SBATCH --partition=gpu,cpu,bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%j.out

set -uo pipefail
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH

SRC=/mnt/reference/rnaseq/GRCh38_p14/fasta/GRCh38_p14.fa
DEST=/mnt/shared/projetos/tese_natalia/referencia_ncbi
FAI_REF="${FAI_REF:-/home/ronald/Documentos/ronald/tese_natalia/evidencia_juncoes/GRCh38_p14.fa.fai}"  # so existe no head; opcional
# M5 independentes (verificados contra o dict do interomics): NC_000001.11 == chr1, NC_012920.1 == chrM
M5_CHR1=6aef897c3d6ff0c78aff06ac189178dd
M5_CHRM=c68f52674c9fb33aef52dcf399755519

GATK="${GATK:-}"
if [ -z "$GATK" ]; then
    for g in /mnt/shared/conda/snakemake_envs/interomics/*/bin/gatk; do [ -x "$g" ] && { GATK="$g"; break; }; done
fi
command -v samtools >/dev/null 2>&1 || { echo "ERRO: samtools ausente em $(hostname -s)" >&2; exit 127; }
[ -n "$GATK" ] && [ -x "$GATK" ] || { echo "ERRO: gatk nao encontrado (defina GATK=...)" >&2; exit 127; }
# O wrapper 'gatk' chama 'java' do PATH: num job SLURM o ambiente conda dele NAO esta ativado
# (falhou aqui com FileNotFoundError: 'java'). Entao o bin do ambiente do gatk entra no PATH.
export PATH="$(dirname "$GATK"):$PATH"
command -v java >/dev/null 2>&1 || { echo "ERRO: java ausente no PATH apos incluir $(dirname "$GATK")" >&2; exit 127; }
[ -f "$SRC" ] || { echo "ERRO: FASTA nao existe: $SRC" >&2; exit 2; }

mkdir -p "$DEST" || exit 1
FA="$DEST/GRCh38_p14.fa"; FAI="$FA.fai"; DICT="$DEST/GRCh38_p14.dict"
[ -L "$FA" ] || ln -s "$SRC" "$FA"
[ "$(readlink -f "$FA")" = "$(readlink -f "$SRC")" ] || { echo "ERRO: $FA nao aponta para $SRC" >&2; exit 1; }

valid() {
    [ -s "$FAI" ] && [ -s "$DICT" ] || return 1
    local n_fai n_dict m1 mm
    n_fai=$(wc -l < "$FAI")
    n_dict=$(grep -c '^@SQ' "$DICT")
    [ "$n_fai" = "$n_dict" ] || { echo "  fai=$n_fai contigs, dict=$n_dict" >&2; return 1; }
    # nome e comprimento batem entre .fai e .dict, na mesma ordem
    paste <(cut -f1,2 "$FAI") <(grep '^@SQ' "$DICT" | awk -F'\t' '{sub("SN:","",$2); sub("LN:","",$3); print $2"\t"$3}') \
        | awk -F'\t' '$1!=$3 || $2!=$4{bad++} END{exit bad?1:0}' || { echo "  nomes/comprimentos divergem entre .fai e .dict" >&2; return 1; }
    m1=$(grep -P '\tSN:NC_000001.11\t' "$DICT" | grep -o 'M5:[0-9a-f]*' | cut -d: -f2)
    mm=$(grep -P '\tSN:NC_012920.1\t'  "$DICT" | grep -o 'M5:[0-9a-f]*' | cut -d: -f2)
    [ "$m1" = "$M5_CHR1" ] && [ "$mm" = "$M5_CHRM" ] || { echo "  M5 inesperado: chr1=$m1 chrM=$mm" >&2; return 1; }
    return 0
}

if valid; then echo "ja existe e e valido: $DICT"; exit 0; fi

echo "[$(date -Is)] gerando .fai e .dict em $DEST"
T0=$(date +%s)
rm -f "$FAI" "$DICT"
samtools faidx "$FA" || { echo "ERRO: faidx falhou" >&2; exit 1; }
# A saida precisa terminar em .dict (o GATK nao decide formato por extensao, mas o nome e o contrato)
"$GATK" --java-options "-Xmx4g" CreateSequenceDictionary -R "$FA" -O "$DICT" || { echo "ERRO: CreateSequenceDictionary falhou" >&2; rm -f "$DICT"; exit 1; }
T1=$(date +%s)
echo "tempo: $((T1 - T0)) s para $(stat -L -c %s "$FA") bytes de FASTA"

# conferencia extra: o .fai gerado aqui e igual ao gerado antes em evidencia_juncoes (se visivel neste no)
if [ -s "$FAI_REF" ]; then
    cmp -s "$FAI" "$FAI_REF" && echo ".fai identico ao de evidencia_juncoes" || echo "AVISO: .fai difere do de evidencia_juncoes" >&2
fi
valid || { echo "ERRO: a validacao final falhou" >&2; exit 1; }
echo "OK: $(grep -c '^@SQ' "$DICT") contigs; M5 de NC_000001.11 e NC_012920.1 conferem"
ls -la "$DEST"
