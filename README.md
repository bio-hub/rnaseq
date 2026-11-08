# RNA-seq Analysis Workflow

Workflow Snakemake para RNA-seq paired-end: trimagem (fastp), alinhamento (STAR),
quantificacao (RSEM + `ReadsPerGene` do STAR), QC (MultiQC) e, opcionalmente,
expressao diferencial (DESeq2 + ReactomePA).

Dois fatos de desenho que explicam o resto deste README:

1. A **referencia e um insumo externo**, compartilhado e somente-leitura, construido
   **uma unica vez** por um workflow separado. Nenhuma analise baixa ou reconstroi a referencia.
2. O trabalho acontece em **duas arvores**: intermediarios volumosos vao para o
   `work_dir` (scratch) e so os **entregaveis** sao exportados para o `final_dir`.
   Os BAMs exportados sao a **entrada do variant calling** (ver "Elo com o variant calling").

> **Referencia = NCBI RefSeq GRCh38.p14 (contigs `NC_000001.11`).** O workflow de
> variant calling (`interomics_snakemake`) usa outra (contigs `chr1`). Leia
> [`documentation/REFERENCIA_UNIFICADA.md`](documentation/REFERENCIA_UNIFICADA.md)
> antes de levar um BAM daqui para la.

## Documentacao

| Arquivo | Para que serve |
|---|---|
| [`documentation/RUNBOOK.md`](documentation/RUNBOOK.md) | execucao **local** no `r660xs`: pre-requisitos, o que ja esta pronto, o que nao fazer, verificacao |
| [`documentation/RUNBOOK_SLURM.md`](documentation/RUNBOOK_SLURM.md) | execucao de um lote via **SLURM**: perfil, recursos, distribuicao entre nos, medicoes |
| [`documentation/REFERENCIA_UNIFICADA.md`](documentation/REFERENCIA_UNIFICADA.md) | a divergencia de genoma com o `interomics_snakemake`, o que ela custa e a proposta de unificacao |

## Requisitos

- Ambiente conda de `environment.yaml` (Snakemake 9, STAR 2.7.11b, RSEM 1.3.3, fastp, samtools, MultiQC, R/Bioconductor)
- FASTQ paired-end em `samples/`, nomeados `{sample}_1.fastq.gz` e `{sample}_2.fastq.gz`.
  O nome da amostra e o que sobra ao tirar `_1.fastq.gz`; arquivo fora do padrao e ignorado
  **em silencio** (se `SAMPLES` sair vazio, o dry-run nao mostra job nenhum)
- Single-end nao e suportado por este Snakefile
- A referencia ja construida (ver abaixo)

## Instalacao

```bash
conda env create -f environment.yaml
conda activate biohub_rnaseq
```

Neste cluster existe um ambiente compartilhado em `/mnt/shared/conda` (ver RUNBOOK).

## Passo 1 - construir a referencia (uma vez por genoma e por `read_length`)

**Neste servidor ela ja existe** em `/mnt/reference/rnaseq/GRCh38_p14/`, com indices
`rsem_star_sjdb{74,75,99,124,149}/`. So rode o build para um `read_length` que nao esteja la.

```bash
# o diretorio raiz precisa existir e ser gravavel
sudo mkdir -p /mnt/reference/rnaseq
sudo chown "$USER":"$USER" /mnt/reference/rnaseq

bash scripts/build_reference.sh
```

Baixa genoma, GTF e o `gene_info` do **NCBI** (URLs em `reference.*_url` do `config.yaml`)
e constroi o indice STAR/RSEM em:

```
/mnt/reference/rnaseq/GRCh38_p14/
├── fasta/        gtf/        gene_info/
├── rsem_star_sjdb149/        # o overhang derivado entra no nome
├── logs/         MANIFEST.txt
└── .done_sjdb149             # sentinela: so existe se o build terminou
```

Reserve ~40 GB de disco e algumas horas (cada indice adicional: ~33 GB). O script segura
um `flock` exclusivo: dois builds ao mesmo tempo, o segundo aborta em vez de corromper o
indice. Ao terminar, os arquivos ficam somente-leitura. Para reconstruir de proposito:

```bash
chmod -R u+w /mnt/reference/rnaseq/GRCh38_p14
```

O `build_reference.sh` le `config.yaml` e, se existir, `config.local.yaml` (ou o arquivo
apontado por `SNAKEMAKE_CONFIG`). `SNAKEMAKE_CORES` (20) e `SNAKEMAKE_MEM_MB` (64000) ajustam recursos.

### O que a referencia e (e o que isso implica)

- Genoma **RefSeq GRCh38.p14** completo: 705 contigs, com alt loci e patches, nomes `NC_/NT_/NW_`.
- `gene_id` do GTF e o **simbolo** do gene. E ele que vira linha em `ReadsPerGene.out.tab`,
  e o `DEG.R` junta com o `gene_info` por simbolo para obter o Entrez do Reactome.
- Isso **nao e** o namespace do `interomics_snakemake` (`chr1`, 195 contigs, sem ALT). Os BAMs
  precisam de conformacao (`rna_conform_bam`) antes do variant calling. Detalhes, custo e
  proposta de unificacao: `documentation/REFERENCIA_UNIFICADA.md`.

## Passo 2 - rodar a analise

Local (padrao):

```bash
conda activate biohub_rnaseq
snakemake --workflow-profile profiles/local
```

SLURM:

```bash
snakemake --workflow-profile profiles/slurm
```

Sempre faca antes um dry-run: `snakemake --workflow-profile profiles/local -n`.

### Local vs SLURM

| | `profiles/local` | `profiles/slurm` |
|---|---|---|
| Executor | processos locais | `snakemake-executor-plugin-slurm` |
| Concorrencia | `cores: 120` e `resources.mem_mb: 900000` (dimensionado para o r660xs) | `jobs: 20`, o SLURM gateia por recurso |
| Particao | - | `bigmem` (r660xs). Para usar os R550, `bigmem,cpu`; ver RUNBOOK_SLURM |
| Falha transitoria | `rerun-incomplete` | `rerun-incomplete` + `restart-times: 1` |
| Orquestrador | precisa de `tmux` | precisa de `tmux` ou `sbatch --wrap` |

O `mem_mb` e o que impede o Snakemake de agendar jobs pesados demais em paralelo: sem
ele, varios `map_reads` de ~70 GB sao disparados juntos. Ajuste o profile para a maquina.
Nao ha `groups:` no profile SLURM, de proposito (ver comentario no proprio arquivo).
O workflow nao depende do SLURM: em qualquer maquina, o profile local basta.

## Configuracao

Tudo em `config.yaml`. Os campos que mais importam:

| Campo | Para que serve |
|---|---|
| `reference.dir` | raiz das referencias compartilhadas |
| `reference.genome_id` | subdiretorio do build; um por genoma |
| `reference.tax_id` | especie filtrada no `gene_info` (9606 = humano) |
| `star.read_length` | comprimento das reads em bp. O workflow deriva `sjdbOverhang = read_length - 1` |
| `star.twopass` / `star.shared_genome` | trade-off de mapeamento; **mutuamente exclusivos** (ver abaixo) |
| `output.project`, `work_dir`, `final_dir` | onde escrever (ver "Saida") |
| `output.keep_bams`, `cleanup_scratch`, `local_scratch` | o que exportar, o que apagar, onde ficam temporarios |
| `fix_read_group` | rede de seguranca do `@RG` (ver "Read group") |
| `threads.*` / `mem_mb.*` / `runtime_min.*` | recursos por rule (`runtime_min` so vale no SLURM) |
| `deg_analysis.*` | liga/desliga a etapa, design experimental, alpha e limiar de log2FC |

Para rodar em outra maquina, caminho ou lote, **nao edite o `config.yaml` versionado**:
crie um `config.local.yaml` (ja no `.gitignore`) so com o que muda, e passe os dois:

```yaml
# config.local.yaml
reference:
  dir: "/caminho/alternativo"
output:
  project: "hs_rnaseq_L5"
```

```bash
snakemake --configfile config.yaml config.local.yaml --workflow-profile profiles/local
```

### Comprimento das reads: `star.read_length`

Informe `star.read_length` em bp; o workflow deriva o `sjdbOverhang` que o STAR espera
(`read_length - 1`) e o usa tanto no build do indice quanto no mapeamento, que precisam
bater. O valor entra no nome do indice e da sentinela, entao projetos de comprimentos
diferentes convivem no mesmo volume:

| `read_length` | indice | sentinela |
|---|---|---|
| 150 | `rsem_star_sjdb149/` | `.done_sjdb149` |
| 125 | `rsem_star_sjdb124/` | `.done_sjdb124` |
| 100 | `rsem_star_sjdb99/` | `.done_sjdb99` |
| 76 | `rsem_star_sjdb75/` | `.done_sjdb75` |
| 75 | `rsem_star_sjdb74/` | `.done_sjdb74` |

(sao os cinco indices existentes hoje). Para um comprimento novo, **nao edite o
`config.yaml` versionado**:

```bash
printf 'star:\n  read_length: 100\n' > config.local.yaml
bash scripts/build_reference.sh        # reaproveita FASTA e GTF ja baixados
snakemake --configfile config.yaml config.local.yaml --workflow-profile profiles/local
```

O Snakemake **nao detecta** o comprimento das reads a partir dos FASTQ: ele confia no
config e exige que o indice correspondente exista. Se nao existir, falha no parse com
mensagem explicita, em vez de deixar o STAR abortar no meio do mapeamento.

`star.sjdbOverhang` continua aceito como override opcional, **com precedencia** sobre
`read_length` (e o workflow avisa em stderr). Serve para reads de comprimento misto (use o
maior - 1) ou para reaproveitar um indice construido com outro valor.

### 2-pass ou genoma compartilhado: escolha um

| | `star.twopass: true` | `star.shared_genome: true` |
|---|---|---|
| O que faz | passada 1 descobre juncoes novas, reinseridas antes da passada 2 | carrega o genoma uma vez em memoria compartilhada (`--genomeLoad`) |
| Custo | dobra o `map_reads` (~11 min para ~20-23 min na amostra mediana) | segmento de shm **vaza** (~32 GB) se um job morre; um genoma por no; cada no precisa do proprio load |
| Ganho medido (L1/L3) | so 0,8% a 1,7% das reads unicas apoiam juncoes nao anotadas | a carga do genoma e ~24 s de mediana contra ~653 s do mapeamento, ~3,5% |
| Quando usar | variant calling (GATK Best Practices pedem 2-pass) | execucao local, nunca em fila SLURM |

**Sao mutuamente exclusivos:** o STAR aborta com "2-pass method is not compatible with
--genomeLoad". `common.smk` falha no parse se os dois estiverem ligados. Ambos ficam
`false` por padrao. Ligar `twopass` **remapeia tudo e muda as contagens**: para DEG puro, deixe `false`.
O `_SJ.out.tab` de cada amostra e exportado de qualquer forma, o que mantem a decisao reversivel.

### Read group e MAPQ (o que o STAR escreve e por que)

- `@RG` e escrito pelo proprio STAR com `ID`, `SM`, `LB:TRUSEQ` e `PL:ILLUMINA`, separados por
  espaco no `--outSAMattrRGline` (o STAR os une com tab). Uma versao antiga escrevia `\\t`,
  que o bash virava `t` literal e gerava um `@RG` de campo unico, sem `SM`/`PL`: BAMs antigos
  assim quebram MarkDuplicates/BQSR. A rule `fix_read_group` e a rede de seguranca
  (`fix_read_group: true`): se o `@RG` ja tem `SM:`, faz hardlink; senao reescreve com `samtools addreplacerg`.
  Reparo retroativo de BAMs legados: `misc/repair-bam-readgroups.sh`.
- `--outSAMmapqUnique 60`: sem isso o STAR da MAPQ **255** aos alinhamentos unicos e o GATK
  le 255 como "qualidade indisponivel", descartando exatamente as reads que importam.
  BAMs de lotes anteriores a 2026-11-01 (L1-L3) tem 255; o `SplitNCigarReads` do GATK converte.
- Tambem saem `NM` e `MD` (`--outSAMattributes NH HI AS nM NM MD XS`).

## Saida: `work_dir` vs `final_dir`

```
work_dir  (scratch)  /mnt/scratch-shared/{user}/{project}/
    clean/ mapped/ counts/ logs/ rdata/ tmp/       <- intermediarios
final_dir (duravel)  results/   (relativo ao diretorio de trabalho, ou caminho absoluto)
    counts/ logs/ mapped/ rdata/  .export_done     <- entregaveis
```

`rule all` pede `final_dir/.export_done`: enquanto a sentinela nao existe, o workflow nao
terminou, por mais completo que o `work_dir` pareca. A sentinela tem **zero bytes de proposito**.

O que o export (`scripts/export_results.sh`) copia:

- `counts/`: `*_ReadsPerGene.out.tab`, `*.genes.results`, `*.isoforms.results`, `*_SJ.out.tab`, e `{amostra}.stat/` (modelo do RSEM)
- `logs/`: `multiqc_report.html` (+ `_data`), `*_Log.final.out`, JSON do fastp, `.cnt` do RSEM
- `mapped/`: depende de `keep_bams` (abaixo)
- `rdata/DEG.RData`, se a etapa de DEG rodou

e verifica por **conteudo** (nao existencia): nenhum arquivo de zero bytes, mesma contagem de
`genes.results`/`SJ.out.tab` na origem e no destino, e `samtools quickcheck` em 5 BAMs
sorteados (precisa de `samtools` no PATH, isto e, dentro do ambiente conda).

| `output.keep_bams` | Exporta |
|---|---|
| `none` | nenhum BAM; so contagens e relatorios |
| `coord` (padrao) | `*_Aligned.sortedByCoord.out.bam`, na verdade o `_rg.bam` (com `@RG` utilizavel) renomeado. **Permite variant calling** |
| `all` | acima + `*_Aligned.toTranscriptome.out.bam` (so serve ao RSEM) |

O BAM exportado **nao vem indexado**. Indexe a mao para o IGV (`misc/index-bams.sbatch`
faz isso em lote por SLURM).

### `cleanup_scratch`

Com `cleanup_scratch: true` (padrao), o `onsuccess` roda `scripts/cleanup_scratch.sh`, que
**apaga o `work_dir` inteiro** so depois de reconferir os entregaveis no `final_dir`. Na
duvida, preserva e sai 0. So aceita `work_dir` sob `/mnt/scratch-shared/`. Em falha, o
`work_dir` e **preservado** (e onde esta o material para diagnosticar e retomar).
A limpeza roda no no do orquestrador e **nao** alcanca o `local_scratch` de cada no.

### `local_scratch`: temporarios de ferramenta no disco do no

`output.local_scratch` (vazio por padrao) manda para o disco local do no apenas os
**temporarios de ferramenta** (`tmp/`: `_STARtmp` do STAR e `.rsemtmp` do RSEM, ~30 GB por
amostra so do RSEM). **Saidas declaradas de regra continuam no `work_dir` compartilhado**,
porque o Snakemake confere saidas a partir do no orquestrador: um arquivo no disco de um
compute node nao existe para o head (`MissingOutputException`, que ja apagou 48 GB de trabalho valido).

Exemplo em `config.local.yaml`: `output: {local_scratch: "/scratch/{user}/{project}"}`.
No r660xs (head), `/scratch/ronald` e **tmpfs** (RAM, 300 GB): la o local_scratch consome memoria.

> **Agrupamento por amostra:** o comentario em `common.smk` (bloco `if LOCAL:` e a mensagem em
> stderr) ainda diz que o `local_scratch` exige agrupamento por amostra via profile. **Isso deixou
> de ser verdade** quando
> `clean` voltou para o `work_dir` (commit `968555a`, que tambem removeu o `groups:` do profile):
> hoje so `tmp/` segue o `local_scratch`, e um temporario nao cria dependencia entre etapas.
> O codigo (`W` em `common.smk`) e o que vale. Se alguem religar `clean` no disco local, o
> agrupamento volta a ser obrigatorio.

## Elo com o variant calling

Os BAMs de `final_dir/mapped/` (`keep_bams: coord`) sao a **ENTRADA** do modo `rnaseq_mode`
do `interomics_snakemake` (`/home/ronald/Documentos/ronald/interomics_snakemake`):

- `rnaseq.bam_dirs` aponta para `.../hs_rnaseq_Lx/results/mapped`, **somente-leitura**.
  A amostra e o nome do BAM menos `rnaseq.bam_suffix` (`_Aligned.sortedByCoord.out.bam`), e tem de ser o `SM` do `@RG`.
- A rule `rna_conform_bam` renomeia os contigs `NC_` para `chr*`, reordena e descarta os 514
  contigs sem par, via `rna_contig_map.tsv`. Ver `documentation/REFERENCIA_UNIFICADA.md`.
- O BAM precisa ter `@RG` com `SM` e `PL`, e MAPQ 60 (ou 255, que o `SplitNCigarReads` converte).
- Config de producao que amarra os dois:
  `/home/ronald/Documentos/ronald/tese_natalia/corridas/vc_L1_L3/config_rnaseq_L1.yaml`
  (derivados para o L4 em `.../corridas/vc_L4/`).
- chrM esta fora do alvo de chamada de variantes (HaplotypeCaller e diploide).

## Estrutura

```
.
├── Snakefile                  # analise (nao constroi referencia)
├── reference.smk              # construcao da referencia, rodado 1x
├── common.smk                 # caminhos de referencia, work_dir/final_dir, guardas do config
├── config.yaml
├── scripts/
│   ├── build_reference.sh     # wrapper com flock
│   ├── export_results.sh      # work_dir -> final_dir, com verificacao
│   └── cleanup_scratch.sh     # apaga o work_dir so apos conferir o export
├── profiles/local/            # execucao local (padrao)
├── profiles/slurm/            # execucao via SLURM
├── metadata/                  # exp_design.csv (dado do projeto)
├── DEG.R  ETL_gene_info.R
├── misc/                      # SRA, @RG, BAMs e insumos do variant calling (ver abaixo)
├── documentation/             # runbooks e referencia unificada
├── samples/                   # FASTQ de entrada
└── pipeline.sh  pipeline_single_ended.sh   # legado, anteriores ao Snakefile
```

## Etapas

1. **trim_reads** - fastp remove adaptadores e bases de baixa qualidade
2. **map_reads** - STAR alinha ao genoma e ao transcriptoma, com read group, MAPQ unico 60
3. **fix_read_group** - garante `@RG` completo no BAM coordenado (hardlink se ja valido)
4. **quantify_transcripts** - RSEM quantifica por gene e transcrito
5. **run_multiqc** - consolida QC de fastp, STAR e RSEM
6. **run_deg_analysis** - DESeq2 + ReactomePA (desligue com `deg_analysis.enabled: false`)
7. **export_results** - copia entregaveis para o `final_dir` e escreve `.export_done`

> **Atencao ao `run_deg_analysis`.** O `deg_analysis.enabled` e `true` no `config.yaml`, mas o
> `DEG.R` ainda le `./output/counts` e grava `./output/rdata/DEG.RData` (caminhos relativos ao
> diretorio atual), enquanto o workflow agora escreve contagens no `work_dir`. Com o `work_dir`
> padrao a etapa **nao encontra as contagens**. Todos os lotes de producao rodaram com
> `deg_analysis.enabled: false` (o DEG so faz sentido com os lotes juntos). Trate a etapa
> como pendente de ajuste em `DEG.R` antes de liga-la.

## Design experimental

`metadata/exp_design.csv`, no formato de `metadata/exp_design.example.csv`:

```csv
sample,condition
SRR00000001,control
SRR00000003,case
```

Os valores de `condition` precisam ser `control` e `case` - o contraste esta fixo no
`DEG.R`. Cada grupo precisa de pelo menos 2 amostras. O arquivo e dado do projeto e esta no `.gitignore`.

## Scripts auxiliares (`misc/`)

```bash
bash misc/download_sra.sh lista.csv [dir_saida] [dir_temp]
bash misc/download_sra_single_ended.sh lista.csv [dir_saida] [dir_temp]
bash misc/merge_lanes.sh amostras_master.csv [dir_samples]   # funde lanes para {amostra}_{1,2}.fastq.gz
bash misc/fix_rg_flag.sh [dir_samples] [dir_mapped]          # legado: BAMs sem @RG algum
bash misc/repair-bam-readgroups.sh <dir_bams> [check|repair] [paralelismo]   # legado: @RG colapsado
bash misc/recompress-bams.sh <dir_bams> [check|run]          # BAM escrito sem compressao
sbatch misc/index-bams.sbatch                                # .bai em lote (job array)
```

O `dir_temp` do `download_sra.sh` default aponta para `/mnt/scratch-shared/$USER/sra_tmp`,
nao para `/tmp`: neste cluster `/tmp` e tmpfs em RAM e o `fasterq-dump` estoura a memoria la.

Insumos do variant calling (ver `documentation/REFERENCIA_UNIFICADA.md` §2):
`make-ncbi-dict.sbatch`, `make-contig-map.sbatch`, `compare-masked-contigs.sbatch`,
`make-exon-intervals.sbatch`. `monitor-array.sbatch` monitora o array de armazenamento durante medicoes de vazao.

## Problemas comuns

- **"Referencia ausente ou incompleta"** - falta a sentinela. Rode `scripts/build_reference.sh`,
  ou confira se `reference.dir` / `genome_id` / `star.read_length` batem com o que existe em disco.
- **"twopass e shared_genome nao podem estar ambos ligados"** - guarda do parse; desligue um.
- **"Outro build em andamento"** - o `flock` esta pego. Confirme que nao ha outro build rodando
  antes de remover `.build.lock`.
- **Permission denied ao escrever na referencia** - e proposital. Veja a nota de rebuild acima.
- **`MissingOutputException` num job que "terminou bem"** - saida declarada em disco local de
  compute node. Veja `local_scratch` acima.
- **`[cleanup] ABORTADO`** no log - alguma conferencia nao fechou; o `work_dir` foi preservado de proposito.
- Logs de cada etapa ficam em `{work_dir}/logs/` (e, apos o export, em `{final_dir}/logs/`); os do
  build, em `<reference>/<genome_id>/logs/`.
