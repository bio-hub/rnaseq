# RNA-seq Analysis Workflow

Workflow Snakemake para RNA-seq paired-end: trimagem, alinhamento, quantificacao,
QC e analise de expressao diferencial.

A **referencia e um insumo externo**, compartilhado e somente-leitura. Ela vive fora
do diretorio de trabalho, em `config["reference"]["dir"]` (default `/mnt/reference/rnaseq`),
e e construida **uma unica vez** por um workflow separado. Nenhuma execucao de analise
baixa ou reconstroi a referencia.

## Requisitos

- Ambiente conda de `environment.yaml`
- FASTQ paired-end em `samples/`, nomeados `{sample}_1.fastq.gz` e `{sample}_2.fastq.gz`
- A referencia ja construida (ver abaixo)

## Instalacao

```bash
conda env create -f environment.yaml
conda activate biohub_rnaseq
```

## Passo 1 - construir a referencia (uma vez por genoma)

```bash
# o diretorio raiz precisa existir e ser gravavel
sudo mkdir -p /mnt/reference/rnaseq
sudo chown "$USER":"$USER" /mnt/reference/rnaseq

bash scripts/build_reference.sh
```

Isso baixa genoma, GTF e o `gene_info` do NCBI e constroi o indice STAR/RSEM em:

```
/mnt/reference/rnaseq/GRCh38_p14/
├── fasta/        gtf/        gene_info/
├── rsem_star_sjdb100/        # o sjdbOverhang entra no nome
├── logs/         MANIFEST.txt
└── .done_sjdb100             # sentinela: so existe se o build terminou
```

Reserve ~40 GB de disco e algumas horas. O script segura um `flock` exclusivo: se
voce disparar dois builds ao mesmo tempo, o segundo aborta em vez de corromper o
indice. Ao terminar, os arquivos ficam somente-leitura. Para reconstruir de proposito:

```bash
chmod -R u+w /mnt/reference/rnaseq/GRCh38_p14
```

## Passo 2 - rodar a analise

```bash
conda activate biohub_rnaseq
snakemake --workflow-profile profiles/local
```

O profile local ja define `cores` e `resources.mem_mb`. O `mem_mb` e o que impede o
Snakemake de agendar jobs pesados demais em paralelo - sem ele, varios `map_reads`
de ~70 GB sao disparados juntos. Ajuste `profiles/local/config.yaml` para a maquina.

Para planejar sem executar:

```bash
snakemake --workflow-profile profiles/local -n
```

Opcionalmente ha um perfil SLURM em `profiles/slurm/` (particao `bigmem`, que e o
unico no que enxerga `/mnt/reference` enquanto nao houver NFS). Nao e necessario.

## Configuracao

Tudo em `config.yaml`. Os campos que mais importam:

| Campo | Para que serve |
|---|---|
| `reference.dir` | raiz das referencias compartilhadas |
| `reference.genome_id` | subdiretorio do build; um por genoma |
| `reference.tax_id` | especie filtrada no `gene_info` (9606 = humano) |
| `star.sjdbOverhang` | ideal = comprimento da read - 1. **Precisa ser o mesmo no build e no mapeamento** |
| `threads.*` / `mem_mb.*` | recursos por rule |
| `deg_analysis.*` | design experimental, alpha e limiar de log2FC |

Para rodar em outra maquina ou caminho, nao edite o `config.yaml` versionado - crie um
`config.local.yaml` (ja no `.gitignore`) com apenas o que muda:

```yaml
reference:
  dir: "/caminho/alternativo"
```

```bash
snakemake --configfile config.yaml config.local.yaml --workflow-profile profiles/local
```

O `sjdbOverhang` entra no nome da sentinela (`.done_sjdb100`). Se voce mudar o valor no
config sem reconstruir o indice, o workflow falha no parse com mensagem explicita, em
vez de deixar o STAR abortar no meio do mapeamento.

## Estrutura

```
.
├── Snakefile                  # analise (nao constroi referencia)
├── reference.smk              # construcao da referencia, rodado 1x
├── common.smk                 # caminhos derivados de config.reference
├── config.yaml
├── scripts/build_reference.sh # wrapper com flock
├── profiles/local/            # execucao local (padrao)
├── profiles/slurm/            # opcional
├── metadata/                  # exp_design.csv (dado do projeto)
├── DEG.R  ETL_gene_info.R
├── misc/                      # download do SRA, correcao de read group
├── samples/                   # FASTQ de entrada
└── output/
    ├── clean/ mapped/ counts/ logs/ rdata/
```

## Etapas

1. **trim_reads** - fastp remove adaptadores e bases de baixa qualidade
2. **map_reads** - STAR alinha ao genoma e ao transcriptoma, com read group
3. **quantify_transcripts** - RSEM quantifica por gene e transcrito
4. **run_multiqc** - consolida QC de fastp, STAR e RSEM
5. **run_deg_analysis** - DESeq2 + ReactomePA (desligue com `deg_analysis.enabled: false`)

## Design experimental

`metadata/exp_design.csv`, no formato de `metadata/exp_design.example.csv`:

```csv
sample,condition
SRR00000001,control
SRR00000003,case
```

Os valores de `condition` precisam ser `control` e `case` - o contraste esta fixo no
`DEG.R`. Cada grupo precisa de pelo menos 2 amostras.

## Scripts auxiliares

```bash
bash misc/download_sra.sh lista.csv [dir_saida] [dir_temp]
bash misc/download_sra_single_ended.sh lista.csv [dir_saida] [dir_temp]
bash misc/fix_rg_flag.sh [dir_samples] [dir_mapped]
```

O `dir_temp` default aponta para `/mnt/scratch-shared`, nao para `/tmp`: neste cluster
`/tmp` e tmpfs em RAM e o `fasterq-dump` estoura a memoria la.

## Problemas comuns

- **"Referencia ausente ou incompleta"** - falta a sentinela. Rode `scripts/build_reference.sh`,
  ou confira se `reference.dir` / `genome_id` / `star.sjdbOverhang` batem com o que existe em disco.
- **"Outro build em andamento"** - o `flock` esta pego. Confirme que nao ha outro build rodando
  antes de remover `.build.lock`.
- **Permission denied ao escrever na referencia** - e proposital. Veja a nota de rebuild acima.
- Logs de cada etapa ficam em `output/logs/`; os do build, em `<reference>/<genome_id>/logs/`.
