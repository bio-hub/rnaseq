# RUNBOOK — rodar o workflow RNA-seq no r660xs

Instruções operacionais para executar este workflow neste servidor.
Leia tudo antes de executar qualquer comando.

---

## 1. Contexto

- **Host:** `r660xs` (head node + bigmem do cluster iLIKA). 128 threads, ~1 TB RAM.
- **Repo:** este arquivo mora em `documentation/` dentro do repo; a raiz do repo é
  um nível acima (`cd "$(dirname "$(readlink -f documentation/RUNBOOK.md)")/.."`).
  Cada lote roda num clone próprio em **`/mnt/shared/projetos/<projeto>/`**, não
  em `~`: o home do administrador é local em cada nó e some de vista dos workers.
  Tenants têm home em `/mnt/home/<user>`, que é `700` — **não assuma caminhos de
  outro usuário**, você não os lê.
- **Branch:** `main`.

  > Até 31/10 este runbook dizia que a `main` não tinha o workflow e que o
  > Snakefile dela nem parseava. **Isso deixou de valer**: o refactor foi
  > mergeado na `main` (PR #1, 12 commits) e o Snakefile de lá carrega —
  > conferido com `snakemake --list` em 01/11. As branches
  > `refactor/reference-path` e `feat/scratch-workdir-and-rg-fix` continuam no
  > repositório, mas são históricas. Use a `main`.
- **Ambiente conda:** `biohub_rnaseq`, na instalação **compartilhada**
  `/mnt/shared/conda` (ver §5.1).
- **Referência:** `/mnt/reference/rnaseq/GRCh38_p14/` — **já construída**, 40 GB,
  somente-leitura. É o **NCBI RefSeq GRCh38.p14** (contigs `NC_000001.11`, 705
  contigs), **não** o genoma do `interomics_snakemake` (`chr1`, 195 contigs).
  Leia `REFERENCIA_UNIFICADA.md` antes de levar BAMs deste workflow ao variant calling.

O workflow faz: fastp → STAR → RSEM → MultiQC → (DESeq2/ReactomePA) → export dos
entregáveis. Intermediários vão para o `work_dir` (scratch) e só os entregáveis
para o `final_dir` (§5.5).

Outros documentos: `RUNBOOK_SLURM.md` (execução via SLURM, medições, distribuição
entre nós), `REFERENCIA_UNIFICADA.md` (divergência de genoma com o variant calling),
`../README.md` (visão geral e configuração).

---

## 2. O que já está pronto — NÃO refaça

| Item | Estado |
|---|---|
| Ambiente conda `biohub_rnaseq` | criado, 14 binários e 10 pacotes R validados |
| Referência GRCh38_p14 | construída; sentinelas `.done_sjdb{74,75,99,124,149}` |
| Índices STAR/RSEM | `rsem_star_sjdb{74,75,99,124,149}/`, ~33 GB cada (um por `read_length`: 75, 76, 100, 125, 150) |
| `gene_info.tsv` | 193.952 genes (193.953 linhas com cabeçalho) |
| sra-tools | configurado (`~/.ncbi/user-settings.mkfg`) |

---

## 3. Proibições

**Não rode `scripts/build_reference.sh`** sem necessidade. A referência existe, e o
script leva horas e ocupa ~33 GB por índice. Só se justifica para criar um índice
de um `read_length` que **não** esteja entre os cinco existentes (ver §7).

**Não faça `chmod` na referência.** Os diretórios em
`/mnt/reference/rnaseq/GRCh38_p14/` estão `dr-xr-xr-x` de propósito, para impedir
que uma execução acidental apague o índice. Se algo falhar com *Permission denied*
ao escrever lá, isso é o mecanismo funcionando — investigue a causa, não contorne.

**SLURM não é proibido.** Até 31/10 este runbook dizia que só o r660xs montava
`/mnt/reference`. Isso deixou de valer: `/mnt/reference`, `/mnt/shared`,
`/mnt/scratch-shared` e `/mnt/home` são NFSv4 nos quatro nós. Este runbook cobre
a execução **local** (`profiles/local`); para fila, use `RUNBOOK_SLURM.md`.

**Não use `/tmp` para arquivos grandes.** Aqui `/tmp` é **tmpfs em RAM** (502 GB).
Escrever dezenas de GB lá consome memória real.

**Não faça commit nem push** sem pedir ao usuário.

---

## 4. Pré-requisitos de dados (fornecidos pelo usuário)

### 4.1 FASTQ

Em `samples/`, paired-end, exatamente neste padrão:

```
samples/{amostra}_1.fastq.gz
samples/{amostra}_2.fastq.gz
```

O nome da amostra é inferido removendo `_1.fastq.gz`. Arquivos fora desse padrão
são ignorados silenciosamente — se `SAMPLES` sair vazio, o workflow não faz nada
e não avisa. **Confira a contagem antes de rodar** (§5.3).

Single-end não é suportado por este Snakefile.

### 4.2 Design experimental

`metadata/exp_design.csv`:

```csv
sample,condition
SRR00000001,control
SRR00000002,control
SRR00000003,case
SRR00000004,case
```

Requisitos rígidos:

- Coluna `sample` deve bater **exatamente** com os nomes inferidos dos FASTQ.
- Coluna `condition` aceita **apenas** os literais `control` e `case` — o contraste
  está fixo em `DEG.R` (`contrast = c("condition","case","control")`).
- **Mínimo 2 amostras por grupo.** Com menos, o filtro de contagens baixas e o
  DESeq2 falham.

Se ainda não houver design, desligue a etapa de DEG em vez de improvisar um CSV:

```yaml
# config.yaml
deg_analysis:
  enabled: false
```

---

## 5. Execução

### 5.1 Ativar o ambiente

Em shell não-interativo, `conda activate` só funciona após carregar o profile.
Use a instalação **compartilhada**, que é a única acessível a todos os usuários:

```bash
source /mnt/shared/conda/etc/profile.d/conda.sh
conda activate biohub_rnaseq

cd /mnt/shared/projetos/<projeto> # clone do repo, na main (ver §1)
git branch --show-current         # esperado: main
```

Confirme que o env veio do lugar certo — se `which snakemake` apontar para o
home de outro usuário, você não terá permissão de leitura:

```bash
which snakemake                   # esperado: /mnt/shared/conda/envs/biohub_rnaseq/bin/snakemake
snakemake --version               # esperado: 9.x
```

### 5.2 Conferir o comprimento das reads

`config.yaml` traz `star.read_length: 150`, que deriva `sjdbOverhang = 149`.
Existem índices para `read_length` 75, 76, 100, 125 e 150 (`sjdb74`, `75`, `99`,
`124`, `149`). Se as reads do projeto têm outro comprimento, pare e leia §7. Para
os comprimentos existentes, ponha `star.read_length` no `config.local.yaml` (nunca
no `config.yaml` versionado). O workflow **não detecta** o comprimento dos FASTQ.

Confira também as duas otimizações do STAR, ambas `false` por padrão e
**mutuamente exclusivas** (o parse falha se as duas estiverem ligadas):
`star.twopass` (recomendado para variant calling; dobra o `map_reads`) e
`star.shared_genome` (ganho de ~3,5%, deixa segmento de shm órfão se um job morre;
não use em fila). Para DEG puro, deixe as duas `false`. Detalhes no `../README.md`.

### 5.3 Dry-run (sempre antes de executar)

```bash
snakemake --configfile config.yaml config.local.yaml \
          --workflow-profile profiles/local -n
```

(sem `config.local.yaml`, omita `--configfile`.) Confirme na saída:

- o número de jobs bate com o número de amostras esperado;
- não há `MissingInputException`;
- não há `WorkflowError` sobre referência ausente.

### 5.4 Rodar

É longo (horas, conforme o número de amostras). Use `tmux`:

```bash
REPO=$(pwd)                       # estando no diretório do repo
tmux new-session -d -s rnaseq-run
tmux send-keys -t rnaseq-run \
  "source /mnt/shared/conda/etc/profile.d/conda.sh && \
   conda activate biohub_rnaseq && \
   cd $REPO && \
   snakemake --configfile config.yaml config.local.yaml \
     --workflow-profile profiles/local 2>&1 | tee run_\$(date +%Y%m%d_%H%M%S).log" C-m
```

> **Use `tmux` sempre.** Sessões de SSH caem e o workflow leva horas. O `tmux`
> sobrevive à queda da conexão — mas **não** sobrevive a reboot do host.

Acompanhar:

```bash
tmux attach -t rnaseq-run     # Ctrl+B D para sair sem matar
```

O profile já define `cores: 120` e `resources.mem_mb: 900000`. **Não aumente o
`mem_mb`**: ele é o que limita a concorrência — sem teto, vários `map_reads` de
~70 GB são agendados juntos.

---

### 5.5 Onde ficam os arquivos

| | Caminho | Conteúdo |
|---|---|---|
| `work_dir` | `/mnt/scratch-shared/<user>/<project>/` (`output.project`, default `rnaseq`) | `clean/ mapped/ counts/ logs/ rdata/ tmp/` — intermediários |
| `final_dir` | `results/` no diretório do repo (ou caminho absoluto) | `counts/ logs/ mapped/ rdata/` — entregáveis, mais a sentinela `.export_done` |

- **Dê um `output.project` por lote** em `config.local.yaml`; senão dois lotes
  disputam o mesmo `work_dir`.
- `output.keep_bams`: `coord` (padrão) exporta o BAM coordenado, que é a entrada do
  variant calling; `none` só contagens; `all` inclui o `toTranscriptome`.
- `output.cleanup_scratch: true` (padrão) **apaga o `work_dir`** ao fim com sucesso,
  depois de reconferir o `final_dir`. Em falha o `work_dir` é preservado.
- `output.local_scratch` (vazio por padrão) só afeta temporários de ferramenta
  (STAR, RSEM). Não preencha em execução local: no r660xs `/scratch/ronald` é tmpfs
  de 300 GB, ou seja, RAM.
- O BAM exportado **não vem indexado** (`misc/index-bams.sbatch` indexa em lote).
- `.export_done` tem zero bytes **de propósito**.

> **DEG:** `deg_analysis.enabled` é `true` no `config.yaml`, mas o `DEG.R` ainda lê
> `./output/counts` e grava `./output/rdata/`, caminhos que o workflow não usa mais.
> Os lotes de produção rodaram com `deg_analysis.enabled: false` (o DEG só faz
> sentido com os lotes juntos). Desligue-o em `config.local.yaml` até o `DEG.R`
> ser ajustado.

---

## 6. Verificação de sucesso

Não confie no exit code sozinho. Confira os artefatos:

```bash
ls results/counts/*.genes.results        # 1 por amostra (RSEM)
ls results/counts/*_ReadsPerGene.out.tab # 1 por amostra (STAR)
ls results/counts/*_SJ.out.tab           # 1 por amostra (junções)
ls results/logs/multiqc_report.html      # QC consolidado
ls results/mapped/*.bam                  # só se keep_bams != none
ls results/.export_done                  # sentinela do export
ls results/rdata/DEG.RData               # só se deg_analysis.enabled
```

(Com `cleanup_scratch: true` o `work_dir` já foi apagado; a verificação é no `results/`.
Se o `final_dir` for outro caminho, troque `results/`.)

Taxa de mapeamento, por amostra:

```bash
grep -H "Uniquely mapped reads %" results/logs/*_Log.final.out
```

Abaixo de ~70% para RNA-seq humano merece investigação (contaminação, adaptador
residual, espécie errada) — **reporte ao usuário, não siga adiante em silêncio.**

### 6.1 Arquivo existir não é arquivo estar certo

Houve uma queda de energia neste host em 2026-10-29. O resultado foi uma árvore
de **1867 arquivos com zero bytes** e estrutura de diretórios intacta: o journal
do XFS restaura metadados, não os dados. Uma verificação do tipo
`if [ -d saida/ ]` passa nesse estado.

Então, ao conferir resultados, cheque **conteúdo**, não existência:

```bash
# nenhum output pode ter zero bytes (.export_done tem, de propósito)
find results -type f -size 0 ! -name .export_done | head

# BAMs devem ser legíveis de verdade, não só existir
for b in results/mapped/*_Aligned.sortedByCoord.out.bam; do
    samtools quickcheck "$b" || echo "CORROMPIDO: $b"
done

# contagens devem ter linhas além do cabeçalho
wc -l results/counts/*.genes.results | tail -1

# @RG utilizável pelo variant calling: ID, SM, LB e PL como campos separados
samtools view -H results/mapped/<amostra>_Aligned.sortedByCoord.out.bam | grep '^@RG'
```

O próprio `export_results.sh` já faz parte disso (zero bytes, contagens
origem×destino, `quickcheck` em 5 BAMs sorteados); o `EXPORT FALHOU` no
`work_dir/logs/export.log` vale mais que o exit code do Snakemake.

Se o host reiniciar no meio de uma execução, **não presuma que dá para continuar**:
rode `snakemake --workflow-profile profiles/local -n` e confira o que ele planeja
refazer; o profile já traz `rerun-incomplete: true`. O `work_dir` é preservado em falha.

---

## 7. Reads que não são de 150 bp

O `sjdbOverhang` ideal é `read_length - 1`, e ele entra no nome do índice e da
sentinela, então índices de comprimentos diferentes coexistem.

**Existem cinco hoje:** `sjdb74`, `75`, `99`, `124` e `149` (reads de 75, 76, 100,
125 e 150 bp). Para outro comprimento — por exemplo 51 bp — é preciso construir o
índice (`sjdb50`): leva ~25 min e ~33 GB, e reaproveita FASTA e GTF já baixados:

```bash
printf 'star:\n  read_length: 51\n' > config.local.yaml
bash scripts/build_reference.sh
snakemake --configfile config.yaml config.local.yaml \
          --workflow-profile profiles/local
```

Reads de comprimento misto: `star.sjdbOverhang` (maior − 1) tem precedência sobre
`read_length`. Reaproveitar um índice vizinho (ex.: `sjdb99` para reads de 76 bp)
é possível por esse override, mas não é o que o STAR recomenda.

Nunca edite `star.read_length` no `config.yaml` versionado para isto — use
`config.local.yaml`, que está no `.gitignore`.

**Confirme com o usuário antes de construir um índice novo.**

---

## 8. Problemas conhecidos

**`WorkflowError: Referencia ausente ou incompleta`**
A sentinela esperada não existe. A mensagem imprime `dir`, `genome_id` e a origem
do overhang. Quase sempre é `read_length` divergente do índice disponível (§7).
Não é para "resolver" criando a sentinela à mão.

**`Outro build ja esta em andamento`**
O `flock` está preso. Confirme que não há build rodando (`pgrep -f reference.smk`)
antes de considerar remover `.build.lock`.

**`Permission denied` escrevendo na referência**
Comportamento intencional (§3). Não contorne.

**Nenhum job no dry-run**
`SAMPLES` saiu vazio. Confira o padrão dos nomes em `samples/` (§4.1).

**Erros de `fasterq-dump` nos scripts de `misc/`**
O temp default é `/mnt/scratch-shared/$USER/sra_tmp`. Passe `-t` apontando para
disco, nunca `/tmp`.

**`WorkflowError: twopass e shared_genome ...`**
As duas opções do STAR estão ligadas. São mutuamente exclusivas (§5.2).

**`MissingOutputException` com SLURM/`local_scratch`**
Saída declarada em disco local de compute node. Ver `RUNBOOK_SLURM.md` §9.1.

**`[cleanup] ABORTADO`**
Uma conferência do `cleanup_scratch.sh` não fechou. O `work_dir` foi preservado
de propósito; investigue antes de apagar à mão.

Logs por etapa ficam em `<work_dir>/logs/` durante a execução e em
`<final_dir>/logs/` depois do export.

---

## 9. Quando escalar para o usuário

- Qualquer rule falhar duas vezes seguidas pelo mesmo motivo.
- Taxa de mapeamento anômala (§6).
- Necessidade de construir índice novo (§7).
- Qualquer impulso de usar `chmod`/`rm` em `/mnt/reference`.
- Qualquer commit ou push.

---

## 10. Depois do RNA-seq: variant calling

Os BAMs de `results/mapped/` (`keep_bams: coord`) são a **entrada** do modo
`rnaseq_mode` do `interomics_snakemake`
(`/home/ronald/Documentos/ronald/interomics_snakemake`). Cuidados:

- O BAM é **NCBI (`NC_...`, 705 contigs)**; o variant calling é `chr*` (195). A rule
  `rna_conform_bam` faz a ponte via `rna_contig_map.tsv` e **descarta** 8,7–19,0% dos
  registros (alt loci, patches, secundários). Ver `REFERENCIA_UNIFICADA.md`.
- `rnaseq.bam_dirs` do config do interomics aponta para `results/mapped` do lote,
  somente-leitura. Exemplo de produção:
  `/home/ronald/Documentos/ronald/tese_natalia/corridas/vc_L1_L3/config_rnaseq_L1.yaml`
  (para o L4: `.../corridas/vc_L4/config_rnaseq_L4{a,b}.yaml`).
- Exija `@RG` com `SM`/`PL` (ver §6) e, para lotes novos, MAPQ 60 (o workflow já
  passa `--outSAMmapqUnique 60`; L1–L3 têm 255, que o `SplitNCigarReads` converte).
- `star.twopass` é recomendado para variant calling (GATK Best Practices), mas
  muda as contagens; decida **antes** de mapear o lote.
- Não combine, numa mesma coorte, amostras mapeadas em referências diferentes.
