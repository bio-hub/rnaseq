# RUNBOOK — rodar o workflow RNA-seq no r660xs

Instruções operacionais para executar este workflow neste servidor.
Leia tudo antes de executar qualquer comando.

---

## 1. Contexto

- **Host:** `r660xs` (head node + bigmem do cluster iLIKA). 128 threads, ~1 TB RAM.
- **Repo:** o diretório que contém este arquivo. Descubra onde você está:
  ```bash
  cd "$(dirname "$(readlink -f RUNBOOK.md)")" 2>/dev/null || cd ~/rnaseq
  pwd
  ```
  Cada usuário tem a própria cópia (ex.: `/mnt/home/cliente01/rnaseq`).
  **Não assuma caminhos de outro usuário** — os homes são `700` e você não os lê.
- **Branch:** `refactor/reference-path` — **a `main` NÃO tem este workflow.**
  O Snakefile da `main` tem erro de sintaxe e nem parseia. Confira a branch antes
  de qualquer coisa.
- **Ambiente conda:** `biohub_rnaseq`, na instalação **compartilhada**
  `/mnt/shared/conda` (ver §5.1).
- **Referência:** `/mnt/reference/rnaseq/GRCh38_p14/` — **já construída**, 40 GB,
  somente-leitura.

O workflow faz: fastp → STAR → RSEM → MultiQC → DESeq2/ReactomePA.

---

## 2. O que já está pronto — NÃO refaça

| Item | Estado |
|---|---|
| Ambiente conda `biohub_rnaseq` | criado, 14 binários e 10 pacotes R validados |
| Referência GRCh38_p14 | construída, sentinela `.done_sjdb149` presente |
| Índice STAR/RSEM | `rsem_star_sjdb149/`, 33 GB |
| `gene_info.tsv` | 193.953 genes |
| sra-tools | configurado (`~/.ncbi/user-settings.mkfg`) |

---

## 3. Proibições

**Não rode `scripts/build_reference.sh`.** A referência existe. O script leva horas
e ocupa 33 GB. Só se justifica para criar um índice de **outro** `read_length`
(ver §7).

**Não faça `chmod` na referência.** Os diretórios em
`/mnt/reference/rnaseq/GRCh38_p14/` estão `dr-xr-xr-x` de propósito, para impedir
que uma execução acidental apague o índice. Se algo falhar com *Permission denied*
ao escrever lá, isso é o mecanismo funcionando — investigue a causa, não contorne.

**Não use SLURM.** Existe `profiles/slurm/`, mas só o r660xs monta
`/mnt/reference` (ainda não há NFS exportando o volume). Jobs em outros nós não
enxergam a referência. Use `profiles/local`.

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

cd ~/rnaseq                       # ou o diretório onde está este RUNBOOK
git switch refactor/reference-path
```

Confirme que o env veio do lugar certo — se `which snakemake` apontar para o
home de outro usuário, você não terá permissão de leitura:

```bash
which snakemake                   # esperado: /mnt/shared/conda/envs/biohub_rnaseq/bin/snakemake
snakemake --version               # esperado: 9.x
```

### 5.2 Conferir o comprimento das reads

`config.yaml` traz `star.read_length: 150`, que deriva `sjdbOverhang = 149`.
**Só existe o índice `sjdb149`.** Se as reads deste projeto não forem de 150 bp,
pare e leia §7 — não altere o config sem construir o índice correspondente.

### 5.3 Dry-run (sempre antes de executar)

```bash
snakemake --workflow-profile profiles/local -n
```

Confirme na saída:

- o número de jobs bate com o número de amostras esperado;
- não há `MissingInputException`;
- não há `WorkflowError` sobre referência ausente.

### 5.4 Rodar

É longo (horas, conforme o número de amostras). Use `tmux`:

```bash
REPO=$(pwd)                       # estando no diretório do repo
mkdir -p output/logs
tmux new-session -d -s rnaseq-run
tmux send-keys -t rnaseq-run \
  "source /mnt/shared/conda/etc/profile.d/conda.sh && \
   conda activate biohub_rnaseq && \
   cd $REPO && \
   snakemake --workflow-profile profiles/local 2>&1 | tee output/logs/run_\$(date +%Y%m%d_%H%M%S).log" C-m
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

## 6. Verificação de sucesso

Não confie no exit code sozinho. Confira os artefatos:

```bash
ls output/counts/*.genes.results        # 1 por amostra (RSEM)
ls output/counts/*_ReadsPerGene.out.tab # 1 por amostra (STAR)
ls output/logs/multiqc_report.html      # QC consolidado
ls output/rdata/DEG.RData               # só se deg_analysis.enabled
```

Taxa de mapeamento, por amostra:

```bash
grep -H "Uniquely mapped reads %" output/logs/*_Log.final.out
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
# nenhum output pode ter zero bytes
find output -type f -size 0 | head

# BAMs devem ser legíveis de verdade, não só existir
for b in output/mapped/*_Aligned.sortedByCoord.out.bam; do
    samtools quickcheck "$b" || echo "CORROMPIDO: $b"
done

# contagens devem ter linhas além do cabeçalho
wc -l output/counts/*.genes.results | tail -1
```

Se o host reiniciar no meio de uma execução, **não presuma que dá para continuar**:
rode `snakemake --workflow-profile profiles/local -n` e confira o que ele planeja
refazer; o profile já traz `rerun-incomplete: true`.

---

## 7. Reads que não são de 150 bp

O `sjdbOverhang` ideal é `read_length - 1`, e ele entra no nome do índice e da
sentinela, então índices de comprimentos diferentes coexistem.

**Só existe `sjdb149` hoje.** Para reads de 100 bp é preciso construir o índice
de 99 — leva ~25 min e ~33 GB, e reaproveita FASTA e GTF já baixados:

```bash
printf 'star:\n  read_length: 100\n' > config.local.yaml
bash scripts/build_reference.sh
snakemake --configfile config.yaml config.local.yaml \
          --workflow-profile profiles/local
```

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

Logs por etapa ficam em `output/logs/`.

---

## 9. Quando escalar para o usuário

- Qualquer rule falhar duas vezes seguidas pelo mesmo motivo.
- Taxa de mapeamento anômala (§6).
- Necessidade de construir índice novo (§7).
- Qualquer impulso de usar `chmod`/`rm` em `/mnt/reference`.
- Qualquer commit ou push.
