# RUNBOOK — rodar um lote via SLURM

Complementa `RUNBOOK.md`, que cobre a execução local. Leia aquele primeiro: as
proibições, a verificação por conteúdo e as armadilhas de dado valem igual aqui.

---

## 1. Leia isto antes de decidir usar SLURM

**Hoje o SLURM neste cluster é uma fila, não um cluster.**

Não há NFS. Apenas o `r660xs` monta `/mnt/reference`, `/mnt/home`,
`/mnt/shared` e `/mnt/scratch-shared`. Um job que caia no `r550-1` ou no
`r760` não enxerga a referência, nem os FASTQ, nem o conda compartilhado — e
falha de formas confusas, porque o caminho existe no `fstab` do head node e
não no do worker.

Por isso o `profiles/slurm/config.yaml` fixa **`slurm_partition: bigmem`**,
que é exatamente o `r660xs`. Ou seja: os jobs rodam na mesma máquina em que
você já rodaria localmente.

**O que você ganha mesmo assim:**

- A fila não depende de uma sessão `tmux` viva
- Contabilidade por job (`sacct`): tempo, memória de pico, exit code
- Limite de concorrência por recurso, não por "quantos couberam"
- O `prolog`/`epilog` cria e remove scratch isolado por job
- O mesmo comando funcionará sem alteração quando houver NFS

**O que você NÃO ganha:** paralelismo além dos 128 threads do head node.

> Para usar `cpu` e `gpu`, é preciso exportar os volumes por NFS e replicar
> usuários e grupos com os **mesmos UID/GID** nos 4 nós. Ver plano de rede,
> §10.4 e §13.6.

---

## 2. Pré-requisitos (verificados em 2026-10-30)

| Item | Estado |
|---|---|
| SLURM controller | ativo no `r660xs` |
| Partição `bigmem` | `r660xs`, 128 CPUs, ~1003 GB, limite 5 dias |
| `snakemake-executor-plugin-slurm` | instalado no env `biohub_rnaseq` |
| Submissão como `ronald` | funciona, sem accounting restritivo |
| `prolog`/`epilog` | criam e removem `/mnt/scratch-shared/jobs/$SLURM_JOB_ID` |
| Volumes no nó alocado | os 4 visíveis |

Confira antes de começar:

```bash
sinfo -s
source /mnt/shared/conda/etc/profile.d/conda.sh && conda activate biohub_rnaseq
snakemake --version                       # 9.x
python -c "import snakemake_executor_plugin_slurm; print('plugin ok')"
```

---

## 3. Rodando o lote 2

O L2 são 8 amostras do GSE115099, reads de **125 bp**, índice `sjdb124`
(já construído).

### 3.1 Preparar o diretório de trabalho

Como `ronald`, use um workdir seu — não o do `cliente01`:

```bash
mkdir -p ~/projetos/hs_rnaseq_L2 && cd ~/projetos/hs_rnaseq_L2
git clone -b main https://github.com/bio-hub/rnaseq.git .
mkdir -p samples metadata
```

### 3.2 Baixar os FASTQ

```bash
source /mnt/shared/conda/etc/profile.d/conda.sh && conda activate biohub_rnaseq

# o CSV do lote está em tese_natalia/workflow_inputs/misc/
bash misc/download_sra.sh \
     ~/Documentos/ronald/tese_natalia/workflow_inputs/misc/runs_L2_sjdb124.csv \
     /mnt/shared/projetos/hs_rnaseq/samples \
     /mnt/scratch-shared/$USER/sra_tmp
```

FASTQ vão para `/mnt/shared`, **nunca para o scratch** — ver `CLAUDE.md` §3.
Depois, symlink para o workdir:

```bash
ln -s /mnt/shared/projetos/hs_rnaseq/samples/*.fastq.gz samples/
```

### 3.3 Config do lote

**Nunca edite o `config.yaml` versionado.**

```bash
cat > config.local.yaml <<'EOF'
star:
  read_length: 125          # deriva sjdbOverhang=124
deg_analysis:
  enabled: false            # DEG só quando os três lotes terminarem
output:
  project: "hs_rnaseq_L2"   # separa o scratch por lote
  keep_bams: "coord"
EOF
```

### 3.4 Dry-run — obrigatório

```bash
snakemake --configfile config.yaml config.local.yaml \
          --workflow-profile profiles/slurm -n
```

Confirme: o número de jobs bate com 8 amostras, e não há `WorkflowError`
sobre referência ausente. Se houver, o `read_length` não corresponde a nenhum
índice construído.

### 3.5 Submeter

O processo do snakemake precisa continuar vivo para orquestrar a fila — ele
submete, acompanha e encadeia. Portanto **ainda precisa de `tmux`**, mas agora
o trabalho pesado está na fila, não no processo.

```bash
tmux new-session -d -s l2
tmux send-keys -t l2 'source /mnt/shared/conda/etc/profile.d/conda.sh && \
  conda activate biohub_rnaseq && cd ~/projetos/hs_rnaseq_L2 && \
  snakemake --configfile config.yaml config.local.yaml \
  --workflow-profile profiles/slurm \
  2>&1 | tee logs_L2_$(date +%Y%m%d_%H%M%S).log' C-m

sleep 15
tmux capture-pane -p -t l2 | tail -15      # NÃO pule: confirme que arrancou
```

> `tmux new-session` **falha se a sessão já existir**, e o `send-keys`
> seguinte vai para a sessão velha. O sintoma é "nada aconteceu".

Alternativa sem `tmux`, para deixar o orquestrador ele mesmo na fila:

```bash
sbatch --partition=bigmem --cpus-per-task=2 --mem=4G --time=3-00:00:00 \
       --job-name=snakemake_L2 --output=orquestrador_%j.out \
       --wrap="source /mnt/shared/conda/etc/profile.d/conda.sh && \
               conda activate biohub_rnaseq && cd ~/projetos/hs_rnaseq_L2 && \
               snakemake --configfile config.yaml config.local.yaml \
               --workflow-profile profiles/slurm"
```

Nesse caso o orquestrador sobrevive a logout **e a reboot do seu terminal** —
mas não ao reboot do nó.

---

## 4. Acompanhar

```bash
squeue -u $USER                      # fila agora
squeue -u $USER -o "%.8i %.20j %.8T %.10M %.6D %R"
sacct -u $USER --starttime today \
      --format=JobID,JobName%22,State,Elapsed,MaxRSS,ExitCode
```

`MaxRSS` é o número que interessa para ajustar `mem_mb` na próxima rodada: se
ficar muito abaixo do reservado, você está serializando o DAG à toa.

Job travado em `PD` (pending) — veja o motivo:

```bash
squeue -u $USER -o "%.8i %.9P %.8T %R"
```

`Resources` significa fila por recurso; `Priority`, que há algo na frente.

---

## 5. Recursos: onde eles são definidos

Três camadas, da mais geral para a mais específica:

| Camada | Arquivo | O que define |
|---|---|---|
| Default | `profiles/slurm/config.yaml` | partição, runtime e memória padrão |
| Por rule | `config.yaml` → `threads`, `mem_mb`, `runtime_min` | o que cada etapa pede |
| Override | `config.local.yaml` | ajuste do lote |

O plugin converte `threads` em `--cpus-per-task`, `resources.mem_mb` em
`--mem` e `resources.runtime` (minutos) em `--time`.

> **`runtime` é um teto rígido.** O SLURM mata o job ao atingi-lo. O
> `map_reads` está em 240 min por amostra; se um lote tiver amostras muito
> profundas, suba antes de submeter — não depois de perder 4 horas.

---

## 6. Quando algo falha

**O job morre e o snakemake continua** — por padrão ele tenta outras rules e
só então aborta. Procure `Error in rule` no log e o `JobID` correspondente:

```bash
sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed
```

| Sintoma | Causa provável |
|---|---|
| `OUT_OF_MEMORY` / `OOM` | `mem_mb` baixo para a rule |
| `TIMEOUT` | `runtime_min` baixo |
| `State=FAILED`, exit 1 | erro da ferramenta — veja o log da rule |
| Job em `PD` eterno | partição errada, ou pediu mais recurso do que o nó tem |
| `command not found` | o env não foi ativado no job |

**Após queda ou cancelamento**, o lock do Snakemake fica órfão:

```bash
snakemake --configfile config.yaml config.local.yaml \
          --workflow-profile profiles/slurm --unlock
```

Depois relance normalmente. O `rerun-incomplete` está ligado no profile local;
para o SLURM, acrescente `--rerun-incomplete` se o profile não o trouxer.

**Nunca use `--touch`.** Ver `CLAUDE.md` §4.3.

---

## 7. Onde os arquivos ficam

Igual à execução local — o SLURM não muda isso:

| | |
|---|---|
| Intermediários | `/mnt/scratch-shared/<user>/<project>/` |
| Entregáveis | `results/` no workdir |
| Scratch por job | `/mnt/scratch-shared/jobs/$SLURM_JOB_ID`, removido pelo epilog |

O scratch por job (criado pelo `prolog`) é para temporários da própria
ferramenta. O `work_dir` do workflow é outra coisa e **sobrevive** ao fim do
job — é ele que guarda os BAMs entre as etapas.

> **`/tmp` é tmpfs em RAM** neste host. Ferramentas que usam `TMPDIR` pesado
> devem apontar para o scratch do job.

---

## 8. Verificação ao terminar

Não confie no exit code. Confira conteúdo:

```bash
ls results/counts/*.genes.results | wc -l        # deve bater com as amostras
find results -type f -size 0 ! -path '*/logs/*'  # deve ser vazio
grep -H "Uniquely mapped reads %" results/logs/*_Log.final.out

# @RG utilizável pelo variant calling
samtools view -H results/mapped/<amostra>_Aligned.sortedByCoord.out.bam | grep '^@RG'
```

Espere `ID:`, `SM:`, `LB:` e `PL:` como campos separados por tab. Um `@RG` com
tudo num campo só é o bug corrigido em 2026-10-30 — ver `CLAUDE.md` §4.1.

E o resumo de recursos, para calibrar a próxima rodada:

```bash
sacct -u $USER --starttime today --format=JobName%22,MaxRSS,Elapsed,State \
  | sort -u
```

---

## 9. O que foi aprendido rodando de verdade (2026-10-31)

Três coisas que só apareceram ao executar um lote real, e nenhuma delas aparece
em dry-run.

### 9.1 🔴 Saída declarada de regra NÃO pode ficar em disco local

O Snakemake valida as saídas **a partir do nó orquestrador**. Um arquivo no
`/scratch` de um compute node não existe para o head:

```
MissingOutputException: Job 5 completed successfully,
but some output files are missing. Missing files after 60 seconds.
Removing output files of failed job map_reads since they might be corrupted
```

O group job tinha terminado bem — SLURM reportou `COMPLETED 0:0` nos quatro
passos — e o Snakemake **apagou 48 GB de trabalho válido** por não enxergá-lo.

**Agrupar não resolve.** O job roda num nó só, mas quem confere continua sendo o
head. Custou 1h45 de computação descobrir.

| | |
|---|---|
| Saída declarada de regra | **sempre** no `work_dir` compartilhado |
| Temporário de ferramenta | pode ser local — o Snakemake nem sabe que existe |

### 9.2 Temporário de ferramenta é o maior consumidor de I/O compartilhado

Medido com 2 amostras, antes da correção:

```
32G  counts/SRR11853682.temp     <- temporário do RSEM
29G  counts/SRR11853683.temp
824K counts/SRR11853682_ReadsPerGene.out.tab   <- o entregável
```

**~30 GB por amostra**, contra ~4 GB do FASTQ aparado. STAR e RSEM criam o
temporário **ao lado da saída** quando não se passa `--outTmpDir` /
`--temporary-folder`. O STAR ainda limpa ao terminar, então não deixava rastro —
só a I/O, invisível.

Depois da correção, mesma amostra:

```
compartilhado:  clean 7,9G   mapped 28G   counts 35M   logs 4,1M
local:          tmp 32G  (removido pelas ferramentas ao fim)
```

Nas 73 amostras do L3 são ~2,3 TB que deixam de disputar os 4 fusos do ME5012.

> **Regra geral: medir onde o byte cai, não onde o workflow diz que ele cai.**
> `du` durante a execução revelou o que nenhuma leitura do DAG revelaria.

### 9.3 `slurm-client` tem de existir nos compute nodes

Com agrupamento, o `snakemake` roda **dentro** do nó e invoca cada regra via
`srun`. Sem o pacote, o group job morre em segundos:

```
/bin/sh: 1: srun: not found
```

Nunca apareceu antes porque, sem agrupamento, `srun` só precisa existir no head.
Os R550 vieram só com `slurmd`. A versão do cliente tem de bater com a do
`slurmctld`.

### 9.4 Tempos medidos — base para calibrar o L3

Amostra de ~76M reads (L3, 100 bp), partição `cpu` no R550 (48 CPUs, 380 GB):

| Etapa | Sozinha no nó | 2 amostras competindo |
|---|---|---|
| `trim_reads` (fastp) | 3m19s | 3m05s |
| `map_reads` (STAR) | **14m48s** | 32m49s |
| `quantify_transcripts` (RSEM) | 1h06 | 1h00–1h07 |
| `fix_read_group` | 3s | 3s |
| **total do group job** | **1h24** | 1h33–1h43 |

O STAR mais que dobra de velocidade sem concorrência; o RSEM quase não muda.
Isso sugere limitar jobs simultâneos por nó em vez de maximizar o paralelismo.

`fix_read_group` em 3 segundos é o esperado, não um erro: a regra verifica se o
`@RG` já tem `SM:` e, em caso afirmativo, faz `ln -f` em vez de reescrever.

### 9.5 Verificado no entregável

```
Uniquely mapped reads %  72.27%
genes.results            50.117 linhas, 8 MB
@RG                      ID: SM: LB:TRUSEQ PL:ILLUMINA  (campos separados)
samtools quickcheck      OK
```

> `.export_done` tem **zero bytes de propósito** — é sentinela criada por
> `touch()`. Uma verificação que procure "arquivos de zero bytes" vai acusá-la.

> O BAM exportado **não vem indexado**. É deliberado: o `interomics_snakemake`
> cria os próprios `.bai` nas etapas dele. Para abrir no IGV, indexe à mão.
