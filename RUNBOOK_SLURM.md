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
