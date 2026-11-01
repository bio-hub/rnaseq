# RUNBOOK — rodar um lote via SLURM

Complementa `RUNBOOK.md`, que cobre a execução local. Leia aquele primeiro: as
proibições, a verificação por conteúdo e as armadilhas de dado valem igual aqui.

---

## 1. Leia isto antes de decidir usar SLURM

> **Atualizado em 2026-11-01.** Até 31/10 esta seção dizia que não havia NFS e
> que o SLURM era "uma fila, não um cluster". **Deixou de ser verdade.** Se você
> encontrar essa afirmação em outro documento, é texto velho.

O cluster tem os quatro nós em produção, com **288 núcleos**:

| Partição | Nós | Núcleos | RAM | Tier 0 (`/scratch`) |
|---|---|---:|---:|---|
| `bigmem` | `r660xs` (head) | 128 | 1003 GB | **HDD** — ver §1.1 |
| `gpu` | `r760` | 64 | 510 GB | NVMe |
| `cpu` | `r550-1`, `r550-2` | 48 cada | 380 GB | SSD |

`/mnt/shared`, `/mnt/reference`, `/mnt/scratch-shared` e `/mnt/home` são
exportados por NFSv4 pelo head e montados em todos os workers. Um job em
qualquer nó enxerga a referência, os FASTQ e o conda compartilhado.

O `profiles/slurm/config.yaml` usa `slurm_partition: "bigmem,cpu,gpu"` — o
SLURM aloca em qualquer uma.

**Validado em produção:** o lote L3 rodou 73 amostras, 9.363 M reads, 295/295
passos, em 10h56 distribuídas pelos quatro nós.

### 1.1 O que esperar de desempenho, e o que não esperar

Com quatro nós em vez de um, o ganho fim a fim medido foi de **1,16×** — não
de 2,25×, que é a razão de núcleos. Até o fim do alinhamento o ganho era de
1,96×; o resto se perdeu numa cauda de 4h29 em que o cluster terminava o RSEM
de um punhado de amostras com três máquinas ociosas.

A causa é o `/scratch` do head ser HDD, e ela produziu **8 jobs mortos por
TIMEOUT** numa única corrida, todos nesse nó. Mitigado em 01/11 montando
`/scratch/ronald` em tmpfs no head. Diagnóstico completo e a ordem dos ajustes
pendentes estão em `superpc_lika/DESEMPENHO.md` — **leia antes de propor
qualquer ajuste de desempenho**, porque dois dos ajustes óbvios pioram as
coisas se aplicados na ordem errada.

**O que você ganha com o SLURM, além do paralelismo:**

- A fila não depende de uma sessão `tmux` viva
- Contabilidade por job (`sacct`): tempo, memória de pico, exit code
- Limite de concorrência por recurso, não por "quantos couberam"
- O `prolog`/`epilog` cria e remove scratch isolado por job
- `restart-times` recupera job morto por TIMEOUT — salvou 5 amostras no L3

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

O diretório do projeto vai em **`/mnt/shared/projetos/`**, nunca em `~`:

```bash
mkdir -p /mnt/shared/projetos/hs_rnaseq_L2 && cd /mnt/shared/projetos/hs_rnaseq_L2
git clone -b main https://github.com/bio-hub/rnaseq.git .
mkdir -p samples metadata
```

> **Por que não `~/projetos/`.** O home do administrador é **local em cada nó**
> — `/home/ronald` no head não é o `/home/ronald` do `r550-1`. Um projeto ali
> fica invisível aos workers, fora do backup e preso a um disco só.
>
> Isto não é hipotético: a versão anterior deste runbook mandava
> `mkdir -p ~/projetos/hs_rnaseq_L2`, e foi exatamente assim que os 135 GB de
> entregáveis do L2 foram parar no disco local do head, onde ficaram até
> 01/11. Tiveram de ser copiados de volta, e a raiz do head estava em 87%.
>
> `/mnt/home/<tenant>` serve para usuários tenant; projetos e entregáveis vão
> para `/mnt/shared/projetos/`, que é legível por todo o grupo `bioinfo`.

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
  conda activate biohub_rnaseq && cd /mnt/shared/projetos/hs_rnaseq_L2 && \
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
               conda activate biohub_rnaseq && cd /mnt/shared/projetos/hs_rnaseq_L2 && \
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

---

## 9. Lote 3 — configuração medida

Medido em 2026-10-31 com uma amostra do L3 (`SRR11853682`, 76M reads) sozinha
no `r550-1`:

| Etapa | Threads | Tempo | Observação |
|---|---|---|---|
| `trim_reads` (fastp) | 8 | 3m19s | |
| `map_reads` (STAR) | 16 | 14m48s | 30+ min quando duas amostras competiam |
| `quantify_transcripts` (RSEM) | 12 | **1h06m** | o gargalo — 79% do tempo total |
| `fix_read_group` | 4 | 3s | hardlink; o STAR já grava `@RG` válido |
| **total** | | **1h24m** | |

Resultado: 72,27% de mapeamento único, `genes.results` com 50.117 linhas,
`@RG` com campos separados por tab.

### 9.1 Preparar o diretório

O workdir precisa estar em **`/mnt/shared`**, não em `~`: os homes são locais a
cada nó e divergentes (ver plano §17.5). Um job num compute node não enxerga
nada de `~` do head.

```bash
mkdir -p /mnt/shared/projetos/hs_rnaseq_L3 && cd /mnt/shared/projetos/hs_rnaseq_L3
git clone -b main https://github.com/bio-hub/rnaseq.git .
mkdir -p samples
ln -s /mnt/shared/projetos/hs_rnaseq/samples/*.fastq.gz samples/
```

> Se os FASTQ dos três lotes estiverem no mesmo diretório, faça o symlink só
> dos do L3 — senão o workflow processa 167 amostras em vez de 73.

### 9.2 Config do lote

**Nunca edite o `config.yaml` versionado.**

```yaml
# config.local.yaml
star:
  read_length: 100          # L3 → sjdbOverhang 99
output:
  project: "hs_rnaseq_L3"
  local_scratch: "/scratch/{user}/{project}"
  keep_bams: "coord"
deg_analysis:
  enabled: false            # só quando os três lotes fecharem
```

`local_scratch` move ~32 GB por amostra de temporário do RSEM e do STAR para o
disco local — ~2,3 TB que deixam de atravessar o array de 4 fusos. **Afeta
apenas temporários de ferramenta**, nunca saída declarada de regra: ver o
comentário em `common.smk` sobre a `MissingOutputException`.

### 9.3 Submeter

```bash
tmux new-session -d -s l3
tmux send-keys -t l3 'source /mnt/shared/conda/etc/profile.d/conda.sh && \
  conda activate biohub_rnaseq && cd /mnt/shared/projetos/hs_rnaseq_L3 && \
  snakemake --configfile config.yaml config.local.yaml \
  --workflow-profile profiles/slurm 2>&1 | tee log_L3_$(date +%Y%m%d_%H%M%S).log' C-m

sleep 20 && tmux capture-pane -p -t l3 | tail -15
```

> `tmux new-session` **falha se a sessão já existir**, e o `send-keys` seguinte
> vai para a sessão velha. O sintoma é "nada aconteceu".

### 9.4 Quanto deve demorar

O profile usa `bigmem` (r660xs, 128 CPUs), onde tudo é acessado localmente por
iSCSI — sem NFS no caminho. O `map_reads` pede 16 CPUs e 70 GB, então cabem 8
simultâneos; as demais etapas se sobrepõem entre amostras.

Estimativa grosseira: **15 a 25 horas** para as 73 amostras. Não é linear —
há contenção de I/O, e o número real depende de quantas amostras estão em RSEM
ao mesmo tempo.

Para acrescentar os R550 quando o `bigmem` saturar, troque no profile:

```yaml
default-resources:
  slurm_partition: "bigmem,cpu"
```

O SLURM tenta na ordem. Antes de o `bigmem` encher, espalhar só acrescenta
tráfego de NFS.

### 9.5 O que conferir ao fim

```bash
ls results/counts/*.genes.results | wc -l          # deve ser 73
find results -type f -size 0 ! -name '.export_done'  # deve ser vazio
grep -h "Uniquely mapped reads %" results/logs/*_Log.final.out
samtools view -H results/mapped/<amostra>*.bam | grep '^@RG'
```

> `.export_done` tem zero bytes **por design** — é sentinela criada por
> `touch()`. Excluí-la da busca evita falso alarme.

`MaxRSS` virá vazio no `sacct`: `jobacct_gather/cgroup` com cgroup v2 não
popula RSS neste cluster. Os limites de memória **funcionam** (job que estoura
é OOM-killed); o que falta é só a contabilidade. Ver plano §17.6.
