# Referência unificada — a divergência entre `rnaseq` e `interomics_snakemake`

> **Status deste documento: PROPOSTA.** As seções 1 e 2 registram fatos
> conferidos em 2026-11-08 (comandos na §6). A recomendação da §4 é um ponto de
> partida para decisão, **não** uma decisão tomada e **nada** aqui foi executado.

Os dois workflows foram escritos contra genomas de referência **diferentes**, e
isso só apareceu quando os BAMs do STAR foram levados ao variant calling. A
conta foi paga em semanas de trabalho extra (Sprint 0 e 2.1 do
`tese_natalia/VARIANT_CALLING.md`). Este documento existe para que a próxima
pessoa não redescubra isso, e para que, se alguém quiser unificar, saiba o que
custa.

---

## 1. A divergência

| | `rnaseq` (este repo) | `interomics_snakemake` |
|---|---|---|
| Origem | NCBI RefSeq **GRCh38.p14** (`GCF_000001405.40`) | GRCh38 "analysis set" (nome do arquivo: `GRCh38_character_genome`) |
| FASTA | `/mnt/reference/rnaseq/GRCh38_p14/fasta/GRCh38_p14.fa` (3,3 GB) | `/mnt/reference/interomics_snakemake/reference_genome/GRCh38_character_genome.fasta` (3,1 GB) |
| Anotação | `.../gtf/GRCh38_p14.gtf` (NCBI, 67.512 genes, 539 contigs com linhas) | não usa GTF para mapear; usa o cache VEP `merged` 116 e `refGene` do ANNOVAR |
| Namespace | **`NC_000001.11`**, `NT_187361.1`, `NW_...` | **`chr1`**, `chr1_KI270706v1_random`, `chrUn_...` |
| Nº de contigs | **705** (`chrName.txt` do índice STAR; 705 `@SQ` no `.dict`) | **195** (`.fai` e `.dict`) |
| Composição | 25 primários + unlocalized/unplaced + **514 alt loci e patches** (`NT_` 189, `NW_` 325) | 25 primários + 42 `_random` + 127 `chrUn` + `chrEBV`. **Sem ALT.** |
| Máscara | sem N extra | **N-mascarado** em 6 primários (chr5, chr14, chr19, chr21, chr22, chrY): 13,9 Mb em 54 regiões onde o NCBI tem base |
| Índices | RSEM + STAR em `rsem_star_sjdb{74,75,99,124,149}/` | `.fai`, `.dict`, BWA (`.amb .ann .bwt .pac .sa`) |
| IDs de gene | `gene_id` do GTF = **símbolo** (`DDX11L1`), com `db_xref "GeneID:100287102"` | — |

Pontos que importam:

- **As sequências são as mesmas** onde os dois existem: 185 dos 191 contigs
  pareados têm `M5` idêntico; os outros 6 diferem só porque o interomics tem `N`
  onde o NCBI tem base. Ou seja, **as coordenadas coincidem** e o problema é de
  *nome, ordem e conjunto de contigs*, não de build.
- **Renomear não basta.** Em BAM o `RNAME` é um índice inteiro na lista `@SQ`.
  Trocar a ordem do cabeçalho sem reescrever os registros desloca coordenadas em
  silêncio; `samtools reheader` não serve.
- **Os recursos do GATK só existem em `chr*`:** `dbsnp_151.vcf.gz`,
  `1000G_phase1.snps.high_confidence.hg38.vcf.gz` e
  `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (em
  `/mnt/reference/interomics_snakemake/resources/gatk_resources/`) trazem
  `chr1…` nos registros. BQSR só é possível **depois** da conformação. O VEP e o
  ANNOVAR (`hg38_refGene.txt`) também são `chr*`.

---

## 2. O que existe hoje só por causa da divergência

### 2.1 Artefatos de dado

Em `/mnt/shared/projetos/tese_natalia/referencia_ncbi/`:

| Arquivo | O que é |
|---|---|
| `rna_contig_map.tsv` | tabela `NC_/NT_/NW_ → chr*`, 2 colunas TAB sem cabeçalho: **191 contigs**, 25 primários, na ordem do `.dict` do interomics |
| `rna_contig_map.relatorio.txt` | contagens, 4 contigs do interomics sem par (`chr11_KI270721v1_random`, `chr22_KI270734v1_random`, `chrUn_KI270752v1`, `chrEBV`), 514 do NCBI sem par, os 6 pares aceitos por mascaramento |
| `GRCh38_p14.dict` (+ `.fa` symlink, `.fa.fai`) | `.dict` do FASTA NCBI (705 `@SQ`). A referência original é somente-leitura, e o GATK exige FASTA+`.fai`+`.dict` no mesmo diretório |
| `interomics_mascarado.bed` | as 54 regiões em que o interomics tem `N` e o NCBI tem base (cruzar com genes de interesse: a callability ali é zero) |
| `intervalos/exons_pad50.{bed.gz,interval_list}` (+ `exons_nopad.bed`, `exons_resumo.txt`) | alvo do HaplotypeCaller/ApplyBQSR: exons do **GTF do NCBI** (contigs `NC_`), renomeados pelo mapa, +50 pb, **sem chrM**: 285.824 intervalos, 176.518.454 pb. Declarado em `refgenome.bed` e `refgenome.interval` |

### 2.2 Scripts neste repo (`misc/`)

Todos rodam por `sbatch`, são idempotentes e validam o resultado (não só o exit code):

| Script | Produz |
|---|---|
| `misc/make-ncbi-dict.sbatch` | `GRCh38_p14.dict`, `.fa.fai` e o symlink |
| `misc/make-contig-map.sbatch` | `rna_contig_map.tsv`, `.relatorio.txt`, `interomics_mascarado.bed`. Pareia por `M5` (camada 1) e, se só difere por `N`, por comprimento único (camada 2); qualquer outra diferença reprova o job |
| `misc/compare-masked-contigs.sbatch` | a investigação que descobriu o mascaramento (base a base nos 6 contigs). Foi diagnóstico; o mapa já embute o resultado |
| `misc/make-exon-intervals.sbatch` | `intervalos/`. Falha se qualquer um dos 10 genes Tier 1 não tiver 100% dos exons no alvo |

### 2.3 No `interomics_snakemake`

- `rules/rna_conform_bam.snakefile` + `scripts/rna_conform_bam.sh`: lê o BAM do
  STAR (705 contigs NCBI) e produz um BAM em 195 contigs `chr*`, na ordem do
  `.dict`. Numa passada: gate de `@RG` (SM, PL, SM == amostra), descarta
  alinhamentos secundários (`-F 0x100`), **descarta os 514 contigs sem par** (199
  Mb), renomeia pelo `contig_map`, reordena (`samtools sort`), indexa e verifica.
- Perde-se informação: nos dados da tese, **8,7% a 19,0% dos registros** caem em
  contigs descartados (~80% deles são secundários). É irreversível sem remapear
  (por isso o BAM original entra somente-leitura).
- Custo por amostra: reescrita completa do BAM, 8 threads, 24 GB
  (`threads.rna_conform`, `ram_memory.rna_conform` no config), mais I/O de leitura
  de um BAM de até ~17 GB e escrita do intermediário `temp()`.
- `rnaseq.contig_map` é um parâmetro **obrigatório** do modo RNA: sem ele o parse falha.
- O config de produção que amarra tudo:
  `tese_natalia/corridas/vc_L1_L3/config_rnaseq_L1.yaml` (`rnaseq.contig_map`,
  `refgenome.{ref,dict,bed,interval}`).

### 2.4 Efeitos colaterais que não são "arquivos"

- **callability:** 114.834 pb de exon (de 146.723.341) caem em regiões
  N-mascaradas do interomics (`exons_resumo.txt`). O STAR mapeou onde havia
  sequência; a chamada de variante, não.
- **o alvo de variant calling vem do GTF do NCBI**, não de um BED do interomics.
  Se a anotação de contagem mudar, o alvo muda junto.
- **quem lê os dois mundos precisa saber o namespace** de cada arquivo: contagens
  e BAMs exportados deste workflow são `NC_`; VCFs, BAMs recalibrados e VEP do
  interomics são `chr*`.

---

## 3. Opções de unificação

### (a) Migrar o `rnaseq` para o FASTA do interomics (`chr*`)

O STAR passaria a mapear contra os mesmos 195 contigs e na mesma ordem do
`.dict`. Por construção, o BAM já sairia na ordem do dict: dispensa contig map,
`GRCh38_p14.dict`, BED de máscara e a reescrita do BAM. O que sobra da
conformação é o gate de `@RG`, descartar secundárias e a verificação — a conferir
com `samtools view -H bam | grep '^@SQ'` contra o `.dict`, antes de confiar. BQSR
rodaria direto.

Precisa de uma anotação **no namespace `chr*`**. Três candidatas:

| | (a1) GENCODE primary assembly | (a2) RefSeq do NCBI com contigs renomeados | (a3) RefSeq "estilo UCSC" pronto |
|---|---|---|---|
| Como | GTF `primary_assembly.annotation` do GENCODE (<https://www.gencodegenes.org/human>) | o mesmo `GRCh38_p14.gtf` de hoje, passado por `awk` com o `rna_contig_map.tsv`, descartando linhas em contigs sem par | `hg38.ncbiRefSeq.gtf.gz` em <https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/> |
| Nomes de contig | `chr*` (confirmado na documentação do GENCODE: FASTA e GTF compartilham os nomes) | `chr*` pelo mapa (já validado por `M5`) | `chr*` (UCSC); **não verifiquei o conteúdo do arquivo** |
| IDs de gene | **muda:** `ENSG...` (com versão) em vez de símbolo | **igual** ao de hoje (símbolo + `GeneID` no `db_xref`) | gerado por `genePredToGtf`; semântica de `gene_id` **não verificada** |
| Anotação | outra: Ensembl/HAVANA, mais lncRNA | **a mesma** dos lotes L1–L4 | RefSeq alinhado pelo UCSC; relato (não verificado hoje) de que o chrM ficou fora no hg38 |
| Coerência com o resto da tese | o variant calling usa alvo e ANNOVAR RefSeq, então contagem e variante ficariam em anotações diferentes | mantém | mantém, mas é derivado diferente do atual |
| Artefato a manter | o GTF baixado | **um script nosso** que gera o GTF (rastreável, com md5) | o GTF baixado |

**Impacto em `ETL_gene_info.R` e `DEG.R` (confirmado lendo o código):**

- `ETL_gene_info.R` gera `gene_id` (GeneID Entrez), `gene_symbol` e `synonyms`
  do `gene_info` do NCBI. `DEG.R` junta o resultado do DESeq2 com ele por
  `gene_symbol == gene_symbol`, porque o `rownames` do `ReadsPerGene.out.tab` do
  STAR **é o `gene_id` do GTF — hoje, o símbolo**. O Entrez sai dessa junção e
  alimenta `ReactomePA::enrichPathway`.
- Com **(a1)** o `rownames` viram `ENSG...`. A junção por símbolo devolve `NA`
  para tudo, e o Reactome roda com Entrez vazio. Seria preciso reescrever
  `DEG.R` e o ETL para mapear ENSG → Entrez (por exemplo via `org.Hs.eg.db`, que
  já está no `environment.yaml`), tratando versão do ENSG, duplicatas `_PAR_Y` e
  ENSG sem Entrez. O RSEM com `--append-names` também muda o formato do ID.
  Qualquer análise a jusante indexada por símbolo precisa ser revisada.
- Com **(a2)** nada disso muda: mesmos `gene_id`, mesmo `gene_info`, mesmo
  `DEG.R`. É a única opção que troca o genoma **sem** trocar a anotação.
- Com **(a3)** é preciso abrir o arquivo e ver o que ele põe em `gene_id` antes
  de decidir; não está confirmado.

**Impacto nas contagens (hipótese, não medida):** sem os 514 contigs ALT/patch,
read que hoje cai em mais de um lugar (primário + ALT) passaria a ser única.
`N_multimapping` do STAR cairia e as contagens de genes com cópia em ALT
(MHC/HLA, KIR, entre outros) tendem a subir. Parte dos 8,7% a 19,0% de registros
descartados pela conformação são exatamente esses. O N-mascaramento também
mexe em PAR do chrY e nas duplicações falsas do GRCh38. Isso é, em tese, *melhor*
para RNA-seq (é prática usual mapear em assembly sem ALT), mas **muda os
números**: qualquer DEG que misture lotes de referências diferentes tem a
referência confundida com o lote.

**Custo de recomputar (se aplicado aos lotes existentes):**

| Item | O que precisaria refazer |
|---|---|
| Índices | 1 índice STAR+RSEM novo por `read_length` (hoje 5: 74, 75, 99, 124, 149), ~33 GB cada (valor medido no sjdb149) |
| Mapeamento | todos os BAMs (L1: 86, L2: 8, L3: 73, L4: 62 amostras). Referência de ordem de grandeza: o L3 (73 amostras) levou **10h56** nos quatro nós |
| Contagens e DEG | refazer tudo; os números mudam |
| Variant calling | os BAMs mudam, então HaplotypeCaller, BQSR, filtro e VEP de todas as amostras seriam refeitos (L1+L3 somaram 3.106 jobs) |
| Código | `reference.smk` baixa por `wget` de URL; apontar para um FASTA local exige mudar `reference.smk`/`config.yaml` (ex.: `reference.genome_path`). Não existe hoje |

### (b) Manter os dois e formalizar a camada de conformação

Nada é remapeado, nada muda de número. O custo é permanente e pequeno: a camada
de §2 continua existindo para qualquer lote novo. Formalizar significa:

- tratar `referencia_ncbi/` como artefato versionado (md5 de `rna_contig_map.tsv`
  e `GRCh38_p14.dict` registrados; hoje só o `MANIFEST.txt` do índice existe);
- manter `misc/make-*.sbatch` e `rna_conform_bam` como a forma suportada e
  testada de ligar os dois workflows;
- documentar (feito aqui, e no README/RUNBOOK) que o BAM exportado é `NC_` e
  que a conformação descarta 8,7–19,0% dos registros.

Riscos que permanecem: depende de `M5` bater (se o FASTA do interomics mudar,
o `contig_map` precisa ser refeito, e o script reprova em vez de passar calado);
cada lote paga a reescrita do BAM; e alguém pode esquecer o namespace de um
arquivo.

### (c) Migrar o `interomics_snakemake` para o namespace NCBI

**Não recomendado.** Os recursos do GATK (dbSNP, 1000G, Mills), o cache VEP, o
`refGene` do ANNOVAR, os `interval_files`, o índice BWA e o `.dict` estão em
`chr*`. Teria de converter cada VCF de recurso (`bcftools annotate --rename-chrs`
com o mesmo mapa), refazer BWA e dict, reconstruir os intervalos e revalidar o
VEP/ANNOVAR — e isso quebraria os workflows de DNA (WES/WGS) que são o uso
principal do interomics. Seria um segundo conjunto de recursos para manter, pior
do que a camada de (b).

---

## 4. Recomendação (PROPOSTA)

1. **Agora, nos lotes L1–L4: (b).** Congele o NCBI + conformação. Os BAMs,
   contagens e VCFs já existem e foram validados; remapear muda números,
   obriga a refazer DEG e variant calling, e não destrava nada que a camada
   atual não já entregue.
2. **Para trabalhos futuros que NÃO serão combinados com L1–L4: (a2).** Gere
   um índice novo (`genome_id` diferente, ex. `GRCh38_p14_chr`) com o FASTA do
   interomics e o **GTF do NCBI renomeado** pelo `rna_contig_map.tsv`. É a única
   variante que preserva os `gene_id`, o `DEG.R`, o ETL e a coerência RefSeq com
   ANNOVAR e com o alvo de exons. Não use GENCODE (a1) sem decidir antes mudar
   toda a anotação da tese. Não adote (a3) sem abrir o arquivo.
3. **Regra de ouro:** nunca combine, num mesmo DEG ou coorte de variantes,
   amostras mapeadas em referências diferentes. Um lote novo que entre numa
   análise junto com L1–L4 **tem de** usar o NCBI de hoje.
4. **(c) está descartado.**

### O que NÃO deve ser refeito agora

- Remapear L1–L4 ou refazer contagens e DEG.
- Refazer `rna_contig_map.tsv`, `GRCh38_p14.dict`, `interomics_mascarado.bed` e
  `exons_pad50.*`. Estão validados e são entrada do variant calling em curso.
- Mexer nos recursos do `interomics_snakemake`.

### Antes de adotar (a2) para um projeto futuro

1. Escrever (e versionar) o script que gera o GTF `chr*` a partir do GTF do NCBI
   e do `rna_contig_map.tsv`, descartando linhas em contigs ausentes do FASTA
   (RSEM e STAR reclamam de contigs do GTF fora do FASTA; o comportamento exato de
   cada um não foi testado aqui).
2. Permitir FASTA/GTF locais em `reference.smk` (hoje só URL).
3. **Piloto de 2 a 3 amostras** já mapeadas no NCBI (por exemplo uma de L2 e uma
   de L3): remapear no índice novo e comparar (i) correlação das contagens por
   gene, (ii) `N_multimapping` e % de únicos, (iii) lista de genes cuja contagem
   mudou mais de X% — definindo X antes —, (iv) variantes nos genes Tier 1.
   Só então decidir.
4. Conferir `samtools view -H` do BAM novo contra o `.dict` do interomics.

---

## 5. O que fazer ao encontrar o namespace "errado"

| Sintoma | Causa |
|---|---|
| GATK: `Input files reference and reads have incompatible contigs` | BAM `NC_` contra `.dict` `chr*` (ou o inverso). Passe pela `rna_conform_bam` |
| `bcftools`/`tabix`: "contig not found" ao cruzar VCF de recurso com BAM/BED | namespace misturado. Use `rna_contig_map.tsv` para traduzir (col. 1 → col. 2) |
| Genes do alvo "somem" no variant calling | estão em região N-mascarada do interomics: confira `interomics_mascarado.bed` |

---

## 6. Como conferir cada afirmação

```bash
# interomics: 195 contigs, todos chr*, sem ALT
F=/mnt/reference/interomics_snakemake/reference_genome/GRCh38_character_genome.fasta
wc -l < $F.fai; head -3 $F.fai; cut -f1 $F.fai | grep -vc '^chr'; grep -c '^@SQ' ${F%.fasta}.dict

# rnaseq: 705 contigs, namespace NC_, IDs de gene = símbolo
G=/mnt/reference/rnaseq/GRCh38_p14
wc -l < $G/rsem_star_sjdb149/chrName.txt; head -2 $G/rsem_star_sjdb149/chrName.txt
grep -v '^#' $G/gtf/GRCh38_p14.gtf | awk '$3=="gene"' | head -1 | cut -f1,9 | cut -c1-120

# camada de conciliação
R=/mnt/shared/projetos/tese_natalia/referencia_ncbi
wc -l $R/rna_contig_map.tsv; cat $R/rna_contig_map.relatorio.txt; wc -l $R/interomics_mascarado.bed
cat $R/intervalos/exons_resumo.txt

# recursos do GATK em chr*
zcat /mnt/reference/interomics_snakemake/resources/gatk_resources/dbsnp_151.vcf.gz | head -n 4000 | grep -v '^#' | head -1 | cut -f1,2
```

Documentos relacionados: `README.md` (raiz), `RUNBOOK.md`, `RUNBOOK_SLURM.md`;
no outro lado, `tese_natalia/VARIANT_CALLING.md` (§0.3, 0.4, 2.1) e
`interomics_snakemake/rules/rna_conform_bam.snakefile`.
