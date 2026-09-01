# pdb_kinase
## 1. 含有Kinase结构域的PDB code列表
SIFTS 是 PDBe 提供的 PDB↔UniProt↔Pfam 残基级映射，`pdb_chain_pfam.tsv.gz` 直接给出**每条 PDB 链命中的 Pfam 域及残基范围**
```bash
# 下载 SIFTS 的链级 Pfam 注释（路径以 quick.html 页面链接为准）
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz
zcat pdb_chain_pfam.tsv.gz | awk -F'\t' '$4=="PF00069" || $4=="PF07714"'   > kinase_chains.tsv
```
用到的关键 accession：

- **PF00069（Pkinase）**：丝氨酸 / 苏氨酸激酶催化域
- **PF07714（Pkinase_Tyr）**：酪氨酸激酶催化域

输出每行含 PDB ID、链、Pfam accession、**start/end 残基**—— 这就是 "仅限 kinase 结构域" 的精确边界。注意 SIFTS 的域范围是 **UniProt 编号**，落到 PDB 残基要经过 SIFTS 的残基映射换算（`pdb_chain_uniprot.tsv` 提供 UniProt↔PDB 残基对应）。

```bash
# 总的条目数
wc -l pdb_chain_pfam.tsv
1061318 pdb_chain_pfam.tsv

# 激酶条目数
wc -l kinase_chains.tsv
12590 kinase_chains.tsv

# 预览
head kinase_chains.tsv
10bl    A       Q4E2L0  PF00069 1
10bl    B       Q4E2L0  PF00069 1
10dj    A       P06241  PF07714 1
10dj    B       P06241  PF07714 1
10hm    A       Q382U0  PF00069 1
10hm    B       Q382U0  PF00069 1
10hm    C       Q382U0  PF00069 1
10hm    D       Q382U0  PF00069 1
10hm    E       Q382U0  PF00069 1
10hm    F       Q382U0  PF00069 1
```
输出格式是 `PDB<TAB>CHAIN<TAB>SP_PRIMARY<TAB>PFAM_ID<TAB>COVERAGE`。
其中`SP_PRIMARY` 列就是 **UniProt accession**，可以作为PDB/UniProt 的注释，同时解决 "识别 kinase 链"。

## 2. 将含激酶的链转为含激酶的PDB code
```bash
# 结构级白名单：含 kinase 链的 PDB 列表（判定单位 = 结构）
awk -F'\t' '{print $1}' kinase_chains.tsv | sort -u > kinase_pdbs.txt
wc -l kinase_pdbs.txt
`8024 kinase_pdbs.txt`

# 这些结构的全部链注释（含 non-kinase 链，header 是第 2 行）
# 建库时整结构保留，这张表用于后续标注"哪条链是 kinase 链"
zcat pdb_chain_pfam.tsv.gz \
  | awk -F'\t' 'NR==FNR{k[$1]=1;next} FNR==2 || k[$1]' kinase_pdbs.txt - \
  > kinase_structures_allchains.tsv

# 与本地*.pdb 比对，得到真正属于 kinase 子集的文件
ls *.pdb | sed 's/\.pdb//' | sort > local_pdbs.txt
comm -12 local_pdbs.txt kinase_pdbs.txt > local_kinase.txt
wc -l local_kinase.txt
cat local_kinase.txt
```
现在得到了含有经典激酶结构域的PDB结构列表，注意，这并不包含非经典激酶。

## 3. "kinase" 的定义决定库的边界：

- **经典真核蛋白激酶（ePK）**：PF00069 + PF07714 全覆盖（S/T 和 Y 激酶）。
- **非典型激酶（atypical）**：PIKK 家族（ATM/ATR/DNA-PK）、PI3K、alpha-kinase 等的催化域**不属于** Pkinase 折叠，SCOP/CATH 归到不同超家族。若只做经典激酶抑制剂研究，PF00069/PF07714 正好；若想含 PI3K 等，要额外加对应 accession（如 PIKK）。

人类非经典激酶数量有限，基因名单见附件`atypical_genes.txt`

## 4. 非典型激酶
- **PIKK 家族**（DNA-PKcs、ATM、ATR、mTOR、SMG1、TRRAP）与 PI3K 属于同一个 "PI3K 超家族"：保守 C 端催化域与脂质 PI3K 同源，但磷酸化 Ser/Thr*
- PIKK 催化域在 **Pfam 里标注为 PI3_PI4_kinase（PF00454）**，并额外带 **FAT（PF02259）+ FATC（PF02260）** 特征域 —— 所以用 PF00454 能一次覆盖 PI3K/PI4K/PIKK 全部催化域，再用 FAT/FATC 把 PIKK 单独区分出来
- **alpha-kinase 家族**（eEF2K、ALPK1-3、TRPM6/7）：催化域与经典 ePK 无序列同源性，是独立折叠
- **RIO 激酶**（RIOK1-3）：RIO 域是 ePK 催化域的 "修剪版"，常归入非典型激酶
```bash
# 非经典激酶链级筛选（核心：PF00454 一次覆盖 PI3K/PI4K/PIKK）
zcat pdb_chain_pfam.tsv.gz \
  | awk -F'\t' '$4=="PF00454" || $4=="PF02816" || $4=="PF01163"' \
  > atypical_chains.tsv
wc -l atypical_chains.tsv

871 atypical_chains.tsv

# 结构级白名单（独立列表）
awk -F'\t' '{print $1}' atypical_chains.tsv | sort -u > atypical_pdbs.txt
wc -l atypical_pdbs.txt

673 atypical_pdbs.txt

# 从里面再分出 PIKK（带 FAT/FATC 特征域的就是 PIKK，不是脂质激酶）
zcat pdb_chain_pfam.tsv.gz \
  | awk -F'\t' '$4=="PF02259" || $4=="PF02260"{print $1}' \
  | sort -u > pikk_pdbs.txt
wc -l pikk_pdbs.txt

191 pikk_pdbs.txt
```

## 5. 人类激酶
我们对人类的激酶感兴趣，因此需要聚焦。
```bash
# 下载链级 taxonomy 文件（先看列结构）
wget https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz
zcat pdb_chain_taxonomy.tsv.gz | head -3

# kinase 链的 (PDB, CHAIN) 键
awk -F'\t' '{print $1"\t"$2}' kinase_chains.tsv | sort -u > kinase_pdb_chain.txt

# 统计 kinase 链的物种分布（看 9606 人类排第几）
zcat pdb_chain_taxonomy.tsv.gz \
  | awk -F'\t' 'NR==FNR{key[$1 FS $2]=1;next} key[$1 FS $2]{print $3}' \
    kinase_pdb_chain.txt - | sort | uniq -c | sort -rn | head -20

# 人类(9606) kinase 链 → 键列表
zcat pdb_chain_taxonomy.tsv.gz \
  | awk -F'\t' 'NR==FNR{key[$1 FS $2]=1;next} key[$1 FS $2] && $3=="9606"{print $1 FS $2}' \
    kinase_pdb_chain.txt - > human_kinase_pdb_chain.txt
wc -l human_kinase_pdb_chain.txt

# 人类 kinase 链去重后的 UniProt 蛋白数（真正的"人类激酶被 PDB 覆盖数"）
awk -F'\t' 'NR==FNR{key[$1 FS $2]=1;next} key[$1 FS $2]{print $3}' \
  human_kinase_pdb_chain.txt kinase_chains.tsv | sort -u | wc -l

`309`
```

PDB 里被覆盖的人类经典激酶（按 UniProt 去重）是 309 个蛋白：

- 12,471 条去重 kinase 链里，人类（9606）占 10,669 条（85.5%）—— 激酶结构高度集中在人类，符合它作为最热门药物靶标家族的事实。
- 人类链去重到 **309 个 UniProt 蛋白**。对照人类经典 ePK 总数 **478**（Manning 2002），即 **PDB 已覆盖约 65%** 。
- **309 里可能混入 "假激酶"（pseudokinase）**：PF00069/PF07714 是序列域判定，含 Pkinase 域但催化活性缺失的蛋白也会进来（如 HER3/ERBB3、STRAD、VRK3、ILK、KSR1/2 等）。所以 "有激酶活性的 ePK" 覆盖率实际略低于 65%。
- **冗余度极高**：10,669 人类链 / 309 蛋白 ≈ **每个蛋白平均 34 条结构**—— 比之前全库的 25 还高。

## 6. 聚焦人类激酶的库

```bash
# ① 人类激酶结构级 PDB 数（你要下载多少个人类激酶结构）
awk -F'\t' '{print $1}' human_kinase_pdb_chain.txt | sort -u | wc -l

6884
# ② 人类激酶 UniProt 列表（后续打 is_classical/is_pseudokinase 标签用）
awk -F'\t' 'NR==FNR{key[$1 FS $2]=1;next} key[$1 FS $2]{print $3}' \
  human_kinase_pdb_chain.txt kinase_chains.tsv | sort -u > human_kinase_uniprot.txt
wc -l human_kinase_uniprot.txt

309
```
<img src="https://github.com/gkxiao/pdb_kinase/blob/main/the-size-of-human-kinase-in-pdb-202608.png" align='middle'>

