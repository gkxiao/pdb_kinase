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
···


