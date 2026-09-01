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
<html style="margin:0;padding:0;">
<div style="background-color:transparent;box-sizing:border-box;">
  <div style="padding:6px 2px;font-family:'Roboto','PingFang SC','Segoe UI',Arial,sans-serif;--accent:#A3D5E8;box-sizing:border-box;">
    <div style="font-size:15px;font-weight:600;color:#1A1B1C;margin-bottom:2px;">制药领域关心的非经典激酶（atypical kinase）</div>
    <div style="font-size:12px;color:#6B7280;margin-bottom:12px;">催化域都不属于 PF00069/PF07714。颜色=制药热度（红=高，橙=中，灰=低）。</div>

    <div style="display:flex;gap:10px;flex-wrap:wrap;align-items:stretch;">
      <!-- PIKK -->
      <div style="flex:1 1 250px;min-width:0;padding:10px 12px;background:linear-gradient(135deg,rgba(226,74,104,0.18),rgba(226,74,104,0.34));border-radius:12px;box-sizing:border-box;">
        <div style="font-size:12px;font-weight:600;color:#1A1B1C;">PIKK 家族（热度高）</div>
        <div style="font-size:11px;color:#1A1B1C;line-height:1.55;margin-top:3px;font-family:Consolas,monospace;">mTOR, ATM, ATR,<br>DNA-PKcs, SMG1, TRRAP*</div>
        <div style="font-size:11px;color:#1A1B1C;margin-top:4px;">催化域=PF00454；特征域 FAT(PF02259)+FATC(PF02260)。<br>*TRRAP 无催化活性（假激酶）</div>
      </div>
      <!-- PI3K -->
      <div style="flex:1 1 250px;min-width:0;padding:10px 12px;background:linear-gradient(135deg,rgba(226,74,104,0.18),rgba(226,74,104,0.34));border-radius:12px;box-sizing:border-box;">
        <div style="font-size:12px;font-weight:600;color:#1A1B1C;">PI3K / PI4K 脂质激酶（热度高）</div>
        <div style="font-size:11px;color:#1A1B1C;line-height:1.55;margin-top:3px;font-family:Consolas,monospace;">PIK3CA/B/D/G, VPS34,<br>PI4KIIIα/β, PI4KII</div>
        <div style="font-size:11px;color:#1A1B1C;margin-top:4px;">催化域=PF00454 (PI3_PI4_kinase)。重磅靶点：alpelisib/idelalisib 等已上市。</div>
      </div>
      <!-- SphK -->
      <div style="flex:1 1 250px;min-width:0;padding:10px 12px;background:linear-gradient(135deg,rgba(233,183,143,0.28),rgba(233,183,143,0.5));border-radius:12px;box-sizing:border-box;">
        <div style="font-size:12px;font-weight:600;color:#1A1B1C;">鞘氨醇激酶（热度中）</div>
        <div style="font-size:11px;color:#1A1B1C;line-height:1.55;margin-top:3px;font-family:Consolas,monospace;">SPHK1, SPHK2</div>
        <div style="font-size:11px;color:#1A1B1C;margin-top:4px;">EC 2.7.1.91，S1P 信号通路靶点。Pfam accession 待查（建议用基因名白名单）。</div>
      </div>
      <!-- Alpha -->
      <div style="flex:1 1 250px;min-width:0;padding:10px 12px;background:linear-gradient(135deg,rgba(233,183,143,0.28),rgba(233,183,143,0.5));border-radius:12px;box-sizing:border-box;">
        <div style="font-size:12px;font-weight:600;color:#1A1B1C;">alpha-kinase（热度中）</div>
        <div style="font-size:11px;color:#1A1B1C;line-height:1.55;margin-top:3px;font-family:Consolas,monospace;">EEF2K, ALPK1/2/3,<br>TRPM6, TRPM7</div>
        <div style="font-size:11px;color:#1A1B1C;margin-top:4px;">独立折叠，与 ePK 无同源。Pfam≈PF02816（建议验证）。eEF2K 肿瘤代谢靶点。</div>
      </div>
      <!-- RIO -->
      <div style="flex:1 1 250px;min-width:0;padding:10px 12px;background:rgba(255,255,255,0.7);border:1px solid rgba(0,0,0,0.08);border-radius:12px;box-sizing:border-box;">
        <div style="font-size:12px;font-weight:600;color:#1A1B1C;">RIO 激酶（热度低）</div>
        <div style="font-size:11px;color:#1A1B1C;line-height:1.55;margin-top:3px;font-family:Consolas,monospace;">RIOK1/2/3</div>
        <div style="font-size:11px;color:#1A1B1C;margin-top:4px;">核糖体生物发生相关。Pfam≈PF01163（建议验证）。</div>
      </div>
    </div>
  </div>
</div>
</html>
