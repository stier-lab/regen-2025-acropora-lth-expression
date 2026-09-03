# LTH RNA-seq literature search log

Updated: 2026-09-02

Question: what coral temperature-expression and regeneration-expression papers
should inform candidate genes, pathway scores, and interpretation for the LTH
heat x wound Acropora pulchra RNA-seq analysis?

This is a documented working search, not a PRISMA-complete systematic review.
It is comprehensive enough to support the current RNA-seq planning files, with
explicit follow-up items for a manuscript-grade literature audit.

## Scope

| Component | Included |
|---|---|
| Population | Reef-building corals, prioritized as Acropora, then other scleractinians |
| Exposure | Temperature stress, sub-bleaching heat, bleaching heat, injury, wound healing, regeneration |
| Outcome | Host gene expression, RNA-seq, scRNA-seq, qPCR biomarkers, gene-expression modules, expression-linked physiology |
| LTH link | Heat x wound, tissue closure, skeletal regeneration, calicoblast/biomineralization, physiology, symbionts, source-patch/genetic effects |

## Local searches run first

Searches used the local Stier Lab literature database (`python3
~/literature-db/tools/litdb ...`) before web expansion.

| Query | Main local hits | Decision |
|---|---|---|
| `coral regeneration gene expression Acropora` | Granados-Cifuentes 2013; Xu 2023; Bay 2009; van de Water 2015; Hemond 2014; Conn 2025; Bay 2015 | Kept Xu, van de Water, Granados, Conn, Bay as relevant |
| `coral wound healing transcriptome regeneration` | Han 2025; Lock 2022; Levy 2021; Mass 2016; Chuang 2021; cnidarian outgroups | Kept Han and Lock; cell atlas/development papers are supporting context only |
| `frontloading coral gene expression` | Barshis 2013 | Kept as the frontloading anchor |
| `transcriptomic resilience symbiont shuffling coral heat gene expression` | Thomas 2019 | Kept as recovery/module/symbiont anchor |
| `Barshis genomic basis coral resilience` | Barshis 2013 | Kept |
| `gene expression coral natural bleaching Seneca` | Seneca and Palumbi 2009/2010 | Kept |
| `Dixon genomic determinants heat tolerance coral` | Dixon et al. 2015 | Kept |
| `Cleves Tinoco heat shock transcription factor coral CRISPR` | Cleves et al. 2020 | Kept |
| `Traylor Knowles heat stress epithelial integrity Acropora` | Traylor-Knowles 2019 | Kept as tissue-barrier context |
| `Exposure elevated sea-surface temperatures below bleaching threshold impairs coral recovery regeneration injury` | Bonesso et al. 2017 | Kept as the closest heat x injury analogue |
| `meta-analysis coral environmental stress response Acropora opposing responses stress intensity` | no local hit | Added by web search |
| `cnidarian biomarkers temperature stress CuZn-SOD c-type lectin FGFR1 MMP SLC26` | no local hit | Added by web search |

## Web expansion

Web searches were used because 2025-2026 coral expression papers may not yet be
in the local database.

| Query theme | New or refreshed source | Decision |
|---|---|---|
| Coral regeneration, wound healing, Acropora, transcriptome | Xu 2023 Frontiers; Han 2025 Communications Biology; Ren 2026 Aquaculture Reports | Xu/Han are core; Ren is supporting and needs local PDF acquisition if used in manuscript |
| Acropora heat stress, frontloading, transcriptomic recovery | Stick et al. 2025 Coral Reefs; Stick et al. 2026 Ecology and Evolution | Kept as recent Acropora recovery/dampening evidence |
| Acropora heat gene-expression meta-analysis | Dixon et al. 2020 Molecular Ecology | Kept as stress-intensity guardrail |
| Cnidarian/coral stress biomarkers | Louis et al. 2017; Molinari et al. 2025 | Kept as review-level candidate-gene guardrails |
| Acropora palmata heat transcriptomics | DeSalvo et al. 2010 | Kept as classic Acropora heat-expression support |
| Global Acropora heat-adaptation genomics | Selmoni et al. 2025/2026 Nature Communications | Kept as genomic association context only |
| Acropora downingi targeted heat biomarkers | Javid et al. 2025 Marine Environmental Research | Kept as low-priority targeted-qPCR context |

## Included evidence set

| Source | Local status | Why included |
|---|---|---|
| Xu et al. 2023, doi:10.3389/fevo.2022.979278 | Local PDF queried | Direct Acropora wound/regeneration RNA-seq |
| Han et al. 2025, doi:10.1038/s42003-025-08089-6 | Local PDF queried | Direct Acropora regeneration scRNA-seq plus bulk RNA-seq |
| van de Water et al. 2015, doi:10.1111/mec.13257 | Local PDF queried | Direct Acropora injury immune time course |
| van de Water et al. 2015, doi:10.1007/s10750-015-2243-z | Local PDF queried | Direct Acropora heat x injury immune response |
| Bonesso et al. 2017, doi:10.7717/peerj.3719 | Local PDF queried | Direct Acropora sub-bleaching heat x apical injury |
| Lock et al. 2022, doi:10.1002/ece3.9345 | Local PDF queried | Fragmentation/regeneration physiology and transcriptomics in Porites |
| Barshis et al. 2013, doi:10.1073/pnas.1210224110 | Local PDF queried | Frontloading and baseline expression in thermal resilience |
| Bay and Palumbi 2015, doi:10.1093/gbe/evv085 | Local PDF queried | Rapid acclimation and transcriptional dampening |
| Dixon et al. 2015, doi:10.1126/science.1261224 | Local PDF queried | Heritable heat tolerance and pre-stress expression |
| Seneca and Palumbi 2009/2010, doi:10.1007/s10126-009-9247-5 | Local PDF queried | Natural bleaching qPCR markers in Acropora millepora |
| Thomas et al. 2019, doi:10.1111/mec.15143 | Local PDF queried | Transcriptomic resilience, recovery, and symbiont shuffling |
| Cleves et al. 2020, doi:10.1073/pnas.1920779117 | Local PDF queried | Functional validation of HSF1 heat tolerance role |
| Granados-Cifuentes et al. 2013, doi:10.1186/1471-2164-14-228 | Local PDF queried | Natural colony-level expression variation |
| Traylor-Knowles 2019, doi:10.7717/peerj.6510 | Local record checked | Heat stress, epithelial integrity, Grainyhead/collagen context |
| DeSalvo et al. 2010, doi:10.3354/meps08372 | Web source checked | Classic Acropora palmata heat/bleaching transcriptomics |
| Dixon et al. 2020, doi:10.1111/mec.15535 | Web source checked | Acropora stress-expression meta-analysis |
| Louis et al. 2017, doi:10.1016/j.cbpc.2016.08.007 | Web source checked | Heat biomarker review and limitations |
| Molinari et al. 2025, doi:10.1111/mec.17753 | Web source checked | Recent cnidarian stress-expression review and cBATS panel |
| Stick et al. 2025, doi:10.1007/s00338-025-02722-w | Web source checked | Acropora sub-bleaching/bleaching transcriptomic recovery |
| Stick et al. 2026, doi:10.1002/ece3.72938 | Web source checked | Acropora heat priming and dampened response |
| Ren et al. 2026, doi:10.1016/j.aqrep.2026.103629 | Web source checked | Recent non-Acropora staged tissue-healing transcriptomics |
| Javid et al. 2025, doi:10.1016/j.marenvres.2025.107102 | Web source checked | Targeted Acropora heat biomarkers |
| Selmoni et al. 2025/2026, Nature Communications | Web source checked | Acropora heat-adaptation genomic windows linked to expression annotations |

## Follow-up before manuscript submission

1. Acquire and index local PDFs for Dixon et al. 2020, Molinari et al. 2025,
   Stick et al. 2025, Stick et al. 2026, Ren et al. 2026, Javid et al. 2025,
   Louis et al. 2017, and DeSalvo et al. 2010 if they are cited in the
   manuscript.
2. Forward-chain from Xu 2023, Han 2025, Bonesso 2017, Barshis 2013, Bay 2015,
   Dixon 2020, and Thomas 2019 using Web of Science, OpenAlex, or Google Scholar.
3. Backward-chain Xu 2023 and Han 2025 for older coral/cnidarian regeneration
   expression studies, including candidate-gene studies in Tubastraea and
   Nematostella as context only.
4. Add a formal inclusion/exclusion table if this becomes a systematic review
   claim. Current language should say "working literature synthesis" rather than
   "all papers in the field."
5. After counts and annotation arrive, map candidate families to Acropora pulchra
   orthologs and keep one file with gene IDs, gene symbols, annotation source,
   and orthology confidence.
