# Coral temperature and regeneration gene-expression literature for LTH

Updated: 2026-09-02

Purpose: give the RNA-seq analyst a compact, documented evidence base for
interpreting LTH heat x wound expression results. This is not a target panel and
not a substitute for genome-wide testing.

Core rule:

**Run genome-wide expression analyses first. Use these genes and pathways to
interpret contrasts, score pre-declared modules, and explain mechanisms after the
main effects are estimated.**

## Bottom line

1. The strongest regeneration evidence points to an early wound program
   involving c-Fos/AP-1, JNK/MAPK, Wnt, FGF, metalloproteases, ECM remodeling,
   and later calicoblast/biomineralization programs.
2. The strongest temperature evidence points to proteostasis, redox balance,
   ubiquitin-proteasome activity, immune/apoptosis signaling, ion transport,
   cell adhesion, metabolism, symbiosis, and pigments.
3. The most LTH-like published comparison is heat plus injury. Bonesso et al.
   found that sub-bleaching heat impaired Acropora apical recovery and skeletal
   regeneration, while van de Water et al. found that some early immune responses
   after injury remained largely intact under mild heat.
4. HSP genes are useful but not sufficient. Heat responses depend on intensity,
   duration, acclimation history, sample timing, colony identity, and whether the
   coral is in an acute, dampened, frontloaded, or recovery state.
5. The right LTH test is process-level: does 31 C preserve early tissue closure
   while suppressing the later skeletal/calicoblast regeneration program,
   especially by Day 10?

## Evidence tiers

| Tier | Meaning | Use in LTH |
|---|---|---|
| A | Direct coral expression evidence from Acropora or an LTH-like heat/injury design | Pre-declare module scores and interpret primary contrasts |
| B | Direct coral expression evidence from another genus, life stage, or related stress design | Use as supporting interpretation |
| C | Review, meta-analysis, genomic association, or phenotype-only support | Use as guardrails and hypothesis context |
| D | Lab-only, contradicted, or not yet source-verified | Do not treat as evidence; keep only as a check-list item |

Machine-readable evidence table: `docs/rnaseq/candidate_genes_reference.csv`.
Search provenance: `docs/rnaseq/literature_search_log.md`.

## Process map for LTH

| LTH question | Candidate processes | Most useful contrast |
|---|---|---|
| Did the sampled margin capture a wound signal? | c-Fos/AP-1, JNK/MAPK, TLR/NOD/NF-kB, ECM remodeling, oxidative response | Wounded vs unwounded at Day 1 and Day 3 |
| Does heat spare closure but block regeneration? | Wnt, FGF, calicoblast fate, galaxin/SAARP/mucin, carbonic anhydrase, Ca2+ handling, skeletal matrix proteins | 31 C vs 28 C in wounded Day 10 margins |
| Is the heat response chronic rather than an acute spike? | HSPs, ubiquitin-proteasome, redox, mitochondria, metabolism, cell adhesion, immune/apoptosis | 31 C vs 28 C in unwounded margins across Day 1, 3, 10 |
| Does source C differ from A/D? | Frontloaded/dampened stress modules, redox/proteostasis, symbiosis, cell adhesion, ion transport | temperature x source-thicket label, with SNP clusters added when available |
| Do organismal traits map to expression? | WGCNA modules for stress, symbiosis, pigment, calcification, growth, regeneration | Module eigengenes vs PAM, color, symbionts, growth, morphology milestones |

## Regeneration evidence

Direct Acropora regeneration work is still sparse. Xu et al. (2023) is the key
bulk RNA-seq study for wound healing and regeneration in Acropora millepora. It
shows that the earliest wound stage has the strongest expression shift and
highlights c-Fos, JNK, Wnt, FGF, sprouty, ADAMTS18-like metalloproteases, and
galaxin. This paper is especially useful for the LTH Day 1/Day 3 wound positive
control, but its first strong sampling point is 6 h post-injury, earlier than our
24 h Day 1 samples.

Han et al. (2025) adds cell-type resolution in Acropora muricata. Its useful
contribution is not only a gene list, but the cell-state idea: wounded tissue
involves epidermal cell changes and calicoblast differentiation, with
biomineralization genes such as carbonic anhydrase 2, galaxin/galaxin2, SAARP,
mucin, calumenin, calmodulin, collagens, and skeletal matrix proteins. This is
the best guide for asking whether heat blocks the later skeletal/calicoblast
phase while leaving earlier closure intact.

The van de Water et al. Acropora aspera injury papers supply the immune and
heat-injury bridge. Physical injury triggered phased TLR/NOD, complement/lectin,
phenoloxidase, and GFP-like/chromoprotein responses over 10 days. Under mild
heat, many 24 h injury immune responses were still present, but NF-kB and TRAF6
were reduced and complement/lectin effectors shifted. This supports a precise
LTH prediction: heat may not erase the early wound signal, but it could change
immune allocation and later skeletal repair.

Lock et al. (2022) is Porites, not Acropora, but it is mechanistically valuable
because microfragmentation induced calcium homeostasis disruption, ER stress,
HSPs, redox enzymes, ubiquitin-proteasome activity, mitochondrial energy
production, and calcium-binding proteins. The overlap with heat-stress pathways
is useful for interpreting shared injury/stress modules.

Ren et al. (2026; Turbinaria peltata) is very recent and not yet in the local
library. It independently reports staged tissue healing: early coenosarc healing
with oxidative stress and MAPK signaling, a middle tissue-remodeling/cytoskeleton
phase, and later immune/polyp formation. Use it as supporting evidence, not as an
Acropora anchor.

## Temperature evidence

The Acropora heat literature makes one point clearly: candidate genes are
context-dependent. Barshis et al. (2013) supports frontloading: more tolerant
corals had higher constitutive expression of stress-related genes and a smaller
acute induction. Bay and Palumbi (2015) supports rapid acclimation and
transcriptional dampening: after 7-11 days at 31 C, corals showed increased heat
tolerance with muted expression responses to acute heat, not broad baseline
changes. Stick et al. (2025, 2026) extends this idea in Acropora by showing fast
transcriptome recovery from sub-bleaching stress and dampened response after heat
priming.

Dixon et al. (2020) is the key guardrail for LTH. Across Acropora stress
experiments, high-intensity stress produced a stereotyped environmental stress
response, while lower-intensity stress often produced a different and sometimes
opposite response. Since LTH is chronic and sublethal at 31 C, absence of a
classic acute HSP spike would not mean absence of heat stress.

The gene families that repeatedly matter across heat papers are HSP70/HSP90 and
small HSPs, HSF1, antioxidants and oxidoreductases, ubiquitin/proteasome genes,
TNF/TRAF/apoptosis genes, lectins and immune genes, ion transporters, cell
adhesion/ECM/cytoskeleton genes, mitochondrial and metabolic genes, and
pigment/symbiosis genes. The best LTH implementation is to score these as
pathways/modules and let the genome-wide model decide which genes actually move.

## Analysis implications

- Treat `source_thicket` A/C/D as a phenotype grouping until the SNP/kinship run
  resolves genetic identity. Say "source C" rather than "genet C" when the claim
  depends on identity.
- Do not require HSP70 or HSP90 to move for the heat treatment to be real. At Day
  1/3/10 after chronic exposure, the heat signal could appear as dampening,
  altered recovery, redox/proteostasis load, metabolism, or symbiosis shifts.
- Score regeneration in at least two stages: early wound/closure and later
  skeletal/calicoblast regeneration. Pooling them into one "wound response" can
  hide the main LTH phenotype.
- Keep candidate-gene checks secondary. The primary evidence should be the
  contrast estimates, module effects, enrichment results, and phenotype-linked
  eigengenes.
- Establish orthology to the Acropora pulchra Conn 2025 annotation before naming
  locus-level candidates in a manuscript.

## High-priority source list

| Citation | Why it matters for LTH | DOI or URL |
|---|---|---|
| Xu et al. 2023 | Direct Acropora wound/regeneration RNA-seq; c-Fos, JNK, Wnt, FGF, ADAMTS, galaxin | https://doi.org/10.3389/fevo.2022.979278 |
| Han et al. 2025 | Direct Acropora regeneration scRNA-seq/bulk RNA-seq; calicoblast fate and biomineralization | https://doi.org/10.1038/s42003-025-08089-6 |
| van de Water et al. 2015 Mol Ecol | Acropora injury immune time course over 10 days | https://doi.org/10.1111/mec.13257 |
| van de Water et al. 2015 Hydrobiologia | Acropora heat x injury immune response | https://doi.org/10.1007/s10750-015-2243-z |
| Bonesso et al. 2017 | Sub-bleaching heat impairs Acropora apical recovery and skeletal regeneration after injury | https://doi.org/10.7717/peerj.3719 |
| Lock et al. 2022 | Fragmentation response links calcium homeostasis, ER stress, redox, proteostasis, and growth | https://doi.org/10.1002/ece3.9345 |
| Barshis et al. 2013 | Frontloading model for thermal resilience | https://doi.org/10.1073/pnas.1210224110 |
| Bay and Palumbi 2015 | Rapid acclimation and dampened acute heat response after short exposure to 31 C | https://doi.org/10.1093/gbe/evv085 |
| Dixon et al. 2015 | Heritable heat tolerance linked to expression before stress | https://doi.org/10.1126/science.1261224 |
| Dixon et al. 2020 | Acropora environmental stress response depends on stress intensity | https://doi.org/10.1111/mec.15535 |
| Thomas et al. 2019 | Transcriptomic recovery and symbiont shuffling predict recurrent bleaching vulnerability | https://doi.org/10.1111/mec.15143 |
| Cleves et al. 2020 | HSF1 functional validation by CRISPR in Acropora millepora | https://doi.org/10.1073/pnas.1920779117 |
| Seneca and Palumbi 2009/2010 | Natural bleaching expression in Acropora millepora; catalase, chromoprotein, C-type lectin | https://doi.org/10.1007/s10126-009-9247-5 |
| Granados-Cifuentes et al. 2013 | Natural expression variation in Acropora millepora; caution on colony-level variability | https://doi.org/10.1186/1471-2164-14-228 |
| Louis et al. 2017 | Heat biomarker review; no universal biomarker | https://doi.org/10.1016/j.cbpc.2016.08.007 |
| Molinari et al. 2025 | Recent cnidarian stress-expression review; cBATS candidate panel | https://doi.org/10.1111/mec.17753 |
| Stick et al. 2025 | Acropora transcriptomic recovery after sub-bleaching vs bleaching heat | https://doi.org/10.1007/s00338-025-02722-w |
| Stick et al. 2026 | Acropora heat priming dampens expression response | https://doi.org/10.1002/ece3.72938 |
| Ren et al. 2026 | Non-Acropora staged tissue-healing transcriptomics | https://doi.org/10.1016/j.aqrep.2026.103629 |

