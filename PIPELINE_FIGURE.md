# Experimental pipeline and tested options

This figure summarizes the implemented workflow for testing how gene-tree model
misspecification changes branch-support reliability and the downstream
difference between unweighted ASTRAL and support-weighted wASTRAL. Every edge
contains at least one bracketed reference number; the numbers map to the papers
listed immediately below the figure.

```mermaid
flowchart TB
    subgraph S0[Experimental design]
        direction LR
        PROFILE["Run profiles<br/>paper subset: 21 taxa, 10 loci, 2 replicates<br/>production: 51 taxa, 200 loci, 50 replicates"]
        ILS["ILS treatments<br/>low: Ne = 100,000<br/>high: Ne = 1,000,000"]
    end

    subgraph S1[Tree simulation]
        direction LR
        SIMPHY["Conditioned SimPhy MSC histories<br/>50 ingroup + 1 outgroup; one individual/species<br/>birth = 1e-7; ingroup height = 2,500,000<br/>substitution rate exponential mean = 1e-7"]
        TRUE_ST["True species tree"]
        TRUE_GT["True gene tree for every locus"]
    end

    PROFILE -->|"[1]"| SIMPHY
    ILS -->|"[1]"| SIMPHY
    SIMPHY -->|"[1]"| TRUE_ST
    SIMPHY -->|"[1]"| TRUE_GT

    subgraph S2[Sequence simulation]
        direction LR
        LENGTH["Sequence-length options<br/>200, 800, 1,600 bp"]
        GENMODEL["Generating model<br/>GTR rates 1,2,1,1,2<br/>F: A/C/G/T = .30/.20/.20/.30<br/>G4: alpha = .20"]
        ALISIM["IQ-TREE AliSim"]
        ALIGN["One alignment per locus<br/>shared across all analysis models"]
    end

    TRUE_GT -->|"[2]"| ALISIM
    LENGTH -->|"[2]"| ALISIM
    GENMODEL -->|"[2]"| ALISIM
    ALISIM -->|"[2]"| ALIGN

    subgraph S3[Gene-tree inference options]
        direction LR
        MODELS["IQ-TREE analysis models<br/>GTR+G4: correct-model reference<br/>GTR: remove gamma heterogeneity<br/>HKY+G4: simpler substitution model<br/>HKY: simpler model without gamma<br/>JC: severe misspecification"]
        IQTREE["Maximum-likelihood gene-tree inference<br/>aBayes branch support; threads = AUTO"]
        EST_GT["Estimated gene trees<br/>topology + aBayes support"]
    end

    ALIGN -->|"[3]"| IQTREE
    MODELS -->|"[3]"| IQTREE
    IQTREE -->|"[3,4]"| EST_GT

    subgraph S4[Stage 1: mechanism]
        direction LR
        ILS_CHECK["Pre-inference ILS gate<br/>low nRF <= .25; high >= .55; gap >= .30"]
        BRANCH["Branch correctness and support calibration<br/>estimated split vs corresponding true-gene split"]
        S1OUT["Stage 1 outputs<br/>branch error; Brier score; ECE; support bins<br/>P(wrong | support >= 0.95/.99)<br/>missing-support counts"]
    end

    TRUE_GT -->|"[7]"| ILS_CHECK
    TRUE_ST -->|"[7]"| ILS_CHECK
    TRUE_GT -->|"[7]"| BRANCH
    EST_GT -->|"[4,7]"| BRANCH
    ILS_CHECK -->|"[7]"| S1OUT
    BRANCH -->|"[4,7]"| S1OUT

    subgraph S5[Stage 2: consequence]
        direction LR
        ASTRAL["ASTRAL-IV baselines<br/>resolved input<br/>aBayes < .90 contracted input"]
        PREP["wASTRAL support preparation<br/>preserve topology; missing internal support = 1/3"]
        WASTRAL["wASTRAL<br/>support-only mode 2; -B<br/>default hybrid; -B"]
        UST["Resolved + contracted<br/>unweighted species trees"]
        WST["Support-only + hybrid<br/>weighted species trees"]
        EVAL["Species-tree evaluation<br/>normalized RF against truth<br/>paired deltas vs both baselines"]
    end

    EST_GT -->|"[5]"| ASTRAL
    EST_GT -->|"[6]"| PREP
    PREP -->|"[6]"| WASTRAL
    ASTRAL -->|"[5]"| UST
    WASTRAL -->|"[6]"| WST
    TRUE_ST -->|"[7]"| EVAL
    UST -->|"[7]"| EVAL
    WST -->|"[7]"| EVAL

    S1OUT -->|"[4,6]"| SYNTHESIS["Cross-stage interpretation<br/>Does model misspecification make support less reliable<br/>and change the benefit of support weighting?"]
    EVAL -->|"[6,7]"| SYNTHESIS
```

## Reference papers

1. Mallo, D., Martins, L. de O., & Posada, D. (2016). **SimPhy:
   Phylogenomic Simulation of Gene, Locus, and Species Trees.** *Systematic
   Biology, 65*(2), 334–344. https://doi.org/10.1093/sysbio/syv082

2. Ly-Trong, N., Naser-Khdour, S., Lanfear, R., & Minh, B. Q. (2022).
   **AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the
   Genomic Era.** *Molecular Biology and Evolution, 39*(5), msac092.
   https://doi.org/10.1093/molbev/msac092

3. Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015).
   **IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
   Maximum-Likelihood Phylogenies.** *Molecular Biology and Evolution, 32*(1),
   268–274. https://doi.org/10.1093/molbev/msu300

4. Anisimova, M., Gil, M., Dufayard, J.-F., Dessimoz, C., & Gascuel, O.
   (2011). **Survey of Branch Support Methods Demonstrates Accuracy, Power,
   and Robustness of Fast Likelihood-based Approximation Schemes.**
   *Systematic Biology, 60*(5), 685–699.
   https://doi.org/10.1093/sysbio/syr041

5. Mirarab, S., Reaz, R., Bayzid, M. S., Zimmermann, T., Swenson, M. S., &
   Warnow, T. (2014). **ASTRAL: Genome-scale Coalescent-based Species Tree
   Estimation.** *Bioinformatics, 30*(17), i541–i548.
   https://doi.org/10.1093/bioinformatics/btu462

6. Zhang, C., & Mirarab, S. (2022). **Weighting by Gene Tree Uncertainty
   Improves Accuracy of Quartet-based Species Trees.** *Molecular Biology and
   Evolution, 39*(12), msac215. https://doi.org/10.1093/molbev/msac215

7. Robinson, D. F., & Foulds, L. R. (1981). **Comparison of Phylogenetic
   Trees.** *Mathematical Biosciences, 53*(1–2), 131–147.
   https://doi.org/10.1016/0025-5564(81)90043-2

## Implemented comparison grid

The same alignment for a locus is reused across all five IQ-TREE analysis
models. Each resulting estimated-gene-tree set is sent to resolved and
0.90-aBayes-contracted ASTRAL-IV, support-only wASTRAL, and default hybrid
wASTRAL. This paired design isolates the effect of support weighting within
each combination of ILS level, sequence length, replicate, and analysis model.

Each Stage 2 weighted-minus-baseline delta is interpreted as follows:

- `delta < 0`: support weighting helped.
- `delta = 0`: support weighting made no topological difference.
- `delta > 0`: support weighting hurt.

The paper-subset profile shown in the figure is the reduced grid exercised in
`work/paper_subset_10x2`. The production profile records the full defaults in
the pipeline; it should not be described as completed until all expected
outputs have been generated and checked.
