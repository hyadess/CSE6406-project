# Gene-Tree Model Misspecification, Branch Support Reliability, and Weighted Species-Tree Inference

**Course Report**  
**Student:** [Your Name]  
**Student ID:** [Your Student ID]  
**Course:** [Course Code and Name]  
**Instructor:** [Instructor Name]  
**Date:** [Submission Date]

## Abstract

Species-tree estimation is important in phylogenetics because genes can have different evolutionary histories. One major reason is incomplete lineage sorting (ILS), where even true gene trees can disagree with the true species tree. ASTRAL is a popular method for estimating species trees from gene trees. wASTRAL is a weighted version of ASTRAL that uses gene-tree branch support and/or branch-length information.

This project studies whether an incorrect evolutionary model used for gene-tree inference reduces the reliability of branch support and changes the final species-tree result. Species trees and gene trees are simulated under low and high ILS conditions. DNA sequences are generated under a GTR+G4 model and analysed under five gene-tree models: GTR+G4, GTR, HKY+G4, HKY, and JC. First, the project measures gene-tree branch correctness and support calibration. Second, it compares resolved ASTRAL-IV, contracted ASTRAL-IV, support-only wASTRAL, and hybrid wASTRAL using normalized Robinson-Foulds error against the known true species tree.

The full production experiment is still running. Completed smaller subsets are exploratory, so they will not be presented as final confirmatory evidence. This report explains the motivation, background, expectations, and experiment setup. The results and conclusion sections will be completed after the ongoing experiments finish.

**Keywords:** phylogenetics, species tree, gene tree, incomplete lineage sorting, model misspecification, ASTRAL, wASTRAL, branch support

## 1. Introduction

Phylogenetic trees represent evolutionary relationships. A species tree describes the evolutionary history of species, while a gene tree describes the history of one gene or locus. In a simple situation, all gene trees would agree with the species tree. In real and simulated biological data, this often does not happen. Different loci can support different evolutionary histories.

One important biological reason for this disagreement is incomplete lineage sorting (ILS). ILS happens when genetic variation from an ancestral population continues through multiple speciation events. The ancestral gene lineages may coalesce in different orders. Therefore, a true gene tree can have a different topology from the true species tree.

This distinction is very important for the present study. A difference between an estimated gene tree and its corresponding true gene tree is a gene-tree inference error. However, a difference between a true gene tree and the true species tree may be a normal result of ILS. This biological discordance should not be counted as an estimation error.

ASTRAL is a coalescent-based method for estimating species trees from many gene trees [1]. It uses quartet relationships induced by the input gene trees and searches for a species tree with high agreement across these quartets. ASTRAL is useful under ILS because it is designed for the multispecies coalescent setting. However, ordinary ASTRAL does not directly use branch-support values.

wASTRAL is a weighted extension of ASTRAL [2]. It can give more influence to gene-tree relationships with strong branch support and less influence to relationships with weak support. This may improve accuracy when support values correctly show which relationships are reliable.

There is an important possible problem. Gene-tree inference needs an evolutionary model. If the selected model does not represent the true sequence evolution process well, the inferred gene tree can be inaccurate. More seriously, an incorrect branch can receive high support. A support-weighted method may then give extra importance to misleading information.

This project asks the following question:

> Does gene-tree evolutionary-model misspecification reduce the reliability of branch support enough to change the benefit of support-weighted ASTRAL, and is this effect influenced by incomplete lineage sorting?

The project has two connected stages. Stage 1 studies branch correctness and support reliability in inferred gene trees. Stage 2 studies whether these gene-tree properties affect final species-tree accuracy.

## 2. Previous Work and Background

### 2.1 Gene-tree discordance and ILS

Gene-tree discordance means that gene trees do not have the same topology. It can happen because a gene tree is estimated incorrectly, but it can also happen because of genuine biology. ILS is an important biological cause of discordance. Under ILS, some ancestral lineages do not coalesce before later speciation events. This can produce true gene trees that do not match the species tree.

The multispecies coalescent provides a framework for understanding this process. It allows many gene-tree histories to arise within one species tree. High ILS produces more biological disagreement among true gene trees. This can make species-tree estimation difficult, especially when estimated gene trees also contain inference errors.

For this reason, the project creates low-ILS and high-ILS conditions. The treatments are checked using normalized RF distance between true gene trees and the true species tree. Thus, ILS is verified from the simulated histories, instead of being assumed only from the input population-size parameter.

### 2.2 ASTRAL

Mirarab et al. introduced ASTRAL as a method for species-tree estimation from a set of gene trees [1]. ASTRAL uses quartets, which are relationships among four taxa. It searches for a species tree that agrees with a high number of gene-tree quartets. ASTRAL is appropriate for this project because it is designed for settings with ILS.

The quality of the input gene trees can still affect practical performance. A fully resolved gene tree may contain low-confidence branches. Therefore, this project uses two unweighted ASTRAL-IV baselines. The first uses fully resolved inferred gene trees. The second contracts branches with aBayes support below 0.90 before running ASTRAL. Contracting a branch means representing it as unresolved instead of treating a weak relationship as certain.

### 2.3 wASTRAL and gene-tree uncertainty

Zhang and Mirarab proposed weighted ASTRAL methods to improve quartet-based species-tree estimation when input gene trees are uncertain [2]. wASTRAL can use support values, branch lengths, or both. The central idea is that all gene-tree relationships may not deserve equal influence.

This study uses two weighted approaches. Support-only wASTRAL is the main mechanistic analysis because it directly tests whether branch support is useful. Default hybrid wASTRAL is also included because it combines support and branch-length information and represents a practical standard configuration.

wASTRAL is not guaranteed to improve every dataset. Its usefulness depends on whether high support is associated with correct branches. If high support is well calibrated, weighting may help. If high support is assigned to incorrect branches, weighting can lose its advantage or possibly be harmful.

### 2.4 IQ-TREE, support, and model misspecification

IQ-TREE is a maximum-likelihood phylogenetic inference program [3]. In this project, it is used to infer gene trees under selected DNA substitution models. It also produces aBayes branch support values. aBayes is a fast likelihood-based branch-support method [4].

Support is often interpreted as confidence, but it should not automatically be treated as a true probability of correctness. Calibration must be measured. For example, branches with average support near 0.95 should ideally be correct about 95 percent of the time. If they are correct much less often, support is overconfident.

An evolutionary model describes how DNA bases change over time. GTR is flexible because it allows different substitution rates. HKY is simpler because it has fewer parameters. JC is simpler again because it assumes equal base frequencies and equal substitution rates. DNA sites can also evolve at different speeds. Gamma-distributed rate heterogeneity models this variation.

Model misspecification happens when the analysis model is too simple or ignores relevant features of the generating process. This can affect inferred topology, branch length, and support values. The project does not assume that all misspecified models must form a perfect best-to-worst order. Instead, it measures their effect empirically under the same simulated DNA data.

### 2.5 Simulation tools

SimPhy can simulate species trees and gene trees under the multispecies coalescent [5]. AliSim is a fast sequence simulator included with IQ-TREE [6]. Together, these tools make a controlled experiment possible. The true species tree, true gene trees, ILS condition, sequence length, and generating model are known. Therefore, estimated trees can be compared directly with the truth.

## 3. Research Intuition and Expectations

The main logic of the project is:

```text
gene-tree inference model
        -> gene-tree branch correctness
        -> reliability of branch support
        -> usefulness of support weighting in wASTRAL
        -> species-tree accuracy
```

GTR+G4 is the correct-model reference condition. It retains the main features of the sequence-generating process: flexible substitution rates and rate heterogeneity among sites. Under this condition, incorrect branches are expected to have lower support more often than correct branches. Therefore, wASTRAL may use support information successfully.

The misspecified models may cause more gene-tree errors. The more important possibility is that some wrong branches receive high support. In that case, wASTRAL may give additional weight to a relationship that is wrong.

High ILS may make this effect more visible. High ILS creates more genuine disagreement among true gene trees. This does not mean true gene trees are incorrect. It means that the species-tree method must combine more competing gene histories. If model misspecification creates confidently incorrect gene-tree branches, the species-tree estimation task may become harder.

The study tests these hypotheses:

**H1: Model misspecification.** Misspecified gene-tree models will generally have higher branch error than GTR+G4.

**H2: Support reliability.** Under GTR+G4, higher aBayes support will be associated with higher branch correctness. Under misspecification, this relationship may weaken and high-support errors may become more frequent.

**H3: Weighting.** wASTRAL will be more useful when support remains informative. Its advantage over ASTRAL may become smaller when support is poorly calibrated.

**H4: ILS interaction.** The negative effect of misleading support may be stronger under high ILS than under low ILS.

These are expectations, not final results. The output tables must decide whether the hypotheses are supported.

## 4. Experiment Setup

### 4.1 Factorial design

The study uses a factorial simulation design. The main factors are shown below.

| Factor | Levels |
|---|---|
| ILS condition | Low ILS, high ILS |
| Sequence length | 200 bp, 800 bp, 1600 bp |
| Gene-tree analysis model | GTR+G4, GTR, HKY+G4, HKY, JC |
| Species-tree method | Resolved ASTRAL, contracted ASTRAL, support-only wASTRAL, hybrid wASTRAL |

The full planned production design uses 51 taxa, 200 loci per replicate, and 50 independent species-tree replicates. The completed and ongoing subsets are smaller. Some subsets use 21 taxa, 2 replicates, and 7 or 10 loci at selected sequence lengths. Other subsets have larger locus or replicate counts but are still running.

The partial subsets should not be pooled as one completed experiment because they have different numbers of taxa, loci, replicates, and sequence lengths. They are exploratory evidence. The final report should state the exact design of every subset before presenting its results.

### 4.2 Species-tree and gene-tree simulation

SimPhy generates true species trees and true gene trees under the multispecies coalescent. The planned full design contains 50 ingroup species and one outgroup, giving 51 taxa. One individual is sampled from each species.

| Parameter | Planned value |
|---|---:|
| Ingroup taxa | 50 |
| Outgroup taxa | 1 |
| Tree height | 2,500,000 generations |
| Birth rate | $1 \times 10^{-7}$ |
| Mean substitution rate | $1 \times 10^{-7}$ substitutions/site/generation |
| Low-ILS population size | $N_e = 100,000$ |
| High-ILS population size | $N_e = 1,000,000$ |

Population size is the manipulated ILS parameter. A larger effective population size makes ILS more likely because ancestral lineages may remain uncoalesced for longer. The low- and high-ILS conditions use paired species-tree replicate seeds, while the population size differs.

ILS is checked before DNA sequence inference by calculating true-gene-tree to true-species-tree normalized RF distance. The planned acceptance criteria are low ILS pooled nRF less than or equal to 0.25, high ILS pooled nRF greater than or equal to 0.55, and a high-minus-low nRF difference of at least 0.30. For small exploratory subsets, the automatic gate can be skipped because two replicates and few loci cannot reliably certify the production distribution. The ILS summaries are still retained and must be reported honestly.

### 4.3 DNA sequence simulation

AliSim generates DNA alignments along each true gene tree. The generating model is:

$$
\text{GTR}\{1,2,1,1,2\} + F\{0.30,0.20,0.20,0.30\} + G4\{0.20\}.
$$

The GTR component allows different nucleotide changes to have different relative rates. The frequency component fixes unequal base frequencies: 30% A, 20% C, 20% G, and 30% T. The gamma component divides sites into four rate categories and represents variation in evolutionary speed among sites. The shape value 0.20 gives strong rate variation.

The main sequence lengths are 200 bp, 800 bp, and 1600 bp. Shorter alignments usually provide less phylogenetic information. Therefore, they may produce more uncertain gene trees than longer alignments. For a given locus, one simulated alignment is reused under every gene-tree analysis model. This paired design makes model comparisons fair because each model analyses exactly the same sequence information.

### 4.4 Gene-tree inference models

IQ-TREE estimates a maximum-likelihood gene tree for each alignment. Each alignment is analysed under five models.

| Model | Role in this study |
|---|---|
| GTR+G4 | Correct-model reference |
| GTR | Removes gamma rate heterogeneity |
| HKY+G4 | Simpler substitution model with gamma |
| HKY | Simpler substitution model without gamma |
| JC | Severe misspecification: equal rates and equal frequencies |

The GTR+G4 condition is the reference because it matches the important generating-model features. GTR removes the gamma component. HKY+G4 simplifies the substitution model while retaining gamma rate heterogeneity. HKY applies both simplifications. JC is the strongest stress-test condition.

IQ-TREE also calculates aBayes support for internal branches. These support values are used in the Stage 1 calibration analysis and in the weighted species-tree analyses.

### 4.5 Stage 1: branch correctness and support calibration

Stage 1 compares each inferred gene tree with the true gene tree for the same locus. It does not compare inferred gene trees directly with the species tree, because true gene trees can differ from the species tree due to ILS.

For every estimated internal branch, the pipeline records the ILS condition, species-tree replicate, locus, sequence length, analysis model, aBayes support, and whether the branch is present in the corresponding true gene tree.

The main Stage 1 outcomes are described below.

**Branch error rate.** This is the proportion of estimated internal branches that are not found in the true gene tree. Lower is better.

**Support calibration.** Branches are grouped into fixed support bins. For each bin, average support is compared with observed branch correctness. A well-calibrated support value should be close to observed correctness.

**Brier score.** This is the average squared difference between a support value and the binary correctness outcome. Lower values indicate better reliability.

**Expected calibration error (ECE).** This summarizes the gap between mean support and observed correctness across support bins. Lower values indicate better calibration.

**High-support error.** This is the error rate among branches with support at least 0.95 and at least 0.99. This outcome is important because highly supported branches may have larger influence in wASTRAL.

### 4.6 Stage 2: species-tree inference

The inferred gene trees from each condition are used to estimate a species tree with four methods.

| Method | Description |
|---|---|
| Resolved ASTRAL-IV | Uses fully resolved inferred gene trees |
| Contracted ASTRAL-IV | Contracts aBayes branches below 0.90 |
| Support-only wASTRAL | Uses support weighting, mode 2, with `-B` |
| Hybrid wASTRAL | Uses default hybrid weighting with `-B` |

Resolved ASTRAL is the first unweighted baseline. Contracted ASTRAL is the second unweighted baseline, and it removes weak branches before inference. Comparing wASTRAL with both baselines helps separate the benefit of weighting from the benefit of simply removing low-support branches.

Support-only wASTRAL is the main mechanistic method. Hybrid wASTRAL is included because it is a commonly useful practical configuration. Weighted and unweighted analyses use the same estimated gene-tree topologies. When an internal branch has no numeric support label, the weighted input preparation assigns the documented local-Bayesian lower bound of 1/3 and records the event.

### 4.7 Evaluation and statistical unit

Every estimated species tree is compared with the true species tree using normalized Robinson-Foulds (nRF) error [7].

$$
\text{nRF} = 0
$$

means that the estimated topology matches the true species-tree topology. A larger nRF value means greater disagreement with the true species tree.

The main paired comparison is:

$$
\Delta = \text{weighted nRF} - \text{baseline nRF}.
$$

- $\Delta < 0$ means weighting helped.
- $\Delta = 0$ means no topological difference.
- $\Delta > 0$ means weighting hurt.

The main comparisons are support-only wASTRAL minus contracted ASTRAL and hybrid wASTRAL minus contracted ASTRAL. Comparisons with resolved ASTRAL are secondary but still useful.

The species-tree replicate is the independent biological unit. Loci and branches are repeated observations within a replicate. Therefore, final uncertainty analysis should use replicate-level paired comparisons, such as a paired bootstrap or permutation test. Pooled branch results are descriptive and should not be treated as many independent biological replicates.

## 5. Results

**Write this section after the running experiments finish.** Use this order:

1. Show that low and high ILS are separated in `ils_summary.csv`.
2. Show Stage 1 branch error, Brier score, ECE, and high-support error.
3. Show calibration plots using `stage1_calibration.csv`.
4. Compare the four species-tree methods using nRF.
5. Report paired delta values from `stage2_species_tree_error.csv`.
6. Separate completed exploratory subsets from still-running experiments.

## References

[1] S. Mirarab, R. Reaz, M. S. Bayzid, T. Zimmermann, M. S. Swenson, and T. Warnow, "ASTRAL: Genome-scale coalescent-based species tree estimation," *Bioinformatics*, vol. 30, no. 17, pp. i541-i548, 2014, doi: 10.1093/bioinformatics/btu462.

[2] C. Zhang and S. Mirarab, "Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees," *Molecular Biology and Evolution*, vol. 39, no. 12, 2022, doi: 10.1093/molbev/msac215.

[3] L.-T. Nguyen, H. A. Schmidt, A. von Haeseler, and B. Q. Minh, "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies," *Molecular Biology and Evolution*, vol. 32, no. 1, pp. 268-274, 2015, doi: 10.1093/molbev/msu300.

[4] M. Anisimova, M. Gil, J.-F. Dufayard, J. Dessimoz, and O. Gascuel, "Survey of Branch Support Methods Demonstrates Accuracy, Power, and Robustness of Fast Likelihood-based Approximation Schemes," *Systematic Biology*, vol. 60, no. 5, pp. 685-699, 2011, doi: 10.1093/sysbio/syr041.

[5] D. Mallo, L. de O. Martins, and D. Posada, "SimPhy: Phylogenomic Simulation of Gene, Locus, and Species Trees," *Systematic Biology*, vol. 65, no. 2, pp. 334-344, 2016, doi: 10.1093/sysbio/syv082.

[6] N. Ly-Trong, S. Naser-Khdour, R. Lanfear, and B. Q. Minh, "AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era," *Molecular Biology and Evolution*, vol. 39, no. 5, 2022, doi: 10.1093/molbev/msac092.

[7] D. F. Robinson and L. R. Foulds, "Comparison of Phylogenetic Trees," *Mathematical Biosciences*, vol. 53, no. 1-2, pp. 131-147, 1981, doi: 10.1016/0025-5564(81)90043-2.