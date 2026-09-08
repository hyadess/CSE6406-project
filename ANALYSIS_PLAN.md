# Confirmatory analysis plan

This plan defines the independent unit, primary contrasts, and uncertainty
analysis before the production results are interpreted. It follows the paired
simulation structure and avoids treating loci or branches as independent
biological replicates.

## Independent unit and blocking

The independent unit is the species-tree replicate (`replicate`, 1-50). The
same replicate seed produces the same conditioned species-tree topology and
branch times for the low- and high-ILS arms; only haploid effective population
size changes. Within each replicate, sequence lengths, inference models, and
species-tree methods are repeated measurements.

Loci and internal branches provide information within a replicate, but they do
not increase the number of independent species-tree histories. The 1,500 rows
in `stage2_species_tree_error.csv` are therefore 50 paired blocks, not 1,500
independent observations.

## Preregistered ILS gate

Before sequence inference, pooled true-gene/species normalized RF must satisfy:

- low ILS <= 0.25;
- high ILS >= 0.55;
- high minus low >= 0.30.

Zhang and Mirarab report approximately 10% and 70% discordance for the S200
low- and high-ILS categories. The wider bounds here allow the intentionally
different 51-taxon design while preventing two merely ordered but biologically
similar treatments from passing. These bounds were fixed before the revised
production run.

Failure stops the pipeline after writing `ils_verification.csv` and
`ils_summary.csv`. Parameters must then be revised and a new output directory
used; the failed gate must not be tuned after inspecting Stage 1 or Stage 2.

## Stage 1 outcomes

Primary Stage 1 outcomes are replicate-level branch error, Brier score, and
expected calibration error. Calibration-bin tables, mean support within bins,
and conditional error among branches with support >=0.95 and >=0.99 are
descriptive mechanism summaries.

Comparisons use replicate-level summaries. Pooled Wilson intervals over
branches are descriptive only because branches and loci within a replicate are
dependent. Confirmatory uncertainty is obtained by resampling whole replicates
or by a repeated-measures model with replicate as a blocking/random effect.

## Stage 2 outcomes and contrasts

The primary outcome is normalized RF error against the true species tree. The
four methods are:

1. ASTRAL-IV using fully resolved gene trees;
2. ASTRAL-IV after contracting aBayes support below 0.90;
3. support-only wASTRAL (`--mode 2 -B`);
4. default hybrid wASTRAL (`-B`).

The two primary paired contrasts are:

- support-only wASTRAL minus contracted ASTRAL;
- hybrid wASTRAL minus contracted ASTRAL.

Negative differences favor wASTRAL. Comparisons with resolved ASTRAL are
secondary mechanism contrasts. Effects are reported for every ILS x sequence
length x analysis-model cell with replicate-level paired confidence intervals.

A repeated-measures mixed model may use nRF as the response, fixed effects for
ILS, sequence length, analysis model, species-tree method, and preregistered
interactions, with replicate as a random intercept/block. Because nRF is
bounded and quantized, conclusions must also be checked using paired
replicate-block bootstrap intervals and permutation tests on paired deltas.

## Multiplicity and effect reporting

The two primary method contrasts are tested within the prespecified factorial
model. Follow-up cellwise tests use Holm correction within each family of
analysis-model comparisons. Report effect sizes and confidence intervals even
when adjusted tests are not significant; do not rank models solely by p-value.

Normalized quartet distance to the true species tree is a desirable secondary
topological outcome because it has finer resolution than RF. It must be
computed as a common truth-distance metric. Raw ASTRAL and wASTRAL optimization
scores must not be subtracted because weighted and unweighted objectives are
on different scales.

## Scope

GTR+G4 is the correctly specified model-family reference. GTR, HKY+G4, HKY,
and JC are the four misspecified analysis conditions. Support-only wASTRAL
isolates the influence of branch support; hybrid wASTRAL represents the
default practical configuration. Neither result should be described as an
exact replication of the published S100 or S200 datasets.
