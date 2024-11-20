# nepobabies
Code and requisite data to replicate analyses in Orkney et al., 2025
This repository describes a manuscript in-prep, and the individual codes will be updated with full run-time estimates prior to publication

Code repository for Orkney, Rothier and Hedrick 2025 "Great expectations: How parental investment in early development supercharges bird evolution.",
 run in R version 4.2.3 
This repository contains commented codes to enable readers to reproduce the analyses conducted in Orkney, Rothier and Hedrick, 2025.
These scripts were run on a machine with the following specifications: 
Machine specs: 12th Gen Intel(R) Core(TM) i9-12900K 3.20 GHz 128 GB (128 GB usable) 64-bit operating system, x64-based processor
However, no non-standard hardware is required to reproduce the analyses here.
Installation guide: See https://www.r-project.org/ for guidance installing and using R,
 including expected installation times and the installation of dependant packages.


Data preparation: Bird_data_preparation_08_30_2024.R (Run time under 5 minutes) 
This script loads and prepares datasets for subsequent analysis. 
Bird landmark data should be sourced from: https://doi.org/10.18563/journal.m3.125 
A family tree describing bird relationships to each other should be sourced from: 
https://doi.org/10.5281/zenodo.28343
 The script will require the following additional files:
 Bird_landmark_info_23Nov2019.csv ; Eco_meta_data.csv

Bird_landmark_info_23Nov2019.csv This is a metadata file containing information required 
to extract data from the Bjarnason Landmark deposition.

Eco_meta_data.csv This is a metadata file containing bird genus names and masses required
 for data preparation. (it also contains a variety of other information but it is not needed here)

Bird_data_preparation_08_30_2024.R will produce:

Csize.22.10.2022.RData

This is an array of centroid sizes of landmark constellations,
approximating the volumes of several bones in the avian skeleton for the bird species considered here.
The source landmark constellations were originally compiled by Bjarnason et al., 2021
(https://doi.org/10.18563/journal.m3.125)
This file will be generated during data preparation and is not supplied in this GitHub repository.

coords.22.10.2022.RData

This object refers to original landmark constellations prepared from Bjarnason et al., 2021.
(https://doi.org/10.18563/journal.m3.125) This file will be generated during data 
preparation and is not supplied in this GitHub repository.

tree.22.10.2022.RData

This is a phylogenetic tree pruned from Prum et al., 2015 (https://doi.org/10.1038/nature15697),
representing the evolutionary relationships between the bird species considered here. 
This file will be generated during data preparation and is not supplied in this GitHub repository. 
The tree file 'Avian-TimeTree.tre' may be sourced from: 
https://github.com/carolinafishes/Prum_et_al_2015/blob/v.1.0.0/Trees/BEAST/Avian-TimeTree.tre

tree.names.22.10.2022.RData

This is a vector of species names, which matches the birds in this study to congeners on Prum et al.'s phylogeny.

masses.22.10.2022.RData

This is a vector of representative body masses for the avian species considered here, 
aggregated from the Cornell handbook of the birds of the world. (https://birdsoftheworld.org/bow/home)

Brinkworth_builder_11_14_2024.R
This script interacts with the supplementary data file available in Brinkworth et al., 2023, 
in order to compile a specimen-wise matrix of complete wing and leg skeletal proportions across
a wide sample of altricial and precocial birds.
The script depends on Supplementary_Data_1.csv and Prum_merge_hackett_stage2_1K_mcc.tree

Supplementary_Data_1.csv
This file is sourced from the Brinkworth et al., 2023 supplement: https://doi.org/10.1038/s41467-023-41415-2

Prum_merge_hackett_stage2_1K_mcc.tree
This file is sourced from the Brinkworth et al., 2023 supplement: https://doi.org/10.1038/s41467-023-41415-2

The script Brinkworth_builder_11_14_2024.R will produce:

altriciality_scores_brinkworth_11_14_2024.csv
This is a spreadsheet of specimen-specific species-wise wing and leg skeletal proportions and associated developmental
mode classification attributed in line with Starck 1993: https://doi.org/10.1007/978-1-4615-9582-3_6

Other supplied data:

precocial_and_altricial_mass_matched_pairs_10_24_2023.csv
This is a spreadsheet containing carefully selected precocial and altricial bird species from the dataset of Bjarnason et al., 2021.
The altricial and precocial species have a similar body mass distribution, which means that a tendency for precocial birds to have 
larger body masses can be dismissed as a confounding factor in analyses. 
The selected altricial and precocial birds are separated by a similar overall range of evolutionary divergence times, 
and both assemblages are polyphyletic. 

matched_birds_brinkworth_06_03_2024.csv
This is a spreadsheet containing carefully selected precocial and altricial bird species from the dataset of Brinkworth et al., 2023.
As before, body mass has been excluded as a confounding factor by selecting altricial and precocial specimens whose species means
are within 10% of one another. 

flight_masses_22_10_2022_plus_A.csv
species-wise presence-absence of flight-style categorizations sourced from Orkney et al., 2021 (https://doi.org/10.1038/s41559-021-01509-w)

standardised_foot_scores_11_07_2023.csv
species-wise presence-absence of foot-use categroizations, simplified from Brigit Tronrud's contribution to Orkney et al., 2021 
(https://doi.org/10.1038/s41559-021-01509-w)


Once all required data is downloaded and both preparation scripts have been successfully executed it will be possible to replicate
analyses: 
Original illustrations and silhouettes presented in our manuscript are not stored on this GitHub repository, 
and lines that call .jpg or .png files can therefore be disabled. 


Battenberg_06_03_2024.R
This script computes the strength evolutionary correlations between the allometric residuals of limb bone relative sizes across
mass-matched subsets of precocial and altricial birds from the Bjarnason and Brinkworth datasets, expressed as a sample-size 
adjusted Z-score computed under a 2-Block phylogenetic partial least squares approach. 
The analysis demonstrates that forelimb and hindlimb evolution is highly correlated in assemblages of precocial bird species, 
but that forelimb and hindlimb evolution proceed independently of one another in assemblages of altricial bird species. 
Given that precocial developmental strategies are commonly inferred to represent the ancestral avian condition, this suggests
that divergent forelimb and hindlimb evolution is a derived pattern of evolution within the avian crown group.  

Eco_style_diversity_precocial_altricial_w_phylograms_10_30_2024.R
This script seeks to determine whether assemblages of altricial bird species explore a wider repertoire of locomotory ecologies
than do assemblages of precocial bird species. We might expect this to be true if divergent wing and leg evolution in altricial birds
permits the independent adaptation of these structures to distinct new ecologies.
The script specifically engages flight-style and foot-use variety scores that accompany the dataset of Orkney et al., 2021 (https://doi.org/10.1038/s41559-021-01509-w) to achieve this.
The script constructs stacked bar-charts across avian phylogeny, visualizing flight-style and foot-use ecological variety. 
The script produces ecological category richness accumulation curves as a function of sampling effort, subset by developmental strategy and 
locomotor ecology type (flight-style/foot-use). 
A bespoke statistical test is included to compute a 1-sided p-value for the likelihood that the 2-sigma confidence overlaps of
precocial and altricial richness accumulation curves are less than expected under a scenario whereby ecological categories with an 
identical frequency abundance distribution to those observed evolved across the observed phylogeny by a Brownian Walk with no
regard for developmental strategy or autocorrelation-structure among ecological activities. 
The results of these analyses imply that the assemblage of altricial birds has accumulated a larger variety of flight-style and foot-use
ecologies than the precocial assemblage and that this difference is likely to be so great that it cannot be explained by chance.

Birdnest_plots_08_13_2024.R
This script produces ancestral state reconstructions of the relative size of wing and leg skeletons sourced from the dataset of 
Bjarnason et al., 2021. All bird species that can be classified into altricial or precocial developmental strategies are employed 
and body mass is not specifically isolated as a confounding factor; we proceed from the perspective that earlier analyses
isolating this factor have predicted that forelimb and hindlimb evolution should be more divergent in altricial bird lineages. 
Ancestral state reconstruction is computed under a 2-rate Brownian walk that allows altricial and precocial assemblages to evolve
with different variances/rates of change. The exploration of relative wing to leg skeletal size is plotted on the x-axis over time, 
and lineages are manually colorized by their likely developmental strategy. Altricial lineages practicing novel flight-styles and 
foot-uses are identified and colorized distinctly. The distribution of novel flight-styles and foot-uses over the x-axis is illustrated
with a stacked histogram over the time=0 horizon. No attempt is made to reconstruct ecological categorization ancestries. 
The plots which are produced demonstrate that precocial exploration of the x-axis is relatively conservative compared to the 
exploratory evolution of altricial bird lineages. Altricial bird lineages which practice lead to species practicing novel flight-styles
are associated with especially wide-ranging and rapid evolutionary change, but altricial bird lineages which practice novel foot-uses 
tend to exhibit conservative evolution. 
This figure therefore presents qualitative evidence that the divergent evolution of forelimb and hindlimb size that is facilitated by an
altricial developmental strategy has mediated altricial bird species' access to a wide range of flight-styles. No strong comment
can be made about whether adaptation to novel foot-uses might also drive divergent evolution of the wing and leg skeletons in altricial lineages, 
because the absense of a perceived pattern might reflect the semantic categorization of foot-uses across the studied birds. 

Ternary_10_25_2024_updated.R
This script aggregates the grand sum of available bird species from both Brinkworth et al., 2023 and 
Bjarnason et al., 2021. Original landmarks from the Bjarnason dataset are interrogated and long-bone length measurements are extracted.
Altricial and precocial species are then extracted from the pooled data and the script 
produces ternary diagrams representing wing and leg skeletal proprtions and density field estimates for wing and leg proportions subset by developmental strategy.
This is a visualisation, it includes bird species that were excluded from previous analyses, and it mixes the Bjarnason and Brinkworth datasets.
The figure is intended to prompt readers to consider how divergent wing and leg evolution may have fostered the wider exploration of skeletal proportions within-limbs
of altricial bird lineages. It is not a rigorous statistical demonstration of this effect.  

Evolutionary_models_depending_on_developmental_strategy_11_07_2024.R
This script conducts analyses to determine whether divergent evolution of the relative size of the wing and leg skeletons in the bird species
from the Bjarnason et al., 2021 dataset occurs at a significantly different rate depending on altricial or precocial developmental strategy.
A Markov-Chain-Monte-Carlo approach is applied to simulate possible evolutionary histories of binary developmental strategy (altricial/precocial)
across avian phylogeny, with the assumption that precocial development is the ancestral condition and that the evolution of both states occurs at 
equal rates.  
Thereafter, single-rate Brownian walk and multiple-rate Brownian walk models of the relative size of the wing
to the leg skeleton, contingent on inferred developmental strategy across the tree, are computed and compared. 
The results demonstrate that a 2-rate model is overwhelmingly preferred, suggesting that under the assumptions made to conduct the analyses, 
the evolutionary rate of divergent wing and leg evolution is much faster in altricial lineages. We ascribe this result to our earlier
inference that the evolution of individual wing and leg skeletal proportions is divergent in altricial lineages, but synergistic in precocial lineages. 


eco_evo_rates_model_choices_07_11_2024.R
This script computes evolutionary models to determine whether the evolution of relative size of wing and leg in altricial bird lineages
is more rapid in assemblages that practice novel locomotory ecologies, compared to those which do not. 
Novel flight-styles and novel foot-uses are considered independently. 
'Novelty' is assumed when a trait is osberved in altricial species but not observed in precocial species. 
The evolution of individual novel traits across the tree is unknown, and many possible scenarios for their evolution exist.
Instead of trying to accommodate this, we will instead assume that the tendency to acquire novel traits is inherent to lineages in which 
any novel trait originates. 
The tendency for altricial lineages to give rise to novelty can then be mapped as a binary trait across the phylogeny.
Multiple possible stochastic models, such as equal or differenr rates of this process are compared and the favoured model is selected by AIC. 
The apriori assumption of the ancestral condition being either 'novelty' or 'not novelty' is then imposed upon the favoured stochastic model, 
and the assumption that produces the lowest AIC is selected. 
I visually assessed the stability of multiple possible reconstructed histories. 
Thereafter, novel flight-style was assumed to be ancestral to altricial lineages, but novel foot-use was assumed to be derived. 
Propensity to product novelty was stochastically mapped under the above dynamics across phylogeny many times. 
Thereafter evolutionary models of relative wing to leg skeleton size were fitted in the altricial bird assemblage of Bjarnason et al., 2021, 
with a single-rate Brownian and a multiple-rate Brownian depending on lineage propensity to produce novelty. 
We show that propensity to produce novel flight-styles is associated with substantial preference for a higher rate of divergent wing and leg size evolution
in altricial birds. 
We show that propensity to produce novel foot-uses is associated with a subtantial preference for a lower rate of divergent wing and leg size evolution in
alticial birds. 
We interpret the former finding as suggesting that divergent wing and leg evolution is likely actively exploited within altricial lineages during 
adaptation to new flight-styles. 
We interpret the latter findign as either suggesting that hindlimb adaptation does not intersect with divergent wing and leg size evolution in altricial birds, 
or potentially that the semantic categorizations of hindlimb-use are too limited to foster this insight. 
We remain open to the possibility that subtle changes in shape or intramural changes in hindlimb proportion that are permitted by divergent evolution of 
the wing and leg may well facilitate ecological adaptation in altricial birds. 
