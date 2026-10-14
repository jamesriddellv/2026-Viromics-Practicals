# 2026-Viromics-Practicals
 Practical examples for the 8th International Viromics Workshop, Oct 14-16 2026. Contains a guided walkthrough of virus identification, annotation, taxonomy, host prediction, and ecological statistics. The following work is based on results from ["Viruses help shape microbiome response to polyphenol rewiring of methane-suppressed peat microcosms"](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003925)


## Day 1 and 2 exercise
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jamesriddellv/2026-Viromics-Practicals/blob/main/day1_and_2_Bioreactor_viruses_deep_dive_identification.ipynb)

## Day 3 exercise
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jamesriddellv/2026-Viromics-Practicals/blob/main/day3_Bioreactor_viruses_deep_dive_host_prediction.ipynb)


 ### Goals for these hands-on exercises:
 1. Get comfortable using exploratory data analysis in python (with our without AI-assistance) to analyze outputs from geNomad, CheckV, DRAMv, iPHoP, vConTACT3, and relative abundance tables.
 2. Learn dimensionality reduction and statistical analyses for sparse ecological data.
 3. Learn how to leverage multiple viromics tool outputs to build a holistic understanding of the ecological role of viruses.
 4. Learn different data visualization techniques for communicating viromics results.

### Non-goals:
1. We will not cover how to run these tools on the HPC (Dr. Bonnie Hurwitz will cover this topic on the final day of the workshop).
2. We will not learn how to code in python (data structures, algorithms).
3. We will not explore the differences between viromics tools that achieve the same goal (i.e., comparing virus ID tools).


### Context: 
 Human activities are accelerating permafrost thaw and subsequent methane emissions from increased microbial activity, prompting microbiome engineering efforts as an emissions mitigation strategy.

It was recently demonstrated that catechin, a polyphenol subunit of condensed tannin commonly found in tea leaves, could drastically reduce methane emissions (>80%) in thawed permafrost microcosms by enriching catechin-degrading prokaryotes that outcompeted methanogens for hydrogen.

However, viral contributions to such microbiome-level responses remain unexplored.

### Your goal is the answer the following questions:

Does catechin alter virus community activity? If so, by what timepoint is a treatment effect observed?
Which host MAGs are these vOTUs predicted to infect? Do the most active vOTUs infect catechin-degrading MAGs?
What kind of viruses are the most active vOTUs in these samples? What genes do these viruses encode? Are they lytic or lysogenic?

Click through the following notebooks, do some additional exploratory data analysis, and answer the reflection questions to build an intuition for how to interpret the outputs of various viromics tools.

Day 1 and 2 exercises cover virus identification, functional annotation, and ecological statistics. Day 3 covers virus-host prediction, virus taxonomy, and combines these with outputs from day 1 and 2 to build a picture of virus community activity within the bioreactor.

Please open a GitHub issue or email me at riddell.26@osu.edu if you have additional questions or notice anything that can be improved.

