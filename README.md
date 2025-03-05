# LittleBites
an R package to clean FACS enriched bulk RNA-seq datasets

Install with devtools::install_github('https://github.com/alecbarrett/LittleBites')


## this package is designed to aid in studies that intend to examine the expression profiles for a pure cell population, but in practice find that the enrichment methods are not sufficient to eliminate all contaminants

LittleBites is an iterative subtraction method that uses 3 inputs:
1) A bulk RNA-seq counts profile (usually from a gene x sample matrix)
2) A single cell profile for the target cell type (which is expected in the bulk data), and putative contaminant cell types (which should be removed) usually from a gene x cell_types matrix
3) A ground truth table of known gene expression across the cell types in the bulk and single cell datasets -- This should be present in the form of a genes x cell_types table, with binary values (1 = gene is expressed in the cell type, 0 = gene is not expressed in the cell type)

## LittleBites runs an algorithm for conservative iterative subtraction

Each loop in the algorithm has 4 steps
1) Model the sample as a mixture of the single cell reference, extracting the proportion of the sample estimated to be from each contaminant
2) Using a series of learning rates, subtract out the single cell profiles, weighted by the proportions estimates, and a gene level specificity weight (coming from the presumption that the rarer the gene's expression, the more acceptable its subtraction is)
3) For all learning rates tested, calculate an AUROC value for the ground truth genes, and select the learning rate that maximizes the gain in AUROC (in the case of a tie, pick the lowest learning rate possible to subtract less).
4) If no learning rates are better than not subtracting at all, STOP. Else, go back to step 1

<img width="666" alt="subtraction loop" src="https://github.com/user-attachments/assets/6a7a3ca9-3438-4fb2-b50f-2d7e4d6a9b9b" />

## General recommendations

1) For ground truth genes, use generic genes, not highly specific genes. This process makes few assumptions about the genes that are being put into it, and so it ends up working best if the genes used to guide the subtraction are as broad as possible to avoid overfitting.
2) Use multiple sets of testing genes, and here you *should* include highly specific genes, make sure you're not inadvertently removing cell type markers because of an issue with your reference set for instance

As an example, in [this paper currently in review](https://www.biorxiv.org/content/10.1101/2025.01.26.634951v1) we worked with FACS sorted neuronal samples, and subtracted out non-neuronal profiles (intestine, muscle, etc...) and we used ubiquitous and non-neuronal genes as the training ground truth. All this did was ensure the algorithm's subtraction stopped when it achieved a good balance between removing genes from entirely different tissues, and preserving genes expressed in every cell type. We reserved a subset of those genes for testing, and validated the approach using a totally different set of genes that are only expressed in certain neurons, so that we could be sure that the little bites being taken out were not having unforeseen side effects.

<img width="1684" alt="AUCPR values for various deconvolution methods using neuron and generic markers" src="https://github.com/user-attachments/assets/0d6f347c-2b06-4f82-a27c-cabd2f1a8633" />


