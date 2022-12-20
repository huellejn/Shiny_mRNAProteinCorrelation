# Shiny app for visualization of the correlation between mRNA and protein expression

The app requires a .RData object as input, that contains  
- coldata (metadata on the samples)  
- rowdata (metadata on the genes)  
- assay_rna (a matrix of mRNA expression values as FPKM). Row names are Ensembl gene IDs and column names are sample IDs used for RNA sequencing.  
- assay_proteomics (a matrix with protein expression values as intensities). Row names are gene names. If peptides cannot be distinguished, multiple gene names are concatenated with ';'. Column names are sample IDs used in the proteomics assay.  

The file path in the config.R file has to be adjusted to point to the RData object. 
