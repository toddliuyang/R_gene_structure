library("GeneStructurePlot")

##Read the information of gene "stambpl1" from file "gene_structure_table_stambpl1.txt",recalculate the position and plot.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_stambpl1.txt"))
gene_info<-gene_data_format(gene_info)
plot_gene_structure_table (gene_info)

## Plot the structure of gene "gene-stambpl1" in file "gff3_sample.gff3" after re-calculating the position.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_info_extract(file.path(dir,"gff3_sample.gff3"),"gene-stambpl1")
gene_info<-gene_data_format(gene_info)
plot_gene_structure_table (gene_info)

##Obtain the information of gene "gene-stambpl1" for plotting its structure.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_info_extract(file.path(dir,"gff3_sample.gff3"),"gene-stambpl1")

##Read the information of gene "stambpl1"  from file "gene_structure_table_stambpl1.txt" for plotting .
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_stambpl1.txt"))
plot_gene_structure_table (gene_info)

##Read the multi genes information from text, format the data, and then plot. 
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_multi_gene_genomic_pos.txt"))
gene_info<-multi_genes_data_format(gene_info,FALSE)
plot_multi_gene_structure(gene_info)

##Read multi genes information and format, then plot.
dir <- system.file("extdata", package="GeneStructurePlot")
multi_genes_info<-multi_genes_data_merge(file.path(dir,"gff3_sample.gff3"),file.path(dir,"gene_ids.txt"))
plot_multi_gene_structure(multi_genes_info)

##Read one gene information and format, then plot.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-multi_genes_data_merge(file.path(dir,"gff3_sample.gff3"),file.path(dir,"gene_ids3.txt"))
plot_gene_structure_table (gene_info)

## Plot the structure of gene "gene-LOC101070858" in file "gff3_sample.gff3".
dir <- system.file("extdata", package="GeneStructurePlot")
plot_gene_structure(file.path(dir,"gff3_sample.gff3"),"gene-LOC101070858")


##Read the information of gene "stambpl1"  from file "gene_structure_table_stambpl1.txt" and plot.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_stambpl1.txt"))
plot_gene_structure_table (gene_info)

## Plot the structure of gene "gene-stambpl1" in file "gff3_sample.gff3" after re-calculating the position.
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_info_extract(file.path(dir,"gff3_sample.gff3"),"gene-stambpl1")
gene_info<-gene_data_format(gene_info)
plot_gene_structure_table (gene_info)

#Plot multi genes structure based on annotation file (GFF3 or GTF),IDs were saved in "gene_ids.txt".
dir <- system.file("extdata", package="GeneStructurePlot")
genes_info<-multi_genes_data_merge(file.path(dir,"gff3_sample.gff3"),file.path(dir,"gene_ids.txt"))
plot_multi_gene_structure(genes_info)

#Plot multi genes structure based on 9-column file "gene_structure_table_multi_gene.txt".
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_multi_gene.txt"))
plot_multi_gene_structure(gene_info)

#Plot multi genes structure based on 9-column file "gene_structure_table_multi_gene_genomic_pos.txt"
gene_info<-gene_structure_table_read(file.path(dir,"gene_structure_table_multi_gene_genomic_pos.txt"))
gene_info<-multi_genes_data_format(gene_info,FALSE)
plot_multi_gene_structure(gene_info)
plot_multi_gene_structure(gene_info,FALSE)


##Obtain the number and IDs of transcripts from the specific gene "gene-stambpl1" in dataframe "gene_info".
dir <- system.file("extdata", package="GeneStructurePlot")
gene_info<-gene_info_extract(file.path(dir,"gff3_sample.gff3"),"gene-stambpl1")
transcripts_info<-trancsript_analysis_df(gene_info)

##Obtain the number and IDs of transcripts from the specific gene "gene-stambpl1" in file "gff3_sample.gff3".
dir <- system.file("extdata", package="GeneStructurePlot")
transcripts_info<-trancsript_analysis_file(file.path(dir,"gff3_sample.gff3"),"gene-stambpl1")











