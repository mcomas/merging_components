TEX = tex/partition-example-pars.tex

all : merging_components.pdf $(TEX)

merging_components.pdf : $(TEX)

figures/partition-example-mixture.pdf : script.R/partition_example.R figures/partition-example-part6.pdf

figures/partition-example-part6.pdf : script.R/partition_example.R figures/partition-example-part3a.pdf 

figures/partition-example-part3a.pdf : script.R/partition_example.R figures/partition-example-part3b.pdf

figures/partition-example-part3b.pdf : script.R/partition_example.R figures/partition-example-part5.pdf 

figures/partition-example-part5.pdf : script.R/partition_example.R figures/partition-example-part4.pdf

figures/partition-example-part4.pdf : script.R/partition_example.R figures/partition-example-part2.pdf

figures/partition-example-part2.pdf : script.R/partition_example.R figures/partition-example-part1.pdf

figures/partition-example-part1.pdf : script.R/partition_example.R tex/partition-example-pars.tex

tex/partition-example-pars.tex : script.R/partition_example.R
	Rscript $<
