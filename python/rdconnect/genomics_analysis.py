def aggregate_by_gene(hl, sourcePath, destinationPath, nPartitions):
	dataset= hl.read_matrix_table(sourcePath)
	#filter on samples
	new_vcf=new_vcf.filter_entries((new_vcf.sample.gt=='0/1') & (new_vcf.sample.dp > 7) & (new_vcf.sample.gq > 19))
	#filter on root level fields (chrom,pos, populations,clinvar)
	dataset=dataset.explode_rows(dataset.effs)
	#filter on effects
	dataset=dataset.distinct_by_row()
	dataset=dataset.group_rows_by(dataset.effs.gene_name).aggregate(mac = hl.agg.counter(dataset.sample.sample))
	dataset.write(destinationPath,overwrite=True)