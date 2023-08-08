# script for some of the testing data used in tests/testthat/test-utils.R
# produced by using dput() and chatGPT 

# output of dput(filteringMulti(head(contig_list[[1]])))
filteringMulti_expected <- structure(list(
	barcode = c("AAACCTGAGAGCTGGT", "AAACCTGAGCATCATC", "AAACCTGAGCATCATC", "AAACCTGAGTGGTCCC"),
	is_cell = c(TRUE, TRUE, TRUE, TRUE),
	contig_id = c("AAACCTGAGAGCTGGT-1_contig_1", "AAACCTGAGCATCATC-1_contig_2", "AAACCTGAGCATCATC-1_contig_1", "AAACCTGAGTGGTCCC-1_contig_1"),
	high_confidence = c(TRUE, TRUE, TRUE, TRUE),
	length = c(705L, 567L, 693L, 593L),
	chain = c("TRB", "TRA", "TRB", "TRB"),
	v_gene = c("TRBV20-1", "TRAV12-1", "TRBV5-1", "TRBV7-9"),
	d_gene = c("TRBD1", "None", "TRBD2", "TRBD1"),
	j_gene = c("TRBJ1-5", "TRAJ37", "TRBJ2-2", "TRBJ2-5"),
	c_gene = c("TRBC1", "TRAC", "TRBC2", "TRBC2"),
	full_length = c(TRUE, TRUE, TRUE, TRUE),
	productive = c("TRUE", "TRUE", "TRUE", "TRUE"),
	cdr3 = c("CSASMGPVVSNQPQHF", "CVVNDEGSSNTGKLIF", "CASSWSGAGDGELFF", "CASSPSEGGRQETQYF"),
	cdr3_nt = c("TGCAGTGCTAGCATGGGACCGGTAGTGAGCAATCAGCCCCAGCATTTT",
				"TGTGTGGTGAACGATGAAGGCTCTAGCAACACAGGCAAACTAATCTTT",
				"TGCGCCAGCAGCTGGTCAGGAGCGGGAGACGGGGAGCTGTTTTTT",
				"TGTGCCAGCAGCCCCTCCGAAGGGGGGAGACAAGAGACCCAGTACTTC"),
	reads = c(16718L, 18297L, 26719L, 11218L),
	umis = c(6L, 6L, 11L, 6L),
	raw_clonotype_id = c("clonotype96", "clonotype97", "clonotype97", "clonotype98"),
	raw_consensus_id = c("clonotype96_consensus_1", "clonotype97_consensus_1", "clonotype97_consensus_2", "clonotype98_consensus_1")),
	class = c("grouped_df", "tbl_df", "tbl", "data.frame"),
	row.names = c(NA, -4L),
	groups = structure(list(
		barcode = c("AAACCTGAGAGCTGGT", "AAACCTGAGCATCATC", "AAACCTGAGCATCATC", "AAACCTGAGTGGTCCC"),
		chain = c("TRB", "TRA", "TRB", "TRB"),
		.rows = structure(list(1L, 2L, 3L, 4L),
						  ptype = integer(0),
						  class = c("vctrs_list_of", "vctrs_vctr", "list"))
	), row.names = c(NA, -4L), .drop = TRUE, class = c("tbl_df", "tbl", "data.frame"))
)

makeGenes_T_input <- structure(list(
	barcode = c(
		"AAACCTGAGAGCTGGT", "AAACCTGAGCATCATC", "AAACCTGAGCATCATC", 
		"AAACCTGAGTGGTCCC", "AAACCTGCAAACGCGA", "AAACCTGCAAACGCGA"
	),
	is_cell = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	contig_id = c(
		"AAACCTGAGAGCTGGT-1_contig_1", "AAACCTGAGCATCATC-1_contig_1", 
		"AAACCTGAGCATCATC-1_contig_2", "AAACCTGAGTGGTCCC-1_contig_1", 
		"AAACCTGCAAACGCGA-1_contig_1", "AAACCTGCAAACGCGA-1_contig_3"
	),
	high_confidence = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	length = c(705L, 693L, 567L, 593L, 659L, 514L),
	chain = c("TRB", "TRB", "TRA", "TRB", "TRB", "TRA"),
	v_gene = c("TRBV20-1", "TRBV5-1", "TRAV12-1", "TRBV7-9", "TRBV2", "TRAV29DV5"),
	d_gene = c("TRBD1", "TRBD2", "None", "TRBD1", "TRBD1", "None"),
	j_gene = c("TRBJ1-5", "TRBJ2-2", "TRAJ37", "TRBJ2-5", "TRBJ1-6", "TRAJ22"),
	c_gene = c("TRBC1", "TRBC2", "TRAC", "TRBC2", "TRBC1", "TRAC"),
	full_length = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	productive = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"),
	cdr3 = c(
		"CSASMGPVVSNQPQHF", "CASSWSGAGDGELFF", "CVVNDEGSSNTGKLIF", 
		"CASSPSEGGRQETQYF", "CASRVQGNRGSPLHF", "CAASGYGSARQLTF"
	),
	cdr3_nt = c(
		"TGCAGTGCTAGCATGGGACCGGTAGTGAGCAATCAGCCCCAGCATTTT", 
		"TGCGCCAGCAGCTGGTCAGGAGCGGGAGACGGGGAGCTGTTTTTT", 
		"TGTGTGGTGAACGATGAAGGCTCTAGCAACACAGGCAAACTAATCTTT", 
		"TGTGCCAGCAGCCCCTCCGAAGGGGGGAGACAAGAGACCCAGTACTTC", 
		"TGTGCCAGCAGGGTACAGGGTAATAGGGGTTCACCCCTCCACTTT", 
		"TGTGCAGCAAGCGGTTACGGTTCTGCAAGGCAACTGACCTTT"
	),
	reads = c(16718L, 26719L, 18297L, 11218L, 10891L, 874L),
	umis = c(6L, 11L, 6L, 6L, 4L, 1L),
	raw_clonotype_id = c(
		"clonotype96", "clonotype97", "clonotype97", 
		"clonotype98", "clonotype99", "clonotype99"
	),
	raw_consensus_id = c(
		"clonotype96_consensus_1", "clonotype97_consensus_2", 
		"clonotype97_consensus_1", "clonotype98_consensus_1", 
		"clonotype99_consensus_1", "clonotype99_consensus_2"
	)
), row.names = c(1L, 3L, 4L, 6L, 8L, 9L), class = "data.frame")

makeGenes_T_expected <- structure(list(
	barcode = c("AAACCTGAGAGCTGGT", "AAACCTGAGCATCATC", "AAACCTGAGCATCATC", "AAACCTGAGTGGTCCC", "AAACCTGCAAACGCGA", "AAACCTGCAAACGCGA"),
	is_cell = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	contig_id = c("AAACCTGAGAGCTGGT-1_contig_1", "AAACCTGAGCATCATC-1_contig_1", "AAACCTGAGCATCATC-1_contig_2", "AAACCTGAGTGGTCCC-1_contig_1", "AAACCTGCAAACGCGA-1_contig_1", "AAACCTGCAAACGCGA-1_contig_3"),
	high_confidence = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	length = c(705L, 693L, 567L, 593L, 659L, 514L),
	chain = c("TRB", "TRB", "TRA", "TRB", "TRB", "TRA"),
	v_gene = c("TRBV20-1", "TRBV5-1", "TRAV12-1", "TRBV7-9", "TRBV2", "TRAV29DV5"),
	d_gene = c("TRBD1", "TRBD2", "None", "TRBD1", "TRBD1", "None"),
	j_gene = c("TRBJ1-5", "TRBJ2-2", "TRAJ37", "TRBJ2-5", "TRBJ1-6", "TRAJ22"),
	c_gene = c("TRBC1", "TRBC2", "TRAC", "TRBC2", "TRBC1", "TRAC"),
	full_length = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	productive = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"),
	cdr3 = c("CSASMGPVVSNQPQHF", "CASSWSGAGDGELFF", "CVVNDEGSSNTGKLIF", "CASSPSEGGRQETQYF", "CASRVQGNRGSPLHF", "CAASGYGSARQLTF"),
	cdr3_nt = c("TGCAGTGCTAGCATGGGACCGGTAGTGAGCAATCAGCCCCAGCATTTT", "TGCGCCAGCAGCTGGTCAGGAGCGGGAGACGGGGAGCTGTTTTTT", "TGTGTGGTGAACGATGAAGGCTCTAGCAACACAGGCAAACTAATCTTT", "TGTGCCAGCAGCCCCTCCGAAGGGGGGAGACAAGAGACCCAGTACTTC", "TGTGCCAGCAGGGTACAGGGTAATAGGGGTTCACCCCTCCACTTT", "TGTGCAGCAAGCGGTTACGGTTCTGCAAGGCAACTGACCTTT"),
	reads = c(16718L, 26719L, 18297L, 11218L, 10891L, 874L),
	umis = c(6L, 11L, 6L, 6L, 4L, 1L),
	raw_clonotype_id = c("clonotype96", "clonotype97", "clonotype97", "clonotype98", "clonotype99", "clonotype99"),
	raw_consensus_id = c("clonotype96_consensus_1", "clonotype97_consensus_2", "clonotype97_consensus_1", "clonotype98_consensus_1", "clonotype99_consensus_1", "clonotype99_consensus_2"),
	TCR1 = c(NA, NA, "TRAV12-1.TRAJ37.TRAC", NA, NA, "TRAV29DV5.TRAJ22.TRAC"),
	TCR2 = c("TRBV20-1.TRBD1.TRBJ1-5.TRBC1", "TRBV5-1.TRBD2.TRBJ2-2.TRBC2", NA, "TRBV7-9.TRBD1.TRBJ2-5.TRBC2", "TRBV2.TRBD1.TRBJ1-6.TRBC1", NA)),
	row.names = c(1L, 3L, 4L, 6L, 8L, 9L),
	class = "data.frame"
)
