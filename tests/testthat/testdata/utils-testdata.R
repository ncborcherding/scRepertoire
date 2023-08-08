# output of dput(filteringMulti(head(contig_list[[1]])))
filteringMulti_expected <- structure(list(barcode = c("AAACCTGAGAGCTGGT", "AAACCTGAGCATCATC", 
						   "AAACCTGAGCATCATC", "AAACCTGAGTGGTCCC"), is_cell = c(TRUE, TRUE, 
						   													 TRUE, TRUE), contig_id = c("AAACCTGAGAGCTGGT-1_contig_1", "AAACCTGAGCATCATC-1_contig_2", 
						   													 						   "AAACCTGAGCATCATC-1_contig_1", "AAACCTGAGTGGTCCC-1_contig_1"), 
			   high_confidence = c(TRUE, TRUE, TRUE, TRUE), length = c(705L, 
			   														567L, 693L, 593L), chain = c("TRB", "TRA", "TRB", "TRB"), 
			   v_gene = c("TRBV20-1", "TRAV12-1", "TRBV5-1", "TRBV7-9"), 
			   d_gene = c("TRBD1", "None", "TRBD2", "TRBD1"), j_gene = c("TRBJ1-5", 
			   														  "TRAJ37", "TRBJ2-2", "TRBJ2-5"), c_gene = c("TRBC1", "TRAC", 
			   														  											"TRBC2", "TRBC2"), full_length = c(TRUE, TRUE, TRUE, TRUE
			   														  											), productive = c("TRUE", "TRUE", "TRUE", "TRUE"), cdr3 = c("CSASMGPVVSNQPQHF", 
			   														  																										"CVVNDEGSSNTGKLIF", "CASSWSGAGDGELFF", "CASSPSEGGRQETQYF"
			   														  											), cdr3_nt = c("TGCAGTGCTAGCATGGGACCGGTAGTGAGCAATCAGCCCCAGCATTTT", 
			   														  														   "TGTGTGGTGAACGATGAAGGCTCTAGCAACACAGGCAAACTAATCTTT", "TGCGCCAGCAGCTGGTCAGGAGCGGGAGACGGGGAGCTGTTTTTT", 
			   														  														   "TGTGCCAGCAGCCCCTCCGAAGGGGGGAGACAAGAGACCCAGTACTTC"), reads = c(16718L, 
			   														  														   															   18297L, 26719L, 11218L), umis = c(6L, 6L, 11L, 6L), raw_clonotype_id = c("clonotype96", 
			   														  														   															   																		 "clonotype97", "clonotype97", "clonotype98"), raw_consensus_id = c("clonotype96_consensus_1", 
			   														  														   															   																		 																   "clonotype97_consensus_1", "clonotype97_consensus_2", "clonotype98_consensus_1"
			   														  														   															   																		 )), class = c("grouped_df", "tbl_df", "tbl", "data.frame"
			   														  														   															   																		 ), row.names = c(NA, -4L), groups = structure(list(barcode = c("AAACCTGAGAGCTGGT", 
			   														  														   															   																		 															   "AAACCTGAGCATCATC", "AAACCTGAGCATCATC", "AAACCTGAGTGGTCCC"), 
			   														  														   															   																		 												   chain = c("TRB", "TRA", "TRB", "TRB"), .rows = structure(list(
			   														  														   															   																		 												   	1L, 2L, 3L, 4L), ptype = integer(0), class = c("vctrs_list_of", 
			   														  														   															   																		 												   												   "vctrs_vctr", "list"))), row.names = c(NA, -4L), .drop = TRUE, class = c("tbl_df", 
			   														  														   															   																		 												   												   																		 "tbl", "data.frame")))
