#!/usr/bin/env Rscript
# this is version 4.2.2
{
  rm(list=ls())
  graphics.off()
  gc()

  for(lib in c("seqinr", "jackalope", "optparse", "Biostrings", "dplyr", "ape")){
    suppressPackageStartupMessages(library(lib, character.only = T))
  }
  cat(paste0("Packages loaded...", "\n"))
} #start
?jackalope::haps_theta()

{
  args = commandArgs(trailingOnly=TRUE)
  option_list = list(
    make_option(c("--reference"), type="character", default=NULL, 
                help="reference sequence to mutate (fasta)", metavar="character"),
    make_option(c("--transposon"), type="character", default=NULL, 
                help="transposon positions in the reference sequence (bed)", metavar="character"),
    make_option(c("--nonhomoseq_pool"), type="character", default=NULL, 
                help="a fasta file with the sequence to pick from to generate a non homologuous sequence", metavar="character"),
    make_option(c("--duplication_region_size"), type="numeric", default=150000, 
                help="size of the duplicated region", metavar="character"),
    make_option(c("--haplotype_count"), type='numeric', default=10,
                help="number of haplotype to compute (>=1)"),
    make_option(c("--pop_count"), type='numeric', default=2,
                help="number of haplotype to compute (>=1)"),
    make_option(c("--structural_variants"), type='numeric', default=2,
                help="number of structural variants"),
    make_option(c("--deletion_count"), type='numeric', default=3,
                help="number of deletion region to remove"),
    make_option(c("--unmuted"), action = "store_true", default = TRUE,
                help="number of haplotype to compute (>=1)"),
    make_option(c("-p", "--prefix"), type='character', default="my_run",
                help="artificial genome name prefix"),
    make_option(c("--mutation_rate"), type='character', default=as.numeric(1.23e-7),
                help="the genome wide mutation rate"),
    make_option(c("--effective_pop"), type='character', default=as.numeric(50000e10),
                help="the effective population size"),
    make_option(c("--TCAGI"), type = "character", default = "0.172,0.329,0.172,0.329,0.99978",
                help="the rate of the different nucleotide (TCAG) + invariants (I)"),
    make_option(c("--ABCDEF"), type = "character", default = "0.65,0.05,0.21,0.27,0.02,0.65",
                help="the parameters of the GTR model"),
    make_option(c("--indel_scaling_factor"), type = "character", default=as.numeric(5),
                help="indel/snp scaling factor"),
    make_option(c("--slope"), type='numeric', default=300,
                help="Size of the slope between two structural variants"),
    make_option(c("--threads"), type='numeric', default=4,
                help="Number of threads used for fastq generation"),
    make_option(c("-o", "--output"), type='character', default="maketube_results",
                help="output folder name")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  text2numeric <- function(x) {eval(parse(text=as.character(x)))}
  
  opt$mutation_rate = text2numeric(opt$mutation_rate)
  opt$effective_pop = text2numeric(opt$effective_pop)
  opt$indel_scaling_factor = text2numeric(opt$indel_scaling_factor)
  
  opt$TCAGI <- sapply(strsplit(opt$TCAGI, ",")[[1]], text2numeric, USE.NAMES = F)
  opt$ABCDEF <- sapply(strsplit(opt$ABCDEF, ",")[[1]], text2numeric, USE.NAMES = F)
  
  if (is.null(opt$reference)){
    print_help(opt_parser)
    stop("A reference sequence must be supplied", call. = FALSE)
  } else if (is.null(opt$transposon)){
    stop("A bed file with transposons positions must be provided", call. = FALSE)
  }
  
  dir.create(opt$output, showWarnings = FALSE)
} #input parameters

debug = F

if(debug){
  SV_set = 1
  opt = list()
  opt$reference = "REF/H37Rv.fasta"
  opt$transposon = "REF/H37Rv_transposon.bed"
  opt$nonhomoseq_pool = "REF/nonH37Rv_pool_sequence_cp.fasta"
  opt$duplication_region_size = 10000
  opt$haplotype_count = 10
  opt$pop_count = 2
  opt$structural_variants = 3
  opt$deletion_count = 3
  opt$unmuted = F
  opt$prefix = "no_prefix"
  opt$mutation_rate = 1.23e-7
  opt$effective_pop = 2000
  opt$TCAGI = c(0.172, 0.329, 0.172, 0.329, 0)
  opt$ABCDEF = c(0.65, 0.05, 0.21, 0.27, 0.02, 0.65)
  opt$indel_scaling_factor = 0.125
  opt$slope = 250
  opt$threads = 1
  opt$output = "debug"
}


{
  reference_sequence <- list()
  
  reference_sequence.seqinr <- read.fasta(file = opt$reference, set.attributes = T)
  reference_sequence$sequence <- reference_sequence.seqinr[[1]]
  
  #setting a double coordinate system with vector names
  #names == original position != actual nucleotide position
  names(reference_sequence$sequence) <- as.character(1:length(reference_sequence$sequence))
  chromosome_name <- unlist(strsplit(names(reference_sequence.seqinr), '\t'))[1]
  
  transposons_loci <- apply(read.table(opt$transposon, header = F, sep = "\t", col.names = c("chrom", "start", "stop")), 1, function(x) as.list(sapply(x, trimws) ))
  
  one_big_pool <- readLines(opt$nonhomoseq_pool)
  kmer_size <- length(strsplit(one_big_pool[1], "")[[1]])
  
  #seq2pickFromNonHom <- sapply(readLines(opt$nonhomoseq), strsplit, "", USE.NAMES = F)
  
  #convert 0 based bed file to 1 based VCF
  for(i in 1:length(transposons_loci)){
    transposons_loci[[i]]$start <- min(as.numeric(transposons_loci[[i]]$start), as.numeric(transposons_loci[[i]]$stop)) + 1
    transposons_loci[[i]]$stop <- max(as.numeric(transposons_loci[[i]]$start), as.numeric(transposons_loci[[i]]$stop))
  }
  
  
} #initialisation

for(SV_set in 1:opt$structural_variants){
  SV_name <- paste0("SV", SV_set)
  SV_folder_path = paste0(opt$output, "/", SV_name)
  
  cat(paste0("Starting transposition process for ", SV_name, "\n"))
  dir.create(SV_folder_path, recursive = T, showWarnings = F)
  {
    #setting a vector of position that can be picked out
    #F = free //T = either used by a deletion region or a transposon jump
    neo_sequence <- list() ; neo_sequence$sequence <- reference_sequence$sequence
    range_to_pick <- rep(0, length(neo_sequence$sequence))
    
    #setting the initial position of transposons as unpickable
    for(i in 1:length(transposons_loci)){range_to_pick[seq(transposons_loci[[i]]$start, transposons_loci[[i]]$stop)] <- 1}
    
    #setting the post-jump transposon positions
    #object transposon jump is a list with :
    #transposon_jumps.partition[[index]]$start // transposon_jumps.partition[[index]]$stop => initial position
    #transposon_jumps.partition[[index]]$new_start // transposon_jumps.partition[[index]]$new_stop => post jump position
    
    {
      transposon_jumps.partition <- list()
      for(i in 1:length(transposons_loci)){
        #until post-jump position does not overlap with an unpickable nucleotide, retry
        #does the trick for a low number of mutation
        repeat{
          insertion_new_start <- sample(1:length(neo_sequence$sequence), 1)
          insertion_new_stop <- insertion_new_start + abs(transposons_loci[[i]]$stop - transposons_loci[[i]]$start)
          insertion_sites <- seq(insertion_new_start, insertion_new_stop)
          slopped_insertion_site <- seq(insertion_new_start - 300, insertion_new_stop + 300)
          #do{pick a transposon} while{transposon new positions overlaps with unpickable positions}
          if(sum(range_to_pick[slopped_insertion_site] == 1) == 0) {break}
        }
        
        #when done, registers the jump
        transposon_jumps.partition[[i]] <- transposons_loci[[i]]
        transposon_jumps.partition[[i]]$stop <- transposons_loci[[i]]$stop
        transposon_jumps.partition[[i]]$class <- "transposon"
        transposon_jumps.partition[[i]]$new_start <- insertion_new_start
        transposon_jumps.partition[[i]]$new_stop <- insertion_new_stop
        range_to_pick[slopped_insertion_site] <- range_to_pick[slopped_insertion_site] + 1
      }
      
      deletion_count <- opt$deletion_count
      deletions.partition <- list()
      for(i in seq(1, deletion_count)){
        #until deletions.partition are completely independant, pick another (position and size)
        repeat{
          deletion <- sample(1:length(neo_sequence$sequence), size = 1)
          deletion_size <- floor(rgamma(1, shape = 1, scale = 3500))
          deleted_nucleotides <- seq(deletion, deletion + deletion_size)
          slopped_deletion_site <- seq(deletion - opt$slope, deletion + deletion_size + opt$slope)
          if(all(range_to_pick[slopped_deletion_site] == 0)) {break}
        }
        #remove positions of the deletion from the positions to pick from so it does not overlap
        #easier to compute SNP positions after
        range_to_pick[slopped_deletion_site] <- range_to_pick[slopped_deletion_site] + 1
        deletions.partition[[i]] <- list()
        deletions.partition[[i]]$chrom <- chromosome_name
        deletions.partition[[i]]$start <- deletion - 1
        deletions.partition[[i]]$stop <- deletion + deletion_size
        deletions.partition[[i]]$class <- "deletion"
        deletions.partition[[i]]$new_start <- as.character(deletion - 1)
        deletions.partition[[i]]$new_stop <- as.character(deletion)
      }
      
      #this one works a bit differently
      duplication_region_window_start <- sample(1:(length(range_to_pick)-opt$duplication_region_size), size = 1)
      range_to_pick[duplication_region_window_start:duplication_region_window_start+35000] = 1
      
      nohomoseq.list = list(info = list(cumulative_size = 0, wanted_size = sample(4400000*seq(0.002, 0.008, by = 0.00001), 1)))
      
      i = 1
      while(nohomoseq.list[["info"]]$cumulative_size <= nohomoseq.list[["info"]]$wanted_size){
        new_sequence_size <- floor(rgamma(1, shape = 1, scale = 3500))
        insertion_index <- paste("insertion", i, new_sequence_size, sep = "_")
        
        nohomoseq.list[[insertion_index]] = list(
          individual_size = new_sequence_size,
          kmer_count = floor(new_sequence_size / kmer_size)
        )
        nohomoseq.list[[insertion_index]]$kmer_index <- sample(1:length(one_big_pool), size = nohomoseq.list[[insertion_index]]$kmer_count, replace = F)
        nohomoseq.list[[insertion_index]]$sequence <- paste(one_big_pool[nohomoseq.list[[insertion_index]]$kmer_index], collapse = "")
        
        nohomoseq.list[[insertion_index]]$size = nchar(nohomoseq.list[[insertion_index]]$sequence)
        
        nohomoseq.list[["info"]]$cumulative_size = nohomoseq.list[["info"]]$cumulative_size + nohomoseq.list[[insertion_index]]$size
        i = i + 1
      }
      
      nohomoseq.partition <- list()
      for(non_hologuous_index in names(tail(nohomoseq.list, -1))){
        insertion_site <- sample(which(range_to_pick == 0), 1)
        slopped_insertion_site <- seq(insertion_site - opt$slope, insertion_site + opt$slope)
        
        nohomoseq.partition[[non_hologuous_index]] <- list()
        nohomoseq.partition[[non_hologuous_index]]$chrom <- chromosome_name
        nohomoseq.partition[[non_hologuous_index]]$start = insertion_site
        nohomoseq.partition[[non_hologuous_index]]$stop = insertion_site+1
        nohomoseq.partition[[non_hologuous_index]]$new_start = insertion_site
        nohomoseq.partition[[non_hologuous_index]]$new_stop = insertion_site + nchar(nohomoseq.list[[non_hologuous_index]]$sequence)
        nohomoseq.partition[[non_hologuous_index]]$sequence = nohomoseq.list[[non_hologuous_index]]$sequence
        nohomoseq.partition[[non_hologuous_index]]$class = "insertion"
      }
    } #END OF PARTITION CREATION
    
    
    {
      #insertion of non homologuous sequences	
      for(insertion_nohomoseq in nohomoseq.partition){
        insertion_site <- as.numeric(max(which(names(neo_sequence$sequence) %in% insertion_nohomoseq$start)))
        neo_sequence$sequence = c(neo_sequence$sequence[1:insertion_site-1], insertion_nohomoseq$sequence, neo_sequence$sequence[insertion_site:length(neo_sequence$sequence)])
      }

      #removal of deletion regions
      nucleotides_to_delete <- unique(unlist(sapply(deletions.partition, function(del) as.character(seq(del$start, del$stop)))))
      neo_sequence$sequence <- neo_sequence$sequence[!names(neo_sequence$sequence) %in% nucleotides_to_delete]
      
      #removal of transposon initial position
      nucleotides_to_delete <- unique(unlist(sapply(transposons_loci, function(transpo_jump) as.character(seq(min(transpo_jump$start, transpo_jump$stop), max(transpo_jump$start, transpo_jump$stop))))))
      neo_sequence$sequence <- neo_sequence$sequence[!names(neo_sequence$sequence) %in% nucleotides_to_delete]
      
      #insertion of transposons final positions
      transposon_sequence <- NULL
      for(transpo_jump in transposon_jumps.partition){
        #min to selecte the position in the non duplicated area
        new_start_position <- min(which(names(neo_sequence$sequence) == transpo_jump$new_start))
        transposon_sequence <- reference_sequence$sequence[seq(min(transpo_jump$start, transpo_jump$stop), max(transpo_jump$start, transpo_jump$stop))]
        neo_sequence$sequence <- c(neo_sequence$sequence[1:new_start_position-1], transposon_sequence, neo_sequence$sequence[new_start_position:length(neo_sequence$sequence)])
      }

      #duplication region 
      neo_sequence$sequence <- c(neo_sequence$sequence, reference_sequence$sequence[seq(duplication_region_window_start, duplication_region_window_start+opt$duplication_region_size)])
      
    } #PLAYING THE PARTITION
    
    {
      transposons.info <- sapply(transposon_jumps.partition, function(x) {
        previous_start <- min(as.numeric(which(names(neo_sequence$sequence) == x$start - 1)))
        pos_transpo_start <- min(as.numeric(which(names(neo_sequence$sequence) == x$new_start)))
        return(list(
          c(chromosome_name, previous_start - opt$slope, previous_start + opt$slope, "IS_scar_site"),
          c(chromosome_name, pos_transpo_start - opt$slope, pos_transpo_start + opt$slope, "IS_flanking_region"))
        )
      })
      transposons.info <- do.call(rbind, transposons.info)
      
      deletions.info <- as.data.frame(do.call(rbind, deletions.partition))
      deletions.info <- cbind(deletions.info$chrom, as.numeric(deletions.info$new_start) - opt$slope, as.numeric(deletions.info$new_start) + opt$slope, "DR_scar_site")
      
      duplication.info <- rbind(c(chromosome_name, duplication_region_window_start, duplication_region_window_start + opt$duplication_region_size, "duplicated_region"), c(chromosome_name, length(neo_sequence$sequence) - opt$duplication_region_size, length(neo_sequence$sequence), "duplicata_region"))
      
      non_homologuous.info <- do.call(rbind, nohomoseq.partition)
      non_homologuous.info <- apply(non_homologuous.info[, c("chrom", "start", "stop", "class")], 2, as.character)
      #apply((t(sapply(nohomoseq.partition, function(x) c(x$chrom, x$start - opt$slope, x$start + opt$slope, x$class)))), 2, as.character)
      
    } #ANNOTATING THE NEW GENOME
    
    nonhomologsequences.df <- as.data.frame(t(do.call(cbind, nohomoseq.partition)))
    nonhomologsequences.df <- nonhomologsequences.df[, !names(nonhomologsequences.df) %in% c("sequence")]
    
    bed_annotation <- apply(rbind(transposons.info, deletions.info, duplication.info, non_homologuous.info), 2, as.character)
    results <- bind_rows(as.data.frame(t(do.call(cbind, deletions.partition))), as.data.frame(t(do.call(cbind, transposon_jumps.partition))), nonhomologsequences.df)
    
    results$new_start <- sapply(results$new_start, function(x) as.numeric(x) - 1)
    a = length(reference_sequence$sequence)
    b = sum(sapply(deletions.partition, function(del) abs(del$stop - del$start + 1)))
    c = length(neo_sequence$sequence)
    results <- apply(results, 2, as.character)
    cat(paste(SV_name, "\t", "Reference genome length :", a, "\n"))
    cat(paste(SV_name, "\t", "Total DR length:", b, "\n"))
    cat(paste(SV_name, "\t", "Artificial genome length (w/o indels):", c, "\n"))
    
    write.table(bed_annotation, file = paste0(SV_folder_path, "/", SV_name, "_SV.bed"), na = "", append = F, sep = "\t", row.names = F, col.names = F, quote = F)		
    write.table(results, file = paste(SV_folder_path, "/", SV_name, "_equivalence.bed", sep = ""), na = "", append = F, sep = "\t", row.names = F, col.names = F, quote = F)
    cat(paste0("Structural variants added, writing equivalence table...", "\n"))
  } # structural variation
  
  
  write.fasta(neo_sequence$sequence, file.out = paste0(SV_folder_path, "/", opt$prefix, "_unmuted.fasta"), names = chromosome_name)
  ref <- read_fasta(paste0(SV_folder_path, "/", opt$prefix, "_unmuted.fasta"), cut_names = T)
  
  if(!opt$unmuted) {
    unlink(paste(SV_folder_path, "/", opt$prefix, "_unmuted.fasta", sep = ""))
  }
  
  {
    output_population_table <- NULL
    for(population_index in 1:opt$pop_count){
      pop_name <- paste0("pop", population_index)
      
      pop_folder <- paste0(SV_folder_path, "/pop", population_index)
      dir.create(pop_folder, showWarnings = FALSE)
      
      GTR_tuberculosis <- sub_GTR(pi_tcag = opt$TCAGI[1:4], abcdef = opt$ABCDEF, invariant = opt$TCAGI[5], mu = opt$mutation_rate)
      tb_theta <- 2*(opt$effective_pop)*opt$mutation_rate
      
      deletion_model <- indels(rate = opt$indel_scaling_factor*opt$mutation_rate/2, max_length = 90)
      insertion_model <- indels(rate = opt$indel_scaling_factor*opt$mutation_rate/2, max_length = 45)
      
      pop_theta <- haps_theta(theta = tb_theta, n_haps = opt$haplotype_count)
      
      one_population_object <- create_haplotypes(ref,
                                                 haps_info = pop_theta,
                                                 sub = GTR_tuberculosis,
                                                 ins = insertion_model, del = deletion_model)
      
      one_population_object$set_names(paste0("H", 1:opt$haplotype_count))
      cat(paste("Sequence mutated, writing fasta, fastq and vcf for population", population_index, "of", SV_name, "\n"))
      
      output_tree <- pop_theta$phylo()
      output_tree$tip.label <- sapply(output_tree$tip.label, function(y) gsub(pattern = "t", replacement = "H", x = y))
      
      write.tree(output_tree, paste0(pop_folder, "/", SV_name, "_", pop_name, ".nwk"))
      write_vcf(one_population_object, paste(pop_folder, "/", SV_name, "_", pop_name, sep = ""), compress = T, overwrite = T)
      write_fasta(one_population_object, out_prefix=paste(pop_folder, "/FASTA/", SV_name, "_", pop_name, sep = ""), compress = T, text_width = 60, overwrite = T, n_threads = opt$threads)
      
      number_of_read=floor(as.numeric(100*4000000/150))
      
      illumina(obj = one_population_object, out_prefix = paste(pop_folder, "/FASTQ/", SV_name, "_", pop_name, "_", sep = ""),
               seq_sys = "MSv3",
               paired = T, matepair = T,
               n_reads = 5, read_length = 250,
               frag_mean = 400, frag_sd = 30,
               del_prob1 = 0, ins_prob1 = 0, del_prob2 = 0, ins_prob2 = 0, prob_dup = 0, 
               haplotype_probs = NULL,
               sep_files = T, compress = T, overwrite = T, n_threads = opt$threads)
      
      cat(paste("Renaming files for convenience", "\n"))
      
      all_fasta <- list.files(paste0(pop_folder, "/FASTA"), pattern="*.fa.gz", full.names=TRUE, recursive=FALSE)
      all_fastq <- list.files(paste0(pop_folder, "/FASTQ"), pattern="*.fq.gz", full.names=TRUE, recursive=FALSE)
      
      invisible(
        lapply(all_fasta, function(fasta) {
          new_name=gsub(pattern = "__", replacement = "_", x = fasta)
          file.rename(fasta, new_name)
        })
      )
      invisible(
        lapply(all_fastq, function(fastq) {
          new_name=sub(pattern = "__", replacement = "_", x = fastq)
          file.rename(fastq, new_name)
        })
      )
    }
  }
}



# A LOT OF CIRCONVOLUTION TO MAKE SURE EVERYTHING MATCHES NO MATTER THE ORDER OF THE FILES

fasta		<-	list.files(opt$output, pattern = "*.fa.gz", full.names = T, recursive = T)
fasta_names <- sapply(gsub(x = fasta, ".fa.gz", ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
fasta.df <- data.frame(strain = fasta_names, fasta = fasta)

fastq_R1	<-	list.files(opt$output, pattern = "*R1.fq.gz", full.names = T, recursive = T)
fastq_R1_names <- sapply(gsub(x = fastq_R1, ".R1.fq.gz", ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
fastq_R1.df <- data.frame(strain = fastq_R1_names, R1 = fastq_R1)

fastq_R2	<-	list.files(opt$output, pattern = "*R2.fq.gz", full.names = T, recursive = T)
fastq_R2_names <- sapply(gsub(x = fastq_R2, ".R2.fq.gz", ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
fastq_R2.df <- data.frame(strain = fastq_R2_names, R2 = fastq_R2)

newick <- list.files(opt$output, pattern = "*.nwk", full.names = T, recursive = T)
newick_SV_POP <- sapply(gsub(x = newick, pattern = ".nwk", replacement = ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
newick_SV <- sapply(strsplit(newick_SV_POP, "_"), function(x) x[1])
newick_POP <- sapply(strsplit(newick_SV_POP, "_"), function(x) x[2])
newick.df <- data.frame(SV = newick_SV, POP = newick_POP, NWK = newick)

vcf <- list.files(opt$output, pattern = "*.vcf.gz", full.names = T, recursive = T)
vcf_SV_POP <- sapply(gsub(x = vcf, pattern = ".vcf.gz", replacement = ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
vcf_SV <- sapply(strsplit(vcf_SV_POP, "_"), function(x) x[1])
vcf_POP <- sapply(strsplit(vcf_SV_POP, "_"), function(x) x[2])
vcf.df <- data.frame(SV = vcf_SV, POP = vcf_POP, VCF = vcf)

equivalence_file <- list.files(opt$output, pattern = "*_equivalence.bed", full.names = T, recursive = T)
equivalence_SV <- sapply(gsub(x = equivalence_file, pattern = "_equivalence.bed", replacement = ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
equivalence.df <- data.frame(equivalence = equivalence_file, SV = equivalence_SV)

structural_variant_file <- list.files(opt$output, pattern = "*_SV.bed", full.names = T, recursive = T)
SV_SV <- sapply(gsub(x = structural_variant_file, pattern = "_SV.bed", replacement = ""), function(x) {tail(strsplit(x, "/")[[1]], 1)}, USE.NAMES = F)
SV.df <- data.frame(SV = SV_SV, SV_file = structural_variant_file)

output_population_table <- merge(fasta.df, fastq_R1.df, by = "strain")
output_population_table <- merge(output_population_table, fastq_R2.df, by = "strain")
output_population_table$SV <- sapply(strsplit(output_population_table$strain, "_"), function(x) x[1])
output_population_table$POP <- sapply(strsplit(output_population_table$strain, "_"), function(x) x[2])
output_population_table$HAP <- sapply(strsplit(output_population_table$strain, "_"), function(x) x[3])

output_population_table <- merge(output_population_table, newick.df, by = c("SV", "POP"))
output_population_table <- merge(output_population_table, vcf.df, by = c("SV", "POP"))
output_population_table <- merge(output_population_table, equivalence.df, by = "SV")
output_population_table <- merge(output_population_table, SV.df, by = "SV")

order <- c("strain", "HAP", "SV", "POP", "NWK", "equivalence", "SV_file", "VCF", "fasta", "R1", "R2")
output_population_table <- output_population_table[, order]
write(file = paste0(opt$output, "/", opt$prefix, "_arborescence.tsv"), c("strain\tHAP\tSV\tpopulation\ttree\tequivalence\tSV_file\tVCF\tfasta\tR1\tR2"))
write.table(output_population_table, paste0(opt$output, "/", opt$prefix, "_arborescence.tsv"), sep = "\t", col.names = F, row.names = F, quote = F, append = T)


cat(paste("All done :)", "\n"))
