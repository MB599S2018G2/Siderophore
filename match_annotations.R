
library(data.table)

mega_database <- data.table("NA")
no_hit_database <- data.frame("NA")
filecounter <- 0
OF_outputs <- list.files(path="/Users/user/Box Sync/Lab-box/GROUP2_Siderophore/formatted_outORthofinder/", pattern="*.txt", full.names=T, recursive=TRUE)

# Part 1: Replace gene IDs with names from annotation

N_genomes <- ncol(OF_output)

for (this_path in OF_outputs){
  OF_output <- data.frame(fread(this_path))
  for (i in 1:N_genomes){
    print("Load new file")
    filecounter <- filecounter + 1
    this_genome_ID <- colnames(OF_output)[i]
    print(this_genome_ID)
    path_to_GTF <- paste0("/Users/user/Box\ Sync/Lab-box/GROUP2_Siderophore/gtf 2/", as.character(this_genome_ID), ".gtf")
    #this_OF_output <- this_OF_output[-1,]
    gtf_table <- fread(path_to_GTF)
    #print(paste0("Running file no.", filecounter, basename(filein)))
    length(OF_output[,1])
    a_vector <- OF_output[,i]
    
    for (j in 1:length(a_vector)){
      #print(i)
      #if(i == 1000|2000|3000|4000|5000){print((i/nrow(this_OF_output)))}
      gene <- (a_vector)[j]
      if(is.na(a_vector[1]) == "FALSE"){
        if(gene != ""){
          #gene <- as.character((OF_output)[j,i])
          extracted_rows_ortholog <- gtf_table %>% filter(str_detect(gtf_table$V9, gene))
          this_gene_annotation_info <- unlist(strsplit(as.character(extracted_rows_ortholog[9]), ";"))
          if(length(this_gene_annotation_info) < 2){
            #no_hit_database <- rbind(no_hit_database, t(data.frame(c(as.character(this_genome_ID), gene))), fill= TRUE)
          } else {
            to_append <- c(as.character(this_genome_ID), this_gene_annotation_info)
            #mega_database <- rbind(mega_database, t(data.frame(to_append)), fill=TRUE)
            to_append[5] <- gsub('"', "%", to_append[5])
            this_annotation <- (strsplit(to_append[5], "%"))[[1]][2]
            OF_output[j,i] <- this_annotation
            #print(i)
          }
        }
      }
    }
  }
  fwrite(OF_output, this_path, sep="\t")
}

# Part 2: Search for occurences of genes not in expected clusters, but in genome annotations.

Cluster_file <- data.frame(fread("OG0003091- pyoverdine biosynthesis-like protein,  pyoverdine biosynthesis protein PvdN_FAKE.txt", fill = TRUE, header = TRUE))
annotation_file <- read.table("/Users/user/Box Sync/Lab-box/GROUP2_Siderophore/formatted_outORthofinder/PvdN_annotations.txt", fill=TRUE)

has_gene_in_annotation <- (rownames(annotation_file))

strains_without_gene_in_cluster <- data.frame(matrix(NA))

for (i in 1:ncol(Cluster_file)){
  print(i)
  print("Now checking:")
  print(colnames(Cluster_file)[i])
  if (is.na(Cluster_file[1,i]) == "TRUE"){
    print("No gene. Just have:")
    print(Cluster_file[1,i])
    no_gene_in_my_cluster <- colnames(Cluster_file)[i]
    print("Does the annotation say a gene should be there?")
    if (no_gene_in_my_cluster %in% annotation_file[,1]){
      print("It's in the annotation file!")
      strains_without_gene_in_cluster <- rbind(strains_without_gene_in_cluster, no_gene_in_my_cluster)
      print("THIS is what we're adding to output file:")
      print(no_gene_in_my_cluster)
    }
    else{
      print("Gene not in annotation file")
      print("")
    }
  }
  else {
    print("Gene is in the cluster!")
    print("")
  }
  #if (colnames(Cluster_file[i]) == "Pseudomonas_koreensis_D26_4042"){
   # break
  #}
}

fwrite(strains_without_gene_in_cluster, "strains_with_annotated_PvdF_but_not_in_cluster.txt")



