seqDB = fread("/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/2023-05-31_seqDB_exome_genome_in_dragendb.csv")

fam = fread("/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/final_dev/2023-03-24_11-47-14_IGM_n39367_indels_existing.sample.txt")

#### Table S1: Exome Kits
exomekits = as.data.frame(table(fam$V8))
exomekits 

merge = merge(seqDB, fam, by.x = "sample_internal_name", by.y = "V1") 
select_broadphenotype = distinct(merge %>% select("sample_internal_name", "BroadPhenotype"))



#### Table S2: Broad Phenotypes 
broadpheno = as.data.frame(table(select_broadphenotype$BroadPhenotype))
broadpheno


