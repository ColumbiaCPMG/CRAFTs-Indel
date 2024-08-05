seqDB = fread("2023-05-31_seqDB_exome_genome_in_dragendb.csv")
## Using seqDB pull from 5/31/2023

fam = fread("2023-03-24_11-47-14_IGM_n39367_indels_existing.sample.txt")
## Using FAM file generated from 3/24/2023

## supplementary table S1: exome kit
exomekits = as.data.frame(table(fam$V8, fam$V7))
exome_only = exomekits %>% filter(Var2 == "exome")
exome_only
genome_only = exomekits %>% filter(Var2 != "exome")
genome_only


merge = merge(seqDB, fam, by.x = "sample_internal_name", by.y = "V1") 
select_broadphenotype = distinct(merge %>% select("sample_internal_name", "BroadPhenotype"))

## supplementary table S2: broad phenotypes 
broadpheno = as.data.frame(table(select_broadphenotype$BroadPhenotype))
broadpheno


