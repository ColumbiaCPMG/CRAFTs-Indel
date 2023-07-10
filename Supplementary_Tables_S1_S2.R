seqDB = fread("")

fam = fread("")

## supplementary table S1: exome kit
exomekits = as.data.frame(table(fam$V8))
exomekits

merge = merge(seqDB, fam, by.x = "sample_internal_name", by.y = "V1") 
select_broadphenotype = distinct(merge %>% select("sample_internal_name", "BroadPhenotype"))

## supplementary table S2: broad phenotypes 
broadpheno = as.data.frame(table(select_broadphenotype$BroadPhenotype))
broadpheno


