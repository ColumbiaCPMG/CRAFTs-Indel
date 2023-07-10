seqDB = fread("")

fam = fread("")

exomekits = as.data.frame(table(fam$V8))

merge = merge(seqDB, fam, by.x = "sample_internal_name", by.y = "V1") 
select_broadphenotype = distinct(merge %>% select("sample_internal_name", "BroadPhenotype"))


broadpheno = as.data.frame(table(select_broadphenotype$BroadPhenotype))
broadpheno


