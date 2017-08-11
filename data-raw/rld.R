data(dds)
rld <- rlog(dds)
saveData(rld, compress = "xz")
