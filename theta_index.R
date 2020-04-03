source("functions.R")

f1 = unzip("data.zip", list = TRUE)[, "Name"]

r1 = lapply(f1, file_name_and_theta_index)

z2 = "data_3-14-2020.zip"
f2 = unzip(z2, list = TRUE)[, "Name"]
f2 = f2[grepl("spike", f2, ignore.case = TRUE)]
r2 = lapply(f2, file_name_and_theta_index, zipfile = z2)

r = c(r1, r2)
r = do.call(rbind, r)
row.names(r) = NULL
r = r[, c("theta_index", "file")]

write.csv(r, "theta_index.csv")
