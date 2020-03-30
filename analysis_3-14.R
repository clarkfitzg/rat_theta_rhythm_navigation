source("functions.R")

par(mfrow = c(2, 1))

plot_acf_pair(file_pre = "IP21.5-30-17.PreInjection.Standard.Cylinder.2017-05-30_16-07-11.SE9.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP21.5-30-17.PostInjection.Standard.Cylinder.2017-05-30_18-30-23.SE9.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip")

plot_acf_pair(file_pre = "IP24.09-20-17.IP24.PreinjectionStandard.SE7.2017-09-20_17-13-06.SE7.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP24.09-20-17.IP24.Postinjection.Saline.Standard.SE7.2017-09-20_18-59-41.SE7.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip")

plot_acf_pair(file_pre = "IP25.02-06-18.IP25.Preinjection.Standard.SE5.2018-02-06_15-30-43.SE5.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP25.02-06-18.IP25.Postinjection.Muscimol.Standard.SE5.2018-02-06_17-56-13.SE5.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip")
