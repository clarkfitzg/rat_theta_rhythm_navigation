source("functions.R")


d1 = plot_acf_pair(file_pre = "IP21.5-30-17.PreInjection.Standard.Cylinder.2017-05-30_16-07-11.SE9.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP21.5-30-17.PostInjection.Standard.Cylinder.2017-05-30_18-30-23.SE9.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip"
    , plot_func = plot_spike_locs)

d2 = plot_acf_pair(file_pre = "IP24.09-20-17.IP24.PreinjectionStandard.SE7.2017-09-20_17-13-06.SE7.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP24.09-20-17.IP24.Postinjection.Saline.Standard.SE7.2017-09-20_18-59-41.SE7.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip"
    , plot_func = plot_spike_locs)

pdf("IP25.02-06-18.pdf", width = 12, height = 12)
d3 = plot_acf_pair(file_pre = "IP25.02-06-18.IP25.Preinjection.Standard.SE5.2018-02-06_15-30-43.SE5.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP25.02-06-18.IP25.Postinjection.Muscimol.Standard.SE5.2018-02-06_17-56-13.SE5.nse.converted.SpikeSortedForAutocorrelation"
              , zipfile = "data_3-14-2020.zip"
    , plot_func = plot_spike_locs
    , ylab = "spike index locations"
)
dev.off()

pdf("IP18.12-9-16.pdf", width = 12, height = 12)
d4 = plot_acf_pair(file_pre = "IP18.12-9-16.PreInjection.Standard.Cylinder.2016-12-09_14-31-21.SE9.nse.converted.SpikeSortedForAutocorrelation"
              , file_post = "IP18.12-9-16.PostInjection.Standard.Cylinder.2016-12-09_19-18-01.SE15.nse.converted.SpikeSortedForAutocorrelation"
    , plot_func = plot_spike_locs
    , ylab = "spike index locations"
)
dev.off()


# Not seeing any signal in the spectrum for any of these.
# The only thing I see is a spike at 0, which is probably caused by the cells not firing twice in a row.
s = spectrum(d3$pre$spike, spans = 201)


