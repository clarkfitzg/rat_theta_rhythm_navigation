# rat_theta_rhythm_navigation

Work with Jeff Calton analyzing EEG signals


Goal is to address this feedback:

"What would really help is to show whether any of the cells they were recording from were modulated by theta and whether theta within the AD nucleus was affected by the drug infusions.
That seems quite crucial to much of their data interpretation."

## Summary

We found no evidence that these cells were modulated by theta, either before or after the drug infusions.
In total, we tried four different techniques to detect a theta modulation signal.
First, we visually inspected the autocorrelograms for each data set.
__Figure 1__ shows the autocorrelograms of the EEG signal in a typical cell, both before and after the drug infusion.
The autocorrelograms are noisy, and do not show any visually detectable underlying wave structure.
Second, we computed the theta index score on the autocorrelogram using the method described by Tsanov et al. __TODO: cite appropriately__.
In this case, the parameter that represents the magnitude of the periodic component in the autocorrelograms was not significantly different from 0, which means that there's no evidence of theta modulation.
Third, we analyzed the spectral density of the signal by examining the periodogram for each data set.
The periodograms showed no structure indicative of a theta signal.
Finally, we simulated several noisy signals that contained theta modulation, and verified that all of the above methods were able to detect the theta modulation.
Hence, we do not believe these cells were modulated by theta.

__Figure 1__ From IP18.2-17-17.PreInjection.Standard.Cylinder.2017-02-17_13-11-08.SE15.nse.converted.SpikeSortedForAutocorrelation_plot.png
Caption: The autocorrelogram is noisy, and does not show any periodic structure that would indicate theta rhythm.

[right_half_auto.png]()
