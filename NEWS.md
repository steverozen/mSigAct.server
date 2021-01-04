# mSigAct.server 0.0.8.90xx
## Changed
* Only show spectra plots for catalog types which have mutations.

* Only show catalog types for signature attribution which have mutations.

* Changed header information in run-information.txt from generated zip file.

* Preserved the file path extensions for 

<br>

# mSigAct.server 0.0.7
## Added
* Added thumbnail picture for selected sample in "Get signature attributions" tab.

## Changed
* Made it default to round the reconstructed catalog for plotting.

## Fixed
* Fixed a bug when user uploads a VCF with "unknown" region and proceed to do signature 
attribution.

<br>

# mSigAct.server 0.0.6
## Changed
* Don't preselect catalog type for signature attribution.

* Hide the "Choose catalog type" UI after using preloaded spectra.

* Removed widget for use to specify the downloadable zip file name.

* Provided a download button for downloading the .zip file generated from
uploaded VCF.

* Automatically hide the "Show spectra", "Get signature attributions" and
"Results" tab when user uploads new VCF/spectra or run example analysis of VCF
or use preloaded spectra.

* Changed the "Search:" label to "Search in signatures ID and  etiologies:".

* Changed plan to resolve a future to multicore(not supported on Windows).

* More informative errors for malformed catalog spectra.

## Fixed
* Fixed a bug when mutation type change after attribution is not working.

* Fixed a bug of rounding reconstructed catalog.

<br>

# mSigAct.server 0.0.5
## Added
* Added thumbnail pictures for signature attribution results.

* Enabled process VCFs with unknown variant caller.

* Created internal data `COSMIC.v3.sigs` to select the signatures used for 
attribution based on user input.

## Changed
* Create new temporary directory each time when user generate catalogs and do
signature attribution.

## Fixed
* Fixed a bug in generating zip archive when the input region is "unknown".

<br>

# mSigAct.server 0.0.4
## Added
* Added signature presence test for **SBS96** artifact and rare signatures.

* Added new internal data containing the signatures aetiology information.

* Added data table to show the thumbnail of signatures and aetiology information.

* Added thumbnail pictures for SBS96, SBS192, DBS78 and ID.

<br>

# mSigAct.server 0.0.3
## Added
* Added functionality to do signature attribution.

<br>

# mSigAct.server 0.0.2
## Added
* Added new *exported* functions for getting signature assignment
`GetExposureWithConfidence` and `GetExposureAndPlotToPdf`.

* Added new *exported* function for getting the subset of signatures for a
specified cancer type from a specified tumor cohort `CancerTypeToSigSubset`.

<br>

# mSigAct.server 0.0.1
* Initial interface for uploading VCFs and generating and downloading zip archive.
