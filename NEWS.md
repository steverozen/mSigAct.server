# ICAMS.shiny 0.0.4.90**
## Added
* Added thumbnail pictures for signature attribution results.

* Enabled process VCFs with unknown variant caller.

* Created internal data `COSMIC.v3.sigs` to select the signatures used for 
attribution based on user input.

## Changed
* Create new temporary directory each time when user generate catalogs and do
signature attribution.

*

## Fixed
* Fixed a bug in generating zip archive when the input region is "unknown".

<br>

# ICAMS.shiny 0.0.4
## Added
* Added signature presence test for **SBS96** artifact and rare signatures.

* Added new internal data containing the signatures aetiology information.

* Added data table to show the thumbnail of signatures and aetiology information.

* Added thumbnail pictures for SBS96, SBS192, DBS78 and ID.

<br>

# ICAMS.shiny 0.0.3
## Added
* Added functionality to do signature attribution.

<br>

# ICAMS.shiny 0.0.2
## Added
* Added new *exported* functions for getting signature assignment
`GetExposureWithConfidence` and `GetExposureAndPlotToPdf`.

* Added new *exported* function for getting the subset of signatures for a
specified cancer type from a specified tumor cohort `CancerTypeToSigSubset`.

<br>

# ICAMS.shiny 0.0.1
* Initial interface for uploading VCFs and generating and downloading zip archive.
