# Wishlist and TDL

* Add timing/cost
* set memory for tasks from command line option
* rename targetTR to hipstr? it can actually be used to perform non-targeted analysis
* how to extend to other tools (e.g. EH)

* where are the reads in the new release??
* check final batch error
* use new wdl arguments
* `GetFileBatches` function is shared by AoU/UKB. could put in single library to avoid copying code
* Related ambitious goal: could make a launcher super class, and have AoU UKB subclasses that extend that. the launchers basically have same functions but different implementations now (`RunWorkflow`, `RunWorkflow`/`UploadGS`). Some of those are even more general than just targetTR. or at least move commonly used things to a separate module