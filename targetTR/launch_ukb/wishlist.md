# Wishlist and TDL

* where are the reads in the new release??
* check final batch error
* use new wdl arguments
* `GetFileBatches` function is shared by AoU/UKB. could put in single library to avoid copying code
* Related ambitious goal: could make a launcher super class, and have AoU UKB subclasses that extend that. the launchers basically have same functions but different implementations now (`RunWorkflow`, `RunWorkflow`/`UploadGS`). Some of those are even more general than just targetTR