app <- ShinyDriver$new("../")
app$snapshotInit("mutect.test")

app$setInputs(`vcf.type` = "mutect")
app$setInputs(`zipfile.name` = "mutect")
app$setInputs(directory = "click")
app$setInputs(submit = "click")
app$snapshotDownload("download.zipfile")
app$snapshot()
