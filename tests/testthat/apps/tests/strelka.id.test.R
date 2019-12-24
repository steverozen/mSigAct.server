app <- ShinyDriver$new("../")
app$snapshotInit("strelka.id.test")

app$setInputs(`vcf.type` = "strelka.id")
app$setInputs(`zipfile.name` = "strelka.id")
app$setInputs(directory = "click")
app$setInputs(submit = "click")
app$snapshotDownload("download.zipfile")
app$snapshot()
