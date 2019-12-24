app <- ShinyDriver$new("../")
app$snapshotInit("strelka.sbs.test")

app$setInputs(`zipfile.name` = "strelka.sbs")
app$setInputs(directory = "click")
app$setInputs(submit = "click")
app$snapshotDownload("download.zipfile")
