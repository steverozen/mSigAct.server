context("golem tests")

library(golem)

if (FALSE) {
  test_that("app ui", {
    ui <- app_ui()
    expect_shinytaglist(ui)
  })
}


test_that("app server", {
  server <- app_server
  expect_is(server, "function")
})

unlink("Rplots.pdf")









