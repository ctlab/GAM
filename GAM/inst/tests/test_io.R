context("IO")

data("examplesGAM")

test_that("saveModuleToDot works", {
    file <- tempfile(fileext = ".dot")
    saveModuleToDot(module.re, file, name = "module.re")
    
    saveModuleToDot(module.rn, file, name = "module.rn")
})
