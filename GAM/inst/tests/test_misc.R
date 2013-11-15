context("Misc")

test_that("getIdsType works", {
    id.map <- data.frame(ID1=c("a", "b", "c"), ID2=c("a", "a1", "a2"))
    ids <- c("a2", "a")
    expect_equal(getIdsType(ids, id.map), "ID2")
})