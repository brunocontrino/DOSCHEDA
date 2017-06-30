context("Doscheda")
test_that("normalize_data med works", {

    chanNames <- paste0('rep1_C',0:2)

    test <- data.frame(matrix(1:9,ncol = 3), Accession = "P18669", UniquePeps = 1:3)
    colnames(test)[1:3] <- chanNames
    test2<- Doscheda:::normalize_data(dataFrame = test,chans = 3,reps = 1,
                              PD2 = FALSE,channelNames = chanNames,
                              incPDofPD = FALSE,removePeptides = FALSE,
                              dataType = 'LFC',normaliseData = 'median',
                              modelType = 'linear',accessionID = 'Accession',
                              uniquePeptides = 'UniquePeps')[1,1]
  expect_that(test2 , equals(-0.5849625))
})
