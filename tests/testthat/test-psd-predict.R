# 
# Make sure we honor the choice of xmin when predicting values after fitting distributions
# 

test_that("xmin is used when predicting values", { 
  
  # Test for single-element run of patchdistr_sews
  a <- patchdistr_sews(serengeti[[12]], xmin = 1)
  preds <- predict(a)[["pred"]]
  expect_true({ 
    abs(max(preds[ ,"y"], na.rm = TRUE) - 1) < 1e-8
  })
  
  # Test for single-element run of patchdistr_sews, now with xmin
  a <- patchdistr_sews(serengeti[[12]], xmin = 10)
  preds <- predict(a, xmin_rescale = TRUE)[["pred"]]
  # should be rescaled so that different from 1
  expect_true({ 
    abs(max(preds[ ,"y"], na.rm = TRUE) - 1) > 0.1 
  })
  
  # Test for multiple-matrix run of patchdistr_sews
  a <- patchdistr_sews(serengeti[11:12], xmin = 1)
  preds <- predict(a)[["pred"]]
  expect_true({
    abs(max(preds[ ,"y"], na.rm = TRUE) - 1) < 1e-8
  })
  
  # Test for multiple-matrix run of patchdistr_sews
  a <- patchdistr_sews(serengeti[11:12], xmin = 10)
  preds <- predict(a, xmin_rescale = TRUE)[["pred"]]
  expect_true({
    abs(max(preds[ ,"y"], na.rm = TRUE) - 1) > 0.1
  })
  
})
