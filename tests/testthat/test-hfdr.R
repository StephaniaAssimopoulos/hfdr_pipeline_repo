test_that("bonf_minp works", {
  source("../../R/hfdr.R")
  expect_equal(bonf_minp(c(0.01, 0.2, 0.3)), min(1, 3*0.01))
})

test_that("simes_p works on simple case", {
  source("../../R/hfdr.R")
  # p sorted: 0.01,0.02,0.5 => min(3/1*0.01=0.03, 3/2*0.02=0.03, 3/3*0.5=0.5)=0.03
  expect_equal(simes_p(c(0.5,0.02,0.01)), 0.03, tolerance=1e-12)
})
