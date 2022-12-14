library(Rmpfr)
# library(pmcmcDiagnostics)
source("./tests-own/99_helper_diagnostics_model_comparison.R")
load("~/Dropbox/research/usa_energy/04-results/empirical/09_CA.RData")
#
#
#
#
#
# TESTING GENERAL COMPUTATIONS:
# unvec_p_dic <- dic(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9990)
# unvec_p_dic <- dic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9990)
# test_val2 <- test_fun1(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9950)
#
#
#
#
#
# TESTING DIC COMPUTATIONS:
# test_val1 <- dic(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
# test_val2 <- dic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
# print(test_val1, digits = 22)
# print(test_val2, digits = 22)
# print(all.equal(test_val1, test_val2))
# my_bench1 <- microbenchmark::microbenchmark(dic(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9900),
#                                            dic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9900))
# print(my_bench1)
#
#
#
#
#
# TESTING WAIC COMPUTATIONS:
# test_val1 <- waic(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
# test_val2 <- waic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
# print(test_val1, digits = 22)
# print(test_val2, digits = 22)
# print(all.equal(test_val1, test_val2))
# my_bench2 <- microbenchmark::microbenchmark(waic(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9995),
#                                            waic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9995))
# print(my_bench2)
#
#
#
#
#
# TESTING OVERALL FUNCTION FOR LPPD, DIC, WAIC:
ALL <- lppd_dic_waic_test(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
lppd <- test_fun1(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
dic <- dic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)
waic <- waic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9998)

all.equal(ALL[[1]], lppd)
all.equal(ALL[[2]], dic)
all.equal(ALL[[3]], waic)
print(rbind(ALL[[1]], lppd), digits = 22)
print(rbind(ALL[[2]], dic), digits = 22)
print(rbind(ALL[[3]], waic), digits = 22)
#
# my_test_benchmark <- function(burn){
#   a <- dic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, burn);
#   b <- waic_cpp(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, burn);
#   return(list(a, b))
# }
# bench_me3 <- microbenchmark::microbenchmark(lppd_dic_waic_test(y = y_t, X_list = out_pgas_CA_09$xtraj, num_counts = num_counts, 9990),
#                                            my_test_benchmark(9990))
# print(bench_me3)
