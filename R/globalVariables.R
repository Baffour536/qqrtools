# Register non-standard-evaluation variables to satisfy R CMD check:
# ".w" is the weights column used inside quantreg::rq() in .qqr_local_fit();
# "tau", "value", "method" are ggplot2 aesthetics in plot.compare.qqr().
utils::globalVariables(c(".w", "tau", "value", "method"))
