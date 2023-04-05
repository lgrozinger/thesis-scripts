using Pyolin
using Plots

ins = RPUExperiment(Pyolin.search("KT2440", "pSeva221", "1818"))
outs = RPUExperiment(Pyolin.search("KT2440", "pSeva221", "Lmra_n1"))
