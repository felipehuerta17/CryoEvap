import numpy as np
import pandas as pd

a_opt_analytical = lambda ru, lf: (ru) / (2 *(lf + (1-lf)*eta))


r_U_list = [0.25, 1, 4]
LF_list  = [0.05, 0.5, 0.95]
eta      = 0.9

opt_results = []
for r_U in r_U_list:
    for LF in LF_list:
        a_opt = a_opt_analytical(r_U, LF)
        opt_results.append({'r_U': r_U, 'LF': LF, 'a_opt': a_opt})

df_opt_results = pd.DataFrame(opt_results)
df_opt_results.to_csv("../Results_new/Data/analytical_optimal_simplified.csv", index=False)