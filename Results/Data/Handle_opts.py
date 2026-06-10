import pandas as pd
import matplotlib.pyplot as plt

folder_h = "Hydrogen/"

optimos_30d_rq03_h = pd.read_csv(folder_h+'Hydrogen_opti_30d_LFs_rq03_opts.csv')
optimos_30d_rq07_h = pd.read_csv(folder_h+'Hydrogen_opti_30d_LFs_rq07_opts.csv')

folder_m = "Methane/"
optimos_30d_rq03_m = pd.read_csv(folder_m+'Methane_opti_30d_LFs_rq03_opts.csv')
optimos_30d_rq07_m = pd.read_csv(folder_m+'Methane_opti_30d_LFs_rq07_opts.csv')
#LAES
optimos_30d_rq03_l = pd.read_csv('LAES_opti_12h_LFs_rq03_opts.csv')
optimos_30d_rq07_l = pd.read_csv('LAES_opti_12h_LFs_rq07_opts.csv')

dfs = []
optimos_30d_rq03_h["Cryogen"] = "LH2_rq03"
dfs.append(optimos_30d_rq03_h)
optimos_30d_rq07_h["Cryogen"] = "LH2_rq07"
dfs.append(optimos_30d_rq07_h)
optimos_30d_rq03_m["Cryogen"] = "LNG_rq03"
dfs.append(optimos_30d_rq03_m)
optimos_30d_rq07_m["Cryogen"] = "LNG_rq07"
dfs.append(optimos_30d_rq07_m)
optimos_30d_rq03_l["Cryogen"] = "LAES_rq03"
dfs.append(optimos_30d_rq03_l)
optimos_30d_rq07_l["Cryogen"] = "LAES_rq07"
dfs.append(optimos_30d_rq07_l)

# Concatenate into one DataFrame
combined = pd.concat(dfs, ignore_index=True)

# Save to a single CSV
combined.to_csv("all_opts.csv", index=False)

# Example: extract only Exp3
exp3_data = combined[combined["Cryogen"] == "LAES_rq03"]