import pandas as pd
import matplotlib.pyplot as plt

df1 = pd.read_csv('compare_data_with_out_eNRTL.csv')
df2 = pd.read_csv('compare_data_with_eNRTL.csv')
T_range = df2['temperature'].unique()
mfc2 = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
columns = ['x_CO2', 'Henry_CO2', 'act_coeff_CO2', 'fug_CO2', 'Pressure']

for c in columns:
    fig, axs = plt.subplots(figsize=(12, 10))
    for i, T in enumerate(T_range):
        df_cut_1 = df1[df1['temperature'] == T]
        df_cut_2 = df2[df2['temperature'] == T]
        if c == 'Henry_CO2':
            scale = 1e3
        else:
            scale = 1
        axs.plot(df_cut_1['loading'].to_numpy(), df_cut_1[c].to_numpy()/scale,
                 'x', label=f'T = {T}: Without eNRTL', color=mfc2[i])
        axs.plot(df_cut_2['loading'].to_numpy(), df_cut_2[c].to_numpy()/scale,
                 label=f'T = {T}: With eNRTL', color=mfc2[i])
        axs.legend()
        axs.set_xlabel("CO$_{2}$ Loading, mol CO$_{2}$/mol MEA")
        axs.set_ylabel(c)
        if c != 'act_coeff_CO2':
            axs.set_yscale('log')

fig, axs = plt.subplots(figsize=(12, 10))
for i, T in enumerate(T_range):
    df_cut_1 = df1[df1['temperature'] == T]
    axs.plot(df_cut_1['loading'].to_numpy(), df_cut_1['Henry_CO2'].to_numpy() / 1e3, '--',
             label=f'T = {T}: Without eNRTL: Henry', color=mfc2[i])
    axs.plot(df_cut_1['loading'].to_numpy(), df_cut_1['fug_CO2'].to_numpy(),
             label=f'T = {T}: Without eNRTL: fug_CO2', color=mfc2[i])
    axs.legend()
    axs.set_xlabel("CO$_{2}$ Loading, mol CO$_{2}$/mol MEA")
    axs.set_ylabel("Pressure (kPa)")
    axs.set_yscale('log')

fig, axs = plt.subplots(figsize=(12, 10))
for i, T in enumerate(T_range):
    df_cut_2 = df2[df2['temperature'] == T]
    axs.plot(df_cut_2['loading'].to_numpy(), df_cut_2['Henry_CO2'].to_numpy() / 1e3, '--',
             label=f'T = {T}: With eNRTL: Henry', color=mfc2[i])
    axs.plot(df_cut_2['loading'].to_numpy(), df_cut_2['fug_CO2'].to_numpy(),
             label=f'T = {T}: With eNRTL: fug_CO2', color=mfc2[i])
    axs.legend()
    axs.set_xlabel("CO$_{2}$ Loading, mol CO$_{2}$/mol MEA")
    axs.set_ylabel("Pressure (kPa)")
    axs.set_yscale('log')
plt.show()
