import readEqdsk
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.constrained_layout.use'] = True

const = pd.read_csv("dist/constants.csv")
dist0 = pd.read_csv("dist/0dist.csv")
dist1 = pd.read_csv("dist/1dist.csv")
traj = pd.read_csv("traj.csv")
# poink = pd.read_csv("poink.csv")

engn = const.iloc[0]['engn']
ekev = const.iloc[0]['ekev']

numpoints = 80000
skip = 2

eq_dat_dict = readEqdsk.read_eqdsk('eqdsk_t217ms')