import pandas as pd
import matplotlib.pyplot as plt

fname='monthly_relh.csv'

df1 = pd.read_csv(fname, index_col=['month'])
print(df1)
print(df1.columns.values)
print(df1.mean(axis=1))
plt.plot(df1.mean(axis=1), '-c.')

fname='weekly_relh.csv'
df = pd.read_csv(fname, index_col=['week'])
print(df.mean(axis=1))

plt.plot(df.mean(axis=1), '-b.')
plt.show()
