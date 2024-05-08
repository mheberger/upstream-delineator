"""
Plot the distribution of areas of the unit catchments.

"""
import psycopg2
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from util.db import cursor, db_close

id = 77

sql1 = f"SELECT unitarea from merit_basins_{id} WHERE unitarea IS NOT NULL and numup is NOT NULL;"
sql2 = f"SELECT unitarea from merit_basins2_{id} WHERE unitarea IS NOT NULL;"
sql4 = f"SELECT unitarea from merit_basins4_{id} WHERE unitarea IS NOT NULL;"
sql5 = f"SELECT unitarea from merit_basins5_{id} WHERE unitarea IS NOT NULL;"
#sql7 = f"SELECT unitarea from merit_basins7_{id} WHERE unitarea IS NOT NULL;"


cursor.execute(sql1)
rows1 = cursor.fetchall()

cursor.execute(sql2)
rows2 = cursor.fetchall()

cursor.execute(sql4)
rows4 = cursor.fetchall()

cursor.execute(sql5)
rows5 = cursor.fetchall()

#cursor.execute(sql7)
#rows7 = cursor.fetchall()

# Print the number of features
print("Number of unit cathcments:")
print(f"Step 1: {len(rows1)}")
print(f"Step 2: {len(rows2)}")
print(f"Step 4: {len(rows4)}")
print(f"Step 6: {len(rows5)}")
#print(f"Step 7: {len(rows7)}")

# Extract the values from the result
values1 = [row[0] for row in rows1]
values2 = [row[0] for row in rows2]
values4 = [row[0] for row in rows4]
values5 = [row[0] for row in rows5]
#values7 = [row[0] for row in rows7]

db_close()

# Compute the cumulative distribution function (CDF)
sorted_values1 = np.sort(values1)
sorted_values2 = np.sort(values2)
sorted_values4 = np.sort(values4)
sorted_values5 = np.sort(values5)
#sorted_values7 = np.sort(values7)

cumulative_prob1 = np.linspace(0, 1, len(sorted_values1))
cumulative_prob2 = np.linspace(0, 1, len(sorted_values2))
cumulative_prob4 = np.linspace(0, 1, len(sorted_values4))
cumulative_prob5 = np.linspace(0, 1, len(sorted_values5))
#cumulative_prob7 = np.linspace(0, 1, len(sorted_values7))

# Plot the CDF
plt.plot(sorted_values1, cumulative_prob1, label="Unit Catchments", marker='.', linestyle='none')
plt.plot(sorted_values2, cumulative_prob2, label="Merged2", marker='o', linestyle='none')
plt.plot(sorted_values4, cumulative_prob4, label="Merged4", marker='.', linestyle='none')
plt.plot(sorted_values5, cumulative_prob5, label="Merged6", marker='.', linestyle='none')
#plt.plot(sorted_values7, cumulative_prob7, label="Merged7", marker='.', linestyle='none')

plt.xlabel('Unit catchment area, km²')
plt.ylabel('Cumulative Probability')
plt.title('Cumulative Distribution Function (CDF) plot of unit catchment areas')
plt.grid(True)
plt.legend()
plt.show()

# Try the kernel density plot
LOG = False
sns.kdeplot(values1, label='Unit catchments', color='blue',   linestyle='-',  log_scale=LOG)
sns.kdeplot(values2, label='Merged2',         color='red',    linestyle='--', log_scale=LOG)
sns.kdeplot(values4, label='Merged4',         color='red',    linestyle='--', log_scale=LOG)
sns.kdeplot(values5, label='Merged5',         color='green',  linestyle='-',  log_scale=LOG)
#sns.kdeplot(values7, label='Merged7',         color='purple', linestyle='-',  log_scale=LOG)

sns.set_style('whitegrid')  # Optional: Set the style
plt.xlabel('Catchment area, km²')
plt.ylabel('Density')
plt.title('Kernel Density Plot')
plt.legend()
plt.show()


