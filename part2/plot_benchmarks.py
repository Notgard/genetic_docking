import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("scaling_results.csv")

# Compute speedups relative to 1-thread run
base_cpu = df.loc[df['THREADS'] == 1, 'CPU_TIME'].values[0]
base_omp = df.loc[df['THREADS'] == 1, 'OMP_WALLTIME'].values[0]

df['CPU_SPEEDUP'] = base_cpu / df['CPU_TIME']
df['OMP_SPEEDUP'] = base_omp / df['OMP_WALLTIME']

# Plot time vs threads
plt.figure(figsize=(10, 5))
plt.plot(df['THREADS'], df['CPU_TIME'], marker='o', label='CPU Time')
plt.plot(df['THREADS'], df['OMP_WALLTIME'], marker='o', label='OMP Walltime')
plt.xlabel('Number of Threads')
plt.ylabel('Time (s)')
plt.title('Execution Time vs Threads')
plt.grid(True)
plt.legend()
plt.xscale('log', base=2)
plt.xticks(df['THREADS'])
plt.savefig("time_scaling.png", dpi=300)

# Plot speedup vs threads
plt.figure(figsize=(10, 5))
plt.plot(df['THREADS'], df['CPU_SPEEDUP'], marker='o', label='CPU Speedup')
plt.plot(df['THREADS'], df['OMP_SPEEDUP'], marker='o', label='OMP Walltime Speedup')
plt.xlabel('Number of Threads')
plt.ylabel('Speedup')
plt.title('Speedup vs Threads')
plt.grid(True)
plt.legend()
plt.xscale('log', base=2)
plt.xticks(df['THREADS'])
plt.savefig("speedup_scaling.png", dpi=300)

print("Plots saved: time_scaling.png and speedup_scaling.png")
