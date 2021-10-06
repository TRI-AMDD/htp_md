import htpmd
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(precision=6)

results = htpmd.analyze('test_data/test-arash')

print('Time avg. MSD', results['li_msd_curve'])
print('Non avg. MSD', results['li_msd_curve_non_avg'])
print(results.keys())
ts, li_msd = results['li_msd_curve']
ts_non_avg, li_msd_non_avg = results['li_msd_curve_non_avg']

print(ts.shape)
print(ts_non_avg.shape)

def make_plot(x, y, label):
    plt.title(label, fontdict=None, loc='center')
    plt.plot(x, y, marker = 'o')
    plt.show()

make_plot(ts, li_msd, label = 'Time Averaged MSD')
make_plot(ts_non_avg, li_msd_non_avg, label = 'Non-Averaged MSD')
