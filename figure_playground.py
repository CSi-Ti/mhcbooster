import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

lines = list(open('/home/rw762/project/data/batch/JY_Class1_1M_DDA_60min_Slot1-10_1_541.pin'))
xcorrs = []
for i in range(1, len(lines)):
    splits = lines[i].split('\t')
    xcorrs.append(float(splits[8]))

xcorrs = np.array(xcorrs)
plt.hist(xcorrs, bins=100)
plt.show()

print()