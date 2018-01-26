# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

# Get filename
filename = sys.argv[1]

# Read file into pandas df
df = pd.read_csv(filename, header=0)

split_idx = sum(df.iloc[:,0].notnull())

# Split df into known for fitting and unknown for prediction
conc = df.iloc[:split_idx,0].values
known = df.iloc[:split_idx,1:]
to_predict = df.iloc[split_idx:,1:]

pred_df = pd.DataFrame()

# Regress each column in known df on concentration
# X variable (predictor) is concentration
# Y variable (response) is compound ratio
for compound in known.columns:
    # Fit regression model
    X = conc
    Y = known[compound].values
    
    weights = 1 / X
    # Use minimum ratio as cutoff point for predictions
    min_ratio = min(Y)
    
    # Calculate weighted averages
    X_avg = np.sum(X * weights) / np.sum(weights)
    Y_avg = np.sum(Y * weights) / np.sum(weights)
    
    m = np.sum((Y - Y_avg) * (X - X_avg) * weights) / np.sum(((X - X_avg) ** 2) * weights)
    b = Y_avg - X_avg * m
            
    # Predict concentration based on known area ratio
    conc_pred = (to_predict[compound] - b) / m
    # Set prediction to 0 when ratio is below cutoff
    below_cutoff_mask = to_predict[compound][to_predict[compound] < min_ratio].index.values
    conc_pred[below_cutoff_mask] = 0
    pred_df[compound] = conc_pred

    # Calculate r^2
    try:
    	r2 = r2_score(Y, m * X + b)
    except:
    	print compound, m, b
    	print X
    	print Y

    # Plot fit
    plt.title(compound + ', $R^2$: {}'.format(round(r2, 5)) + '\n y = {}x + {}'.format(round(m, 5), round(b, 5)))
    plt.scatter(X, Y)
    plt.xlabel(u'Concentration (µM)')
    plt.ylabel('Area Ratio')
    abline(m, b)
    plt.savefig("{}.png".format(compound))
    plt.close()

    
# Write results to csv
pred_df.to_csv('results.csv', index=False)