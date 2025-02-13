#!  .TSenv/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 03:22:28 2025

@author: savvaspitsi
"""
from datetime import datetime 
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.tsa.stattools import adfuller
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf



def to_time(s):
    return datetime.strptime(s[:-6],"%Y-%m-%d %H:%M:%S")


def get_returns(data:  list, ticker: str):
    returns = np.zeros(len(data)-1)
    for i in range(len(data)-1):
        returns[i] = np.log(data[ticker].iloc[i+1]/data[ticker].iloc[i])
    return returns


def norm_test(data: list, pvalue = 0.05, bins=20):
    p = stats.shapiro(data).pvalue
    p2 = stats.ks_1samp(data, stats.norm.cdf).pvalue
    if p > pvalue:
        print(f"(Shapiro) Fail to reject H0: Data appears to be normally distributed. p-value = {p:.4f}")
    else:
        print(f"(Shapiro) Reject H0: Data does not appear to be normally distributed. p-value = {p:.4f}")

    if p2 > pvalue:
        print(f"(Kolmogorov-Smirnov) Fail to reject H0: Data appears to be normally distributed. p-value = {p2:.4f}")
    else:
        print(f"(Kolmogorov-Smirnov) Reject H0: Data does not appear to be normally distributed. p-value = {p2:.4f}")
    
    fig, (ax1, ax2) = plt.subplots(1,2, squeeze=True)
    stats.probplot(data, dist=stats.norm, plot=ax2)
    ax1.hist(data, bins = bins)
    plt.tight_layout()
    plt.show()



def ADF_test(data: list, pvalue = 0.05):
    
    ADFresult = adfuller(data)
    
    if ADFresult[1] > pvalue:
        print(f"(Aug Dickey-Fuller) Fail to reject H0: Series is non-stationary, p-value = {ADFresult[1]:.4f}")
    else:
        print(f"(Aug Dickey-Fuller) Reject H0: Series is stationary, p-value = {ADFresult[1]:.4f}")


def plot_autocorr(data: list):
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))
    plot_acf(data, ax=ax1)
    plot_pacf(data, ax=ax2)
    plt.show()