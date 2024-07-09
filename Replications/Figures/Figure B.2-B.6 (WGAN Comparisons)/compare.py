####################
# Comparisons between the Keele & Titiunik 2015
# observed samples and the WGAN-generated samples
####################

import wgan
import pandas as pd
from copy import copy
import matplotlib.pyplot as plt
import torch
import numpy as np

torch.manual_seed(0)
np.random.seed(0)

# The two functions below are adapted from Athey et al. (2021)
def save_corr(df_real, df_fake, figsize=5, path=""):
    if "source" in list(df_real.columns): df_real = df_real.drop("source", axis=1)
    if "source" in list(df_fake.columns): df_fake = df_fake.drop("source", axis=1)
    df_real.insert(0, "source", "real"), df_fake.insert(0, "source", "fake")
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    df_joined = pd.concat([df_real[common_cols], df_fake[common_cols]], axis=0, ignore_index=True)
    df_real, df_fake = df_real.drop("source", axis=1), df_fake.drop("source", axis=1)
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    fig1 = plt.figure(figsize=(figsize * 2, figsize * 1))
    s1 = [fig1.add_subplot(1, 2, i) for i in range(1, 3)]
    s1[0].set_xlabel("real")
    s1[1].set_xlabel("fake")
    s1[0].matshow(df_real[common_cols].corr())
    s1[1].matshow(df_fake[common_cols].corr())
    fig1.savefig(path+'_corr.png')

def compare_dfs(df_real, df_fake, scatterplot=dict(x=[], y=[], samples=400, smooth=0),
                table_groupby=[], histogram=dict(variables=[], nrow=1, ncol=1),
                figsize=3,save=False,path=""):
    """
    Diagnostic function for comparing real and generated data from WGAN models.
    Prints out comparison of means, comparisons of standard deviations, and histograms
    and scatterplots.

    Parameters
    ----------
    df_real: pandas.DataFrame
        real data
    df_fake: pandas.DataFrame
        data produced by generator
    scatterplot: dict
        Contains specifications for plotting scatterplots of variables in real and fake data
    table_groupby: list
        List of variables to group mean and standard deviation table by
    histogram: dict
        Contains specifications for plotting histograms comparing marginal densities
        of real and fake data
    save: bool
        Indicate whether to save results to file or print them
    path: string
        Path to save diagnostics for model
    """
    # data prep
    if "source" in list(df_real.columns): df_real = df_real.drop("source", axis=1)
    if "source" in list(df_fake.columns): df_fake = df_fake.drop("source", axis=1)
    df_real.insert(0, "source", "real"), df_fake.insert(0, "source", "fake")
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    df_joined = pd.concat([df_real[common_cols], df_fake[common_cols]], axis=0, ignore_index=True)
    df_real, df_fake = df_real.drop("source", axis=1), df_fake.drop("source", axis=1)
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    # mean and std table

    means = df_joined.groupby(table_groupby + ["source"]).mean().round(2).transpose()
    if save:
        means.to_csv(path+"_means.txt",sep=" ")
    else:
        print("-------------comparison of means-------------")
        print(means)

    stds = df_joined.groupby(table_groupby + ["source"]).std().round(2).transpose()

    if save:
        stds.to_csv(path+"_stds.txt",sep=" ")
    else:
        print("-------------comparison of stds-------------")
        print(stds)
    # covariance matrix comparison
    fig1 = plt.figure(figsize=(figsize * 2, figsize * 1))
    s1 = [fig1.add_subplot(1, 2, i) for i in range(1, 3)]
    s1[0].set_xlabel("real")
    s1[1].set_xlabel("fake")
    s1[0].matshow(df_real[common_cols].corr())
    s1[1].matshow(df_fake[common_cols].corr())
    # histogram marginals
    if histogram and len(histogram["variables"]) > 0:
        fig2, axarr2 = plt.subplots(histogram["nrow"], histogram["ncol"],
                                    figsize=(histogram["nrow"]*figsize, histogram["ncol"]*figsize))
        v = 0
        for i in range(histogram["nrow"]):
            for j in range(histogram["ncol"]):
                plot_var, v = histogram["variables"][v], v+1
                axarr2[i][j].hist([df_real[plot_var], df_fake[plot_var]], bins=8, density=1,
                                  histtype='bar', label=["real", "fake"], color=["blue", "red"])
                axarr2[i][j].legend(prop={"size": 10})
                axarr2[i][j].set_title(plot_var)
        if save:
            fig2.savefig(path+'_hist.png')
        else:
            fig2.show()
            
    # scatterplot grid
    if scatterplot and len(scatterplot["x"]) * len(scatterplot["y"]) > 0:
        df_real_sample = df_real.sample(scatterplot["samples"])
        df_fake_sample = df_fake.sample(scatterplot["samples"])
        x_vars, y_vars = scatterplot["x"], scatterplot["y"]
        fig3 = plt.figure(figsize=(len(x_vars) * figsize, len(y_vars) * figsize))
        s3 = [fig3.add_subplot(len(y_vars), len(x_vars), i + 1) for i in range(len(x_vars) * len(y_vars))]
        for y in y_vars:
            for x in x_vars:
                s = s3.pop(0)
                x_real, y_real = df_real_sample[x].to_numpy(),  df_real_sample[y].to_numpy()
                x_fake, y_fake = df_fake_sample[x].to_numpy(), df_fake_sample[y].to_numpy()
                from math import sqrt,pi
                def fit(xx, yy):
                    xx, yy = torch.tensor(xx).to(torch.float), torch.tensor(yy).to(torch.float)
                    xx = (xx - xx.mean())/ xx.std()
                    bw = 1e-9 + scatterplot["smooth"] # * (xx.max()-xx.min())
                    dist = (xx.unsqueeze(0) - xx.unsqueeze(1)).pow(2)/bw
                    kern = 1/sqrt(2*pi)*torch.exp(-dist**2/2)
                    w = kern / kern.sum(1, keepdim=True)
                    y_hat = w.mm(yy.unsqueeze(1)).squeeze()
                    return y_hat.detach().numpy()
                y_real, y_fake = fit(x_real, y_real), fit(x_fake, y_fake)
                s.scatter(x_real, y_real, color="blue")
                s.scatter(x_fake, y_fake, color="red")
                s.set_ylabel(y)
                s.set_xlabel(x)

        if save:
            fig3.savefig(path+'_scatter.png')
        else:
            fig3.show()


df = pd.read_csv('../../Data/Keele & Titiunik 2015/Observed Samples/data_e2008g.csv')
df = df.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Turnout"})
df_gen = pd.read_csv('../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/turnout_cWGAN.csv')
df_gen = df_gen.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Turnout"})

wgan.compare_dfs(df, df_gen, 
                 scatterplot=dict(x=["Latitude", "Longitude"], y=["Turnout"], samples=200, smooth=0),
                 table_groupby=["W"],
                 histogram=dict(variables=["Turnout", "Turnout", "Latitude", "Longitude"], nrow=2, ncol=2),
                 figsize=5, save=True, path="plots/turnout")
save_corr(df, df_gen, figsize=5, path="plots/turnout")

df = pd.read_csv('../../Data/Keele & Titiunik 2015/Observed Samples/data_price.csv')
df = df.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Price"})
df_gen = pd.read_csv('../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/price_cWGAN.csv')
df_gen = df_gen.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Price"})

wgan.compare_dfs(df, df_gen, 
                 scatterplot=dict(x=["Latitude", "Longitude"], y=["Price"], samples=200, smooth=0),
                 table_groupby=["W"],
                 histogram=dict(variables=["Price", "Price", "Latitude", "Longitude"], nrow=2, ncol=2),
                 figsize=5, save=True, path="plots/price")
save_corr(df, df_gen, figsize=5, path="plots/price")

df = pd.read_csv('../../Data/Keele & Titiunik 2015/Observed Samples/data_age.csv')
df = df.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Age"})
df_gen = pd.read_csv('../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/age_cWGAN.csv')
df_gen = df_gen.rename(columns={"X1": "Latitude", "X2": "Longitude", "Y": "Age"})

wgan.compare_dfs(df, df_gen, 
                 scatterplot=dict(x=["Latitude", "Longitude"], y=["Age"], samples=200, smooth=0),
                 table_groupby=["W"],
                 histogram=dict(variables=["Age", "Age", "Latitude", "Longitude"], nrow=2, ncol=2),
                 figsize=5, save=True, path="plots/age")
save_corr(df, df_gen, figsize=5, path="plots/age")
