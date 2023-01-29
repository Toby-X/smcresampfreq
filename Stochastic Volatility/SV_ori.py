# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.stats import norm
from tqdm import tqdm

# INITIAL VALUES
t = 50  # Dimension
n = 500  # number of particles
n_p = int(1e6)  # to estimate posterior mean
m = int(200)  # parallel experiemnts
p = 100  # of different series
x = np.zeros(t)
y = np.zeros(t)

alpha = .91
sigma = 1.0
beta = .5

w_u = np.ones(n)/n
threshold = np.linspace(0.1, 1, 10)
x = np.zeros(n*t).reshape([n, t])
x0 = np.zeros(int(n_p*t)).reshape([int(n_p), t])
x_estimate_ss = np.zeros(int(t*m*len(threshold))
                         ).reshape([int(m), int(t), int(len(threshold))])
x_estimate_ori = np.zeros(t)
xo = np.zeros(t)
yo = np.zeros(t)
mse_all = np.zeros(p*len(threshold)).reshape([p, int(len(threshold))])

for _ in tqdm(range(p)):
    np.random.seed(_+1)

    # Generate Samples
    v = np.random.normal(size=t)
    u = np.random.normal(size=t)
    xo[0] = np.random.normal(0, sigma**2/(1-alpha**2), 1)
    yo[0] = beta*np.exp(xo[0]/2)*u[0]
    for i in range(1, t):
        xo[i] = alpha*xo[i-1]+sigma*v[i]
        yo[i] = beta*np.exp(xo[i]/2)*u[i]

    # Generate Estimate Mean using bootstrap
    # x0[:,0]=np.random.normal(size=n_p)
    # w_u0 = norm.pdf(x0[:,0],0,sigma**2/(1-alpha**2))*norm.pdf(yo[0],0,beta*np.exp(x0[:,0]/2))/norm.pdf(x0[:,0])
    # w_n0 = np.divide(w_u0,np.sum(w_u0))
    # x_estimate_ori[0] = sum(w_n0*x0[:,0])
    # idx = np.random.choice(range(n_p),size=n_p,replace=True,p=w_n0)
    # x0 = x0[idx,:]
    # w_u0 = np.ones(n_p)/n_p

    # for i in range(1,t):
    #     x0[:,i]=np.random.normal(x0[:,i-1],sigma**2,n_p)
    #     w_u0 = norm.pdf(x0[:,i],loc=alpha*x0[:,i-1],scale=sigma**2)*norm.pdf(yo[i],loc=0,scale=beta*np.exp(x0[:,i]/2))/norm.pdf(x0[:,i],x0[:,i-1]+yo[i])
    #     w_n0 = np.divide(w_u0,np.sum(w_u0))
    #     x_estimate_ori[i] = np.sum(w_n0*x0[:,i])
    #     idx = np.random.choice(range(n_p),n_p,True,w_n0)
    #     x0 = x0[idx,:]

    # Generate estimates
    for l in range(int(m)):
        for j in range(int(len(threshold))):
            x[:, 0] = np.random.normal(size=n)
            w_u = norm.pdf(x[:, 0], 0, sigma**2/(1-alpha**2)) * \
                norm.pdf(yo[0], 0, beta*np.exp(x[:, 0]/2))/norm.pdf(x[:, 0])
            w_n = w_u/np.sum(w_u)
            x_estimate_ss[l, 0, j] = np.sum(w_n*x[:, 0])
            if 1/np.sum(w_n**2) < threshold[j]*n:
                idx = np.random.choice(range(n), n, True, w_n)
                x = x[idx, :]
                w_u = np.ones(n)/n
            for i in range(1, t):
                x[:, i] = np.random.normal(x[:, i-1]+yo[i], size=n)
                lw_tmp = np.log(norm.pdf(x[:, i], loc=alpha*x[:, i-1], scale=sigma**2))+np.log(norm.pdf(
                    yo[i], loc=0, scale=beta*np.exp(x[:, i]/2)))+(-1)*np.log(norm.pdf(x[:, i], x[:, i-1]+yo[i]))
                if np.all(lw_tmp < np.log(1e-100)):
                    lw_tmp += np.log(1e200)
                w_u = np.exp(np.log(w_u)+lw_tmp)
                w_n = w_u/np.sum(w_u)
                # if np.all(w_u) == 0:
                #     x_estimate_ss[l,i,j]=9999
                #     continue
                x_estimate_ss[l, i, j] = np.sum(w_n*x[:, i])
                if 1/np.sum(w_n**2) < threshold[j]*n:
                    idx = np.random.choice(range(n), n, True, w_n)
                    x = x[idx, :]
                    w_u = np.ones(n)/n

    # Estimate
    mse_ss = np.zeros(int(len(threshold)*t)).reshape([int(len(threshold)), t])
    for k in range(int(len(threshold))):
        for i in range(int(m)):
            mse_ss[k, :] += (xo-x_estimate_ss[i, :, k])**2
    mse_ss = mse_ss/m
    mse_sum_ss = np.sum(mse_ss, axis=1)
    mse_all[_, :] = mse_sum_ss

mse_all_df = pd.DataFrame(mse_all)
mse_all_df.to_csv("mse_all_alpha_08.csv")
