import numpy as np
import matplotlib.pyplot as plt
import pyreadr
from scipy.stats import norm, gaussian_kde
import seaborn as sns
import matplotlib.lines as mlines
mu=5
g_p=0.95
sigma=1
alpha_s=1
n=100
fname = "symGSN_plot_mu_em_glmnet_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f.rds" %(mu,g_p,sigma,alpha_s,n)
sname = "symGSN_plot_mu_em_glmnet_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f.png" %(mu,g_p,sigma,alpha_s,n)
print(fname)
res = pyreadr.read_r(fname)
#print(res[None])
res = np.array(res[None])
mu_hat = res[0,0]
sigma_hat=res[0,1]
p_hat=res[0,2]
res1 = np.array(res[:,5:105])
res2 = np.array(res[:,105:205])

def gsnpts(mu, p, sigma, x):
    a = np.arange(1,101)
    pow_p = np.power(1-p, a-1)
    y = np.zeros(shape=(len(x)))
    for i in range(len(x)):
        x1 = x[i]
        y[i] = np.sum([norm.pdf(x1,a[j]*mu,sigma*np.sqrt(a[j]))*p*pow_p[j] for j in range(100)])
    return y

x = np.linspace(-3,28,num=400)
ye = gsnpts(mu_hat, p_hat, sigma_hat, x)
ya = gsnpts(mu, g_p, 1, x)
plt.figure(figsize=(10,20))
ax1 = plt.subplot2grid((9,5),(0,1),3,3)
ax2 = plt.subplot2grid((9,5),(4,0),2,2)
ax3 = plt.subplot2grid((9,5),(4,3),2,2)
ax4 = plt.subplot2grid((9,5),(7,0),2,2)
ax5 = plt.subplot2grid((9,5),(7,3),2,2)


ax1.plot(x,ya,color="red",label="actual")
ax1.plot(x,ye,color="blue",label="estimated")
ax1.set_title("plot of pdf")
ax1.legend()

red_patch = mlines.Line2D([], [], color='r', label='actual')
blue_patch = mlines.Line2D([], [], color='b', label='estimated')
hdl = [red_patch, blue_patch]
sns.kdeplot(res1[0,:], bw=.25,color="r",ax=ax2)
sns.kdeplot(res2[0,:], bw=.25,color="b",ax=ax2)

ax2.legend(handles=hdl)
sns.kdeplot(res1[1,:], bw=.25,color="r",ax=ax3)
sns.kdeplot(res2[1,:], bw=.25,color="b",ax=ax3)
ax3.legend(handles=hdl)
sns.kdeplot(res1[2,:], bw=.25,color="r",ax=ax4)
sns.kdeplot(res2[2,:], bw=.25,color="b",ax=ax4)
ax4.legend(handles=hdl)
sns.kdeplot(res1[3,:], bw=.25,color="r",ax=ax5)
sns.kdeplot(res2[3,:], bw=.25,color="b",ax=ax5)
ax5.legend(handles=hdl)
#plt.show()
plt.savefig(sname)




print(p_hat, mu_hat)
print(len(res1[0,:]))
