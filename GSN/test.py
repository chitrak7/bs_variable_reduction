import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

x = np.linspace(-10,10,num=1000)
y1 = np.linspace(-10,10,num=1000)
y2 = np.linspace(-10,10,num=1000)
y3 = np.linspace(-10,10,num=1000)
y4 = np.linspace(-10,10,num=1000)
a = np.arange(1,101)
p = 0.5
p1 = 0.8
p2 = 0.3
s1 = 1.264
s2 = 0.774
pow_p = np.power(1-p, a-1)
pow_p1 = np.power(1-p1, a-1)
pow_p2 = np.power(1-p2, a-1)

for i in range(len(x)):
    print(i)
    y1[i] = norm.pdf(x[i],0,1/np.sqrt(p))
    y2[i] = np.sum([norm.pdf(x[i],0,1*np.sqrt(a[j]))*p*pow_p[j] for j in range(100)])
    y3[i] = np.sum([norm.pdf(x[i],0,s1*np.sqrt(a[j]))*p1*pow_p1[j] for j in range(100)])
    y4[i] = np.sum([norm.pdf(x[i],0,s2*np.sqrt(a[j]))*p2*pow_p2[j] for j in range(100)])


plt.plot(x,y1,'r')
plt.plot(x,y2,'g')
plt.plot(x,y3,'b')
plt.plot(x,y4,'y')
plt.show()
