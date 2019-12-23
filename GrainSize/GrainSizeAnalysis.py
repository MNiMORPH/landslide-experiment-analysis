import numpy as np
from scipy import stats
from scipy import special
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

data = np.genfromtxt('GrainSize', delimiter='\t', skip_header=1)

def cdfGaussian(x, mean, SD):
    return 0.5 * ( 1 + special.erf( (x-mean)/(SD * 2**.5) ) )

popt, pcov = curve_fit(cdfGaussian, data[:,1], data[:,0]/100., p0=[0.15,0.15])

print 'Mean', popt[0]
print 'SD', popt[1]

plt.plot(data[:,1], data[:,0]/100., 'ko', color='gray')
x = np.linspace(data[:,1][0], data[:,1][-1], 100)
plt.plot(x, cdfGaussian(x, popt[0], popt[1]), 'k-')
plt.xlabel('Grain size [mm]', fontsize=18)
plt.ylabel('Cumulative frequency', fontsize=18)
plt.tight_layout()
plt.savefig('GrainSizeDistribution_AGS_100-140.pdf')
plt.savefig('GrainSizeDistribution_AGS_100-140.png')
plt.savefig('GrainSizeDistribution_AGS_100-140.svg')
plt.show()
