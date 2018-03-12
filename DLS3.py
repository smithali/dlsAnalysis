from __future__ import division
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
from pylab import *
from scipy.optimize import curve_fit
from math import exp, log
from ipykernel import kernelapp as app
# from uncertainties import ufloat

#Enter date & filename of measurement
date = '070318'
filename = '10mMTPE_THF_90.csv'

result_path = '/Users/alicesmith/Documents/Data/DLS PTCL/' + date + '/TPE Sols/'

# Read in the data file
rawdata = pd.read_csv('/Users/alicesmith/Documents/Data/DLS PTCL/'+ date + '/TPE Sols/' + filename, skiprows = 27, header=0, nrows=191)
    #print(rawdata)
    #type(rawdata['Lag [ms]'][1])
    #rawdata['Lag [ms]'] = type(rawdata['Lag [ms]'])
    #rawdata['Average'] = np.array(rawdata['Average'])

# print(rawdata['Lag [ms]'][len(rawdata['Average']) - 1 ])
# print(rawdata['Lag [ms]'][0])

plt.figure(1)
plt.semilogx(rawdata['Lag [ms]'], rawdata['Average'])
    #plt.xlabel(str(header[0]))
    #plt.ylabel(str(header[1]))
plt.show()

# Define fit functions. func = function as is, func2 = log'd functions to remove exponentially high values
def func(t, A, b, gamma, d):
    #return (A + b * np.exp( -1 * gamma * t + ( d ** 2) * ( t ** 2) / 2) )
    return A*np.exp(t)

def func2(t,A,b,gamma, d):
    return np.log(abs(A)) + b * (-1 * gamma * t + (d**2 * t**2)/2)

# Assign data to x,y & test func2
x = rawdata['Lag [ms]']
y = np.log(abs(rawdata['Average']))

# test = func2(y, 0.25,0.25,0.01,0.0001)
# print(test)
# weighting? sigma=100/y, absolute_sigma=True

#Fit to func2
popt, pcov = curve_fit(func2, x, y, p0=[0.25,0.25,0.01,0.0001])
print(popt)
print(pcov)

plt.figure(2)
plt.semilogx(x, np.exp(y), '-o')
x = np.linspace((rawdata['Lag [ms]'][0]), (rawdata['Lag [ms]'][len(rawdata['Average']) - 1 ]), 191)
#print(x)
plt.semilogx(x, np.exp(func2(x, *popt)))
plt.show()

plt.figure(3)
x = np.linspace((rawdata['Lag [ms]'][0]), (rawdata['Lag [ms]'][len(rawdata['Average']) - 1 ]), 191)
plt.semilogx(x, np.exp(func2(x, *popt)))
plt.show()

# Get errors from covariance matrix. Extract diagonals then square root to find error of each A,b,gamma, d.


# Get sample data & experimental parameters
rawdata2 = pd.read_csv('/Users/alicesmith/Documents/Data/DLS PTCL/'+ date + '/TPE Sols/' + filename, nrows=24)
print(rawdata2)

name = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Samplename : '].index[0]]
T = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Temperature [K] :'].index[0]]
V = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Viscosity [cp]  :'].index[0]]
n = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Refractive Index:'].index[0]]
wavelength = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Wavelength [nm] :'].index[0]]
angle = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Angle [∞]       :'].index[0]]
duration = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Duration [s]    :'].index[0]]
runs = rawdata2['Unnamed: 1'][rawdata2.loc[rawdata2['ALV-5000/E-WIN Data'] == 'Runs            :'].index[0]]

# Get theta from angle
angle = np.array(angle)
print(type(angle))
angleflt = angle.astype(np.float)
theta = str(angleflt * np.pi / 180)

# Extract as strings to check everything looks good
print('Sample Name = ' + name + '\nTemperature = ' + T + 'K' + '\nViscosity = ' + V + 'cp' + '\nRefractive Index = ' + n + '\nWavelength = ' + wavelength + 'nm' + '\nAngle = ' + str(angle) + 'deg.' + '\nTheta = ' + theta + 'rad' + '\nDuration = ' + duration + 's' + '\nRuns = ' + runs)

# Convert parameters to floats for calculations.
A = popt[0].astype(np.float)
b = popt[1].astype(np.float)
gamma=(popt[2]).astype(np.float)
d = popt[3].astype(np.float)
kB = 1.3806488*10**(-23)
theta = np.array(theta).astype(np.float)
T = np.array(T).astype(np.float)
V = np.array(V).astype(np.float)
n = np.array(n).astype(np.float)
wavelength = np.array(wavelength).astype(np.float)
q = ((4* np.pi * n )/wavelength) * np.sin(theta/2)

def diffusioncoef(gamma,q):
    return (gamma/(q**2))
    ;
def radius(T,V,D):
    radius = (kB * T) / (6 * np.pi * V * diffusioncoef(gamma,q))
    return radius

def PDI(d,gamma):
    return d**2/gamma**2

print(diffusioncoef(gamma,q))
print(PDI(d, gamma))
print(radius(T,V,diffusioncoef(gamma,q)))

# z=[(1,2,3,4,5,6,7,8,9),(1,2,3,4,5,6,7,8,9)]

    # fit_results = np.array([('A', 'b', 'gamma', 'd'), (A, b, gamma, d)])
    # calc_results = np.array([('Diffusion Coeff.','P.D.I','Radius'), (diffusioncoef(gamma,q),PDI(d, gamma),radius(T,V,diffusioncoef(gamma,q)))])
    # system_parameters = np.array([('Angle (deg)','Theta (rad)','Temp. (K)','Viscosity (cp)','n','Wavelength (nm)'), (angle, theta, T, V, n, wavelength)])

#Export data to same destination as rawdata but appended with _data.
#Uncomment the below chunk if horizonal data is preferable.

    #horiz_data = np.array([('A', 'b', 'gamma', 'd','Diffusion Coeff.','P.D.I','Radius','Angle (deg)','Theta (rad)','Temp. (K)','Viscosity (cp)','n','Wavelength (nm)'), (A, b, gamma, d, diffusioncoef(gamma,q),PDI(d, gamma),radius(T,V,diffusioncoef(gamma,q)),angle, theta, T, V, n, wavelength)])
    #np.savetxt((result_path + name + '_horiz_data'+'.txt'), horiz_data, delimiter=",", fmt="%s")
    #print(type(result_path))

headers = np.array(('Fit parameters:', 'A', 'b', 'gamma', 'd',' ','Caculated Results:','Diffusion Coeff.','P.D.I','Radius','    ','System Parameters:' , 'Angle (deg)','Theta (rad)','Temp. (K)','Viscosity (cp)','Refractive Index','Wavelength (nm)'))
datas = np.array((' ', A, b, gamma, d,' ',' ',diffusioncoef(gamma,q),PDI(d, gamma),radius(T,V,diffusioncoef(gamma,q)),' ',' ', angle, theta, T, V, n, wavelength))
vertical_data = np.column_stack((headers,datas))

result_path = '/Users/alicesmith/Documents/Data/DLS PTCL/' + date + '/TPE Sols/' + name + '/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

np.savetxt((result_path + name + '_data'+'.txt'), vertical_data, delimiter="    ", fmt="%s")

print(name)
# np.savez((result_path + name + '_' + str(angle)), fit_results, calc_results, system_parameters)

# x = np.array(lag)
# y = np.array(ave)
# numlag = x.astype(np.float)
# numave = y.astype(np.float)
