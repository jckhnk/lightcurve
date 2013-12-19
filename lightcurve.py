import sys
import numpy as np
from numpy import pi, sin, cos
from scipy.special import jn


# functions
# ========================

def lommel(n, mu, nu, kmax=100):
	"""
	Approximates the $n^{th}$ order Lommel function $U_n(\mu, \nu)$ by 
	summing the infinite series from k to kmax.
	"""
	sum = 0
	for k in range(kmax+1):
		# val = (-1)**k * (mu/nu)**(n+2*k) * besselj(n+2*k, pi * mu * nu)
		val = (-1)**k * (mu/nu)**(n+2*k) * jn(n+2*k, pi * mu * nu)
		if abs(val) < 1e-15:
			sum += val
			break
		else:
			sum += val
	return sum

def intensity(rho, eta, kmax=100):
	"""
	The intensity of the diffraction pattern caused by the occultation
	of a distant star by a spherical foreground object, $I_{\rho}(\eta)$.
	The parameters $\rho$ and $\eta$ are the radius and distance from the
	line of sight in units of the Fresnel scale $F = (\lambda a / 2)^{1/2}$,
	respectively, i.e. $\rho = r/F$ and $\eta = x/F$.
	"""
	if eta <= rho:
		u0 = lommel(0, eta, rho, kmax=kmax)
		u1 = lommel(1, eta, rho, kmax=kmax)
		# return float(u0**2 + u1**2)
		return u0**2 + u1**2
	else:
		u1 = lommel(1, rho, eta, kmax=kmax)
		u2 = lommel(2, rho, eta, kmax=kmax)
		x = 0.5 * pi * (rho**2 + eta**2)
		# return float(1 + u1**2 + u2**2 - 2*u1*sin(x) + 2*u2*cos(x))
		return 1 + u1**2 + u2**2 - 2*u1*sin(x) + 2*u2*cos(x)

# unit conversion factors
# ========================
km_per_AU = 1.5e8
km_per_nm = 1.0e-12
radians_per_arcsec = 1.0/206265.

# default values
# ========================
NW = 5
# +/- tmax is lightcurve duration (sec)
tmax = 5
# sampling frequency (Hz)
sp = 40.0
# angular size of star (mas)
angsize = 0.01
# distance to oculting object in AU
d = 40
# impact parameter in fraction of stellar radius
b = 0.0
# radius of occulting object (km)
robj = 1
# transverse velocity (km/sec)
vt = 25.0
# number of event widths W to compute
m = 2.0
# number of outputs per Fresnel scale
n_per_F = 100
# name of output file containing lightcurve
outfile = 'lightcurve.txt'

# get input from command line
# ========================
for k in range(1,len(sys.argv)):
	if sys.argv[k] == '-robj':
		robj = sys.argv[k+1]
	if sys.argv[k] == '-d':
		d = sys.argv[k+1]
	if sys.argv[k] == '-star':
		angsize = sys.argv[k+1]
	if sys.argv[k] == '-outfile':
		outfile = sys.argv[k+1]

print "using: -robj",robj,"-d",d,"-star",angsize,"-outfile",outfile

# setup
# ========================
# r filter wavelength limits
w1 = 552
w2 = 689
# mean wavelength
mw = 0.5 * (w1 + w2)
# Fresnel scale (km)
F = np.sqrt(0.5 * mw * km_per_nm * d * km_per_AU)
# projected stellar radius
rstar = 0.5 * angsize / 1000 * radians_per_arcsec * d * km_per_AU
# projected stellar radius in Fresnel units
rho_star = rstar/F
# object radius in Fresnel units
rho = robj/F
# projected stellar radius
rstar = 0.5 * angsize /1000 * radians_per_arcsec * d * km_per_AU
# omega is the characteristic width of the occultation event in Fresnel units,
# Nihei et al. 2007
omega = 2. * (np.sqrt(3.0) ** 1.5 + rho ** 1.5) ** (2/3.) + 2. * rho_star
# characteristic width of occultation event, physical units
W = omega * F
# number of points in half-event
n = m * omega * n_per_F
# time span of half-event
tmax = m * W / vt
# distance moved in time span of half-event
xmax = vt * tmax
# maximum separation between center of occulter and a point on the face of
# the star.
sep_max = np.sqrt(xmax ** 2 + (b * W) ** 2) + rstar
dsep = sep_max / (n-1)
dt = dsep / vt

# allocate arrays
# ========================
g_sum = np.zeros(2 * n + 1)
I_array = np.zeros(n)
lightcurve = np.zeros(2 * n + 1)

# projected area of star (km^2)
# ========================
if 2*rstar < dsep:
	area = 1.
else:
	area = np.pi * rstar ** 2

# Calculate lightcurve values for an array of uniformly spaced values.
# These are averaged over a number of wavelength bins, assuming a neutral
# color.  The values get used in the final lightcurve calculation below.
# ========================
for i in range(int(n)):
	sep = i*dsep
	for j in range(NW):
		w = w1 + j*( (w2- w1) / NW)
		F = np.sqrt(0.5 * w * km_per_nm * d * km_per_AU)
		I_array[i] += intensity(robj/F, sep/F)
	I_array[i] /= NW

for i in range(int(n+1)):
	x = -xmax + xmax * i/n
	y = b * W
	eta = np.sqrt(x ** 2 + y ** 2)
	I = 0.0

	# star is unresolved, linearly interpolate to a single best value
	if 2 * rstar < dsep:
		idx = int(eta/dsep)
		eta_0 = dsep * idx
		I_0 = I_array[idx]
		eta_1 = eta_0 + dsep
		I_1 = I_array[idx+1]
		I_eta = I_0 + ( (I_1 - I_0) / dsep) * (eta - eta_0)
		g_sum[i] += 1.0
		I += I_eta

	else:	# star is resolved, integrate over face of the star
		if eta > rstar:	# center of object is outside of limb of star
			idx_min = int((eta - rstar)/dsep) + 1
			idx_max = int((eta  + rstar)/dsep)
		else:	# center of object is inside of limb of star
			idx_min = 1
			idx_max = int((eta  + rstar)/dsep)
		for j in range(idx_min,idx_max):
			eta_p = dsep * j
			if eta_p < (rstar - eta):
				g = 2 * np.pi
			else:
				g = 2*np.arccos((eta**2 + eta_p**2 - rstar**2)/(2.0*eta*eta_p))
			g_sum[i] += g * eta_p * dsep
			I += g * eta_p * I_array[j] * dsep

	g_sum[2*n-i] = g_sum[i]
	lightcurve[i] = I/g_sum[i]
	lightcurve[2*n-i] = I/g_sum[2*n-i]

# write results to disk
# ========================
x_a = np.array([-xmax + xmax * i/n for i in range(int(2*n+1))])
data_out = np.c_[x_a/vt, g_sum/area, lightcurve, x_a]
header = ' '.join([str(i) for i in [1.0/dt, W, F, rstar]])
np.savetxt(outfile,data_out,header=header)
