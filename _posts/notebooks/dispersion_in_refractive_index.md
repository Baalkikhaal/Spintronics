---
layout: post
title: 
# Dispersion in refractive index

In this notebook, we visualize how the refractive index of Platinum varies with the wavelength of photon. This variation is called dispersion. Characteristic features include

- the general decreasing trend of $\delta$ and $\beta$ with the energy, that respectively suggest that medium behaves nearly as a vacuum for higher energy photons. This is termed **normal dispersion**.
- the real part of the refractive index is slightly less than 1 for photons more energetic than visible photons.
- Spikes in the imaginary part of the refractive index at characteristic energies that correspond to the lines of absorption for the material. This is termed **anomalous dispersion**.
- As we work with metallic ultra-thin films, we shall study the dispersion of refractive in Pt as a representative of metals. Data is sourced from the [atomic scattering files](https://henke.lbl.gov/optical_constants/asf.html) of the Centre for X-ray Optics from Berkeley Lab.
- To explain these features of dispersion in refractive index, we model the material as a collection of classical Lorentz oscillators with a dispersion in the density of these oscillators. The dispersion behaviour is also understood from a quantum mechanical counterpart of this theory.  


```python
import numpy as np
from scipy.constants import pi, Avogadro, physical_constants
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('../../../myMatplotlibStylesheet.mplstyle')
```


```python
from scipy.constants import h, c, e
```


```python
print(h*c/e)
```

    1.2398419843320026e-06



```python
#%% Constants
energy_CuKa = 8.051e3    # energy of photon corresponding to 1.54 Å.

z_pt, m_pt = 78, 195.078 # Atomic number and atomic weight of Pt.
rho_pt = 2.15e+4  # density in kg per m^3

r_e = physical_constants['classical electron radius'][0]          # classical radius of electron (in m)

```


```python
lambda_xrr = h*c/e/ energy_CuKa   # wavelength of photon used for XRR measurement (in m).
```


```python
lambda_xrr
```




    1.539985075558319e-10



## My TeX \newcommands

$$
\newcommand{\myScaSub}[2]{{#1}_\mathrm{#2}}
$$

$$
\newcommand{\myVec}[1]{\mathbf{ { \mathrm{#1} } } }
$$


## Classical theory of refractive index

We consider the interaction of light with matter in a simple model of an electron tied to the nucleus using a spring-like force with a spring constant $k$ and an associated damping characterized by the damping factor $\eta$. The force due to the electric field of light acts only upon the electron. The mechanics of the nucleus can be neglected owing to its larger mass.
The governing equation of motion is

$$
F_{drive} - kx - m \eta \dot{x} = m \ddot{x}
$$

Rewriting the spring constant in terms of the natural frequency of the oscillator $\omega_0$ given by

$$
k = m \omega_0^2
$$

we have

$$
F_{drive} - m \omega_0^2 x - m \eta \dot{x} = m \ddot{x}
$$

Consider a photon with frequency $\omega$ incident on this simple atom. The associated harmonic driving force is given by

$$
F_{drive} = - e E(t) = - e E_{0} e^{i \omega t}
$$

This force drives the electron, whose relative displacement $\delta x$ from the mean position can also be assumed to be of frequency within the linear response approximation so that

$$
\delta x = x_{0} e^{i \omega t}
$$

Substituting the harmonic expressions for the relative displacement and the force gives


$$
x_{0}(\omega) = - \frac{e E_{0}(\omega)}{m} \frac{1}{ \omega_0^2 - \omega^2 + \iota \eta \omega}
$$

This displacement gives rise to a dipole moment with the same frequency given by

$$
\myScaSub{\mu}{e} = -e x_{0}(\omega),
$$

so that if the number density of electrons is $N$, then the polarization, the total dipole moment per unit volume is given by

$$
P(\omega) = N \myScaSub{\mu}{e} = - N ex_{0}(\omega)
$$

This polarization can be expressed as linear response to the applied electric field with constant of proportionality, the susceptibility $\chi$ given by

$$
P(\omega) = \myScaSub{\epsilon}{0} \chi (\omega) E_{0}(\omega),
$$

so that

$$
\chi (\omega) = \frac{N e^2 }{\myScaSub{m \epsilon}{0}} \frac{1}{\omega_0^2 - \omega^2 + \iota \eta \omega}
$$

### Relative permitivity

Now this polarization leads to difference in the electric field $\myVec{E}$ and the displacement $\myVec{D}$ given by the Maxwell's equation

$$
\myVec{D} = \myScaSub{\epsilon}{0}\myVec{E} + \myVec{P}
$$

As a result, the response of a dielectric medium is polarization. This effect is experimentally measured in terms of relative permitivity $\myScaSub{\epsilon}{r}$ given by

$$
\myVec{D} = \myScaSub{\epsilon}{r}\myScaSub{\epsilon}{0}\myVec{E},
$$

so that

$$
\myScaSub{\epsilon}{r}(\omega) = 1 + \chi (\omega) 
$$

If we define the plasma frequency $\myScaSub{\omega}{p}$ as

$$
\myScaSub{\omega}{p}^2 = \frac{N e^2 }{\myScaSub{m \epsilon}{0}},
$$

then

$$
\myScaSub{\epsilon}{r}(\omega) = 1 + \frac{\myScaSub{\omega}{p}^2}{\omega^2 - \omega_0^2 + \iota \eta \omega}
$$

> This important equation says that relativity permitivity depends on the frequency of the field as well as the material's plasma frequency and damping factor.

### Estimate of relative permitivity

Using the typical values of [plasma frequency](../basics/optical_properties_of_metals.ipynb#Estimate-of-plasma-frequency), [damping factor](../basics/optical_properties_of_metals.ipynb#Estimate-of-damping-factor), and $\omega_0 = 0$ for metals, we have

$$
\myScaSub{\omega}{p} \sim 3 \,\mathrm{PHz},
$$

$$
\eta \sim 0.1 \myScaSub{\omega}{p},
$$

so that

$$
\myScaSub{\epsilon}{r} = 1 + \frac{\myScaSub{\omega}{p}^2}{- \omega^2  + \iota \eta \omega}
$$

If we resolve this complex quantity into real and imaginary components as

$$
\myScaSub{\epsilon}{r} = \myScaSub{\epsilon}{real}  + \iota \myScaSub{\epsilon}{imag},
$$

where

$$
\myScaSub{\epsilon}{real} = 1 - \frac{\myScaSub{\omega}{p}^2}{\omega^2  + \eta^2},
$$

and

$$
\myScaSub{\epsilon}{imag} = -  \frac{\myScaSub{\omega}{p}^2}{\omega}\frac{\eta}{\omega^2  + \eta^2}.
$$

In the special case of X-rays, for a wavelength of $1.54 \,\mathrm{A}^\circ$, the frequency is

nearly $2000 \, \mathrm{PHz} \sim 650 \myScaSub{\omega}{p}$. So unless there is anomalous dispersion, we have

$$
\myScaSub{\epsilon}{real} = 1 - \frac{1}{650^2  + 0.1^2} = 1- 2.4\cdot 10^{-6},
$$

and

$$
\myScaSub{\epsilon}{imag} = -  \frac{1}{650}\frac{0.1}{650^2  + 0.1^2} = -3.6\cdot 10^{-10}.
$$


```python
# Frequency of X-rays in PHz

f_xrays = c/lambda_xrr/1e+15
```


```python
f_xrays
```




    1946.7231388025675




```python
# Frequency of X-rays in units of w_p
w_p = 3
f_xrays/3
```




    648.9077129341891




```python
1/(650**2 + 0.1**2)
```




    2.3668638493049977e-06




```python
1/650*0.1/(650**2 + 0.1**2)
```




    3.6413289989307657e-10



### Refractive index



### Case: No damping

In the case of no damping, $\eta = 0$, then the amplitude of displacement is purely real. So that the the actual displacement of electron is oscillatory. For illustration let us consider practical values of $eE_0 = 1$~eV/nm 

## Quantum theory of refractive index

Reference: Paratt Phys Rev [1954]

The quantum theory of radiation accounts for the dispersion of the atomic scattering factors, by assuming a distribution of "dispersion" oscillators. The number of oscillators $\mathrm{d}g$ within a window of $\omega$ and $\omega + \mathrm{d}\omega$ is 

$$
\mathrm{d}g = \Gamma(\omega) \mathrm{d}\omega.
$$

The response of the medium is then governed by the effective number $g$ of oscillators in the vicinity of the radiation, also called the **oscillator strength**. If we know the oscillator density, we can have a fair idea about the optical response of the medium. Now this density is related to the photo electric absorption coefficient $\mu $ by

$$
\Gamma(\omega) = \frac{m_e c}{2 \pi^2 e^2} \mu (\omega),
$$

so that if we know the dispersion of the photo-electric coefficient, we know the optical response of the medium.


Here we make an **assumption** that the dispersion of the photo-electric coefficient follows a power law. Further, we invoke the **quantum nature** of the interaction that forbids absorption of radiation with frequencies smaller than the characteristic frequency $\omega_q$ of the electrons. This is due to the quantization of the electronic motion around the nucleus in the form of shells. If we denote the characteristic frequency of the $q^\mathrm{th}$ shell of atom by $\omega_q$, then

$$
\mu (\omega) = \left(\frac{\omega_q}{\omega}\right)^{p_q} \mu_q \quad (\omega \geq \omega_q)
$$

where $p_q$ is the exponent of the power law distribution of the oscillators and is typically taken to be 3 for calculation purposes. However, it has different values for different shells as


\begin{array}{ll}
K & \frac{11}{4} \\
L & \frac{7}{3}  \\
M & \frac{5}{2}
\end{array}


```python
#%% Extract the data
raw_data_pt = 'pt.nff'     # atomic scattering factors
data = np.loadtxt(raw_data_pt, skiprows=1, unpack=True)
```


```python
energies = data[0]
f1s = data[1]
f2s = data[2]
```

## Anomalous dispersion

The normal dispersion monotonically increases with the energy of the photon. However, in the **left** neighbourhood vicinity of the characteristic frequencies associated with the absorption of photon by the inner electronic shells, there is a decrease of the dispersion with the energy, also called **anomalous** dispersion.

For Pt, the characteristic frequencies of photon absorption correspond to the ejection of electrons from the inner K, L, M shells and are given by

\begin{array}{ll}
\mathrm{line} & eV \\
K_1 & 78394.8 \\
L_1 & 13880.7 \\
L_2 & 13271.9 \\
L_3 & 11562.8 \\
M_1 & 3296.0  \\
M_2 & 3026.5 \\
M_3 & 2645.4 \\
M_4 & 2201.9 \\
M_5 & 2121.6 \\
N_1 & 725.4  \\
N_2 & 609.1  \\
N_3 & 519.4  \\
N_4 & 331.6  \\
N_5 & 314.6  \\
N_6 & 74.5   \\
N_7 & 71.2   \\
O_1 & 101.7  \\
O_2 & 65.3   \\
O_3 & 51.7   \\
\end{array}


```python
lines = [
            [r'$K_1$',
             r'$L_1$', r'$L_2$', r'$L_3$',
             r'$M_1$', r'$M_2$', r'$M_3$', r'$M_4$', r'$M_5$',
             r'$N_1$', r'$N_2$', r'$N_3$', r'$N_4$', r'$N_5$', r'$N_6$', r'$N_7$',
             r'$O_1$', r'$O_2$', r'$O_3$'],
            [78394.8,
             13880.7, 13271.9, 11562.8,
             3296, 3026.5, 2645.4, 2201.9, 2121.6,
             725.4, 609.1, 519.4, 331.6, 314.6, 74.5, 71.2,
             101.7, 65.3, 51.7
            ]
        ]
```


```python
f1s[np.argmin((energies - lines[1][1])**2)]
```




    68.8083




```python
#%% Plot the dispersion of scattering factor with energy.

fig, ax = plt.subplots()

x_min = 5e+1
ax.plot(energies[energies > x_min], f1s[energies > x_min], '-', label=r'$f_1^\mathrm{Pt}$')
ax.plot(energies[energies> x_min], f2s[energies> x_min], '-', label=r'$f_2^\mathrm{Pt}$')

ax.axvline(energy_CuKa, color='C2')
ax.text(energy_CuKa, 0.2, r'$E=8.05$ keV $\equiv 1.54\,\r{A}$ ', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.axhline(z_pt, color='C3', ls = '--')
ax.text(0.2, z_pt, r'$Z_\mathrm{Pt}$ ', transform=ax.get_yaxis_transform(),
        va='top', rotation='horizontal')

# Annotate the anomalous dispersion peaks.
#line_energy = lines[1][1]
#line_name = lines[0][1]
#ax.text(line_energy, f1s[np.argmin((energies - lines[1][1])**2)], line_name,
#        #transform=ax.get_xaxis_transform(),
#        va='top', rotation='horizontal')

for line_name, line_energy in zip (lines[0], lines[1]):
    ax.text(line_energy, f1s[np.argmin((energies - line_energy)**2)], line_name,
            #transform=ax.get_xaxis_transform(),
            va='center', rotation='horizontal')


ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

# Set the ticks in log scale.
major_ticks = mpl.ticker.LogLocator(base=10.0)
minor_ticks = mpl.ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=10)

#ax.xaxis.set_major_locator(major_ticks)
#ax.xaxis.set_minor_locator(minor_ticks)
#ax.yaxis.set_major_locator(major_ticks)
#ax.yaxis.set_minor_locator(minor_ticks)

ax.set_xlabel(r'$E$ (in eV)')
ax.set_ylabel(r'$f_1, f_2$')

```




    Text(0, 0.5, '$f_1, f_2$')




    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_32_1.png)
    



```python
fig.savefig('dispersion_scattering_factor_Pt.png')
fig.savefig('dispersion_scattering_factor_Pt.pdf')
```


```python
# Get the lambdas from energies.

lambdas = np.array([h*c/e/ energy for energy in energies])
```


```python
#%% Plot the dispersion of scattering factor with wavelength.

fig, ax = plt.subplots()

lambda_max = 1e-8
ax.plot(lambdas[lambdas < lambda_max]/1e-10, f1s[lambdas < lambda_max], '-', label=r'$f_1^\mathrm{Pt}$')
ax.plot(lambdas[lambdas < lambda_max]/1e-10, f2s[lambdas < lambda_max], '-', label=r'$f_2^\mathrm{Pt}$')

ax.axvline(lambda_xrr/1e-10, color='C2')
ax.text(lambda_xrr/1e-10, 0.2, r'$E=8.05$ keV $\equiv 1.54\,\r{A}$ ', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.axhline(z_pt, color='C3', ls = '--')
ax.text(0.2, z_pt, r'$Z_\mathrm{Pt}$ ', transform=ax.get_yaxis_transform(),
        va='top', rotation='horizontal')

ax.legend()
ax.set_xscale('log')
#ax.set_yscale('log')

# Set the ticks in log scale.
major_ticks = mpl.ticker.LogLocator(base=10.0)
minor_ticks = mpl.ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=10)

#ax.xaxis.set_major_locator(major_ticks)
#ax.xaxis.set_minor_locator(minor_ticks)
#ax.yaxis.set_major_locator(major_ticks)
#ax.yaxis.set_minor_locator(minor_ticks)

ax.set_xlabel(r'$\lambda$ (in $\r{A}$)')
ax.set_ylabel(r'$f_1, f_2$')

```




    Text(0, 0.5, '$f_1, f_2$')




    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_35_1.png)
    


## Convert atomic scattering factor into refractive index

Using the relation

$$
\delta + \iota \beta =  \frac{r_e \lambda^2 N_A}{2\pi M_a}  \rho  \left( f_1 + \iota f_2 \right)
$$

where $\myScaSub{M}{a}$ is the atomic mass, $\rho$ is the density of material in \si{\gram\per\centi\metre\cubed}, $\myScaSub{r}{e} = 2.818\, \mathrm{nm}$ is the classical radius of electron, $\lambda=$ 1.54 Å is the wavelength of photon. This gives a relation between the atomic scattering factors, the density and the refractive index as

$$
\delta + \iota \beta =  6.406 \cdot 10^{-6}  \rho \frac{  f_1 + \iota f_2 }{M_a}
$$



### Relation between energy and wavelength of photon

$$
E \, \mathrm{ (in eV)} = \frac{1.24}{\lambda \, \mathrm{ (in mum)}}
$$

## Relation between refractive index and the scatering factor

Reference: Paratt, Phys Rev[1954]

As per the quantum theory of radiation, the refractive index for X-rays is given as

$$
n = 1 - (\delta + \iota \beta)
$$

where $\delta$ for a medium containing different atomic species $s$ is

$$
\delta = \sum_{s} \delta_{s},
$$

where

$$
\delta_{s} = \sum_{s} A_{s} \lambda^2 f_1,
$$

where $f_1$ is the atomic scattering factor and

$$
A_{s} = \frac{N_A r_e}{2\pi} \frac{\rho_{s}}{M_{s}}\cdot 10^{3}, 
$$

is the prefactor corresponding to the material, where $N_A$ is the Avogadro number, $r_e = e^2/m_e c^2 = 2.8 \cdot 10^{-15}\, \mathrm{m}$ is the classical electron radius, $\rho_s, M_s$ are the density in $\mathrm{kg/m}^3$ and atomic mass in $\mathrm{Da}$ of the species $s$ 
for very small scattering angle (close to zero), where $f_1$ is the atomic scattering factor.


```python
# Estimate the refractive index from atomic scattering factors.

a_s = Avogadro * r_e /(2*pi) * rho_pt / m_pt * 1e+3
```


```python
a_s/1e+13
```




    2.9766857718838953



The prefactor $A_{s} \sim 3\cdot10^{13}\, \mathrm{m}^{-2}$ so that for X-ray from $K_\alpha^\mathrm{Cu}=1.54 A^{\!0}$, the refractive index in terms of atomic scattering factor $f_1= 73.66 + \iota 7.2356$ as 


```python
f1_pt, f2_pt = 73.66, 7.232
delta_pt, beta_pt = [a_s  * lambda_xrr**2 * f_pt for f_pt in [f1_pt, f2_pt]]
```


```python
delta_pt, beta_pt
```




    (5.1999327869637284e-05, 5.105337213592409e-06)



## $\delta$'s and $\beta$'s from scattering factors

As seen above there is a $\lambda^2$ prefactor in the dependence of $\delta$ on the scattering factor $f_1$. 


```python
deltas = np.array([a_s * lamb**2 * f1 for lamb, f1  in zip(lambdas, f1s)])
betas = np.array([a_s * lamb**2 * f2 for lamb, f2  in zip(lambdas, f2s)])
```


```python
deltas
```




    array([-4.57532804e+03, -4.43087500e+03, -4.29091177e+03, -4.15541702e+03,
           -4.02420143e+03, -3.89708929e+03, -3.77398584e+03, -3.65485996e+03,
           -3.53941808e+03, -3.42763919e+03, -3.31937600e+03, -3.21454499e+03,
           -3.11306148e+03, -3.01474216e+03, -2.91951500e+03, -2.82735103e+03,
           -2.73804246e+03, -2.65156678e+03, -2.56785395e+03, -2.48675997e+03,
           -2.40822419e+03, -2.33215164e+03, -2.25851633e+03, -2.18719472e+03,
           -2.11810321e+03, -2.05121783e+03, -1.98645616e+03, -1.92371548e+03,
           -1.86297336e+03, -1.80413408e+03, -1.74715327e+03, -1.69198589e+03,
           -1.63854714e+03, -1.58679670e+03, -1.53669358e+03, -1.48816239e+03,
           -1.44116623e+03, -1.39565220e+03, -1.35156992e+03, -1.30889936e+03,
           -1.26756400e+03, -1.22753288e+03, -1.18876242e+03, -1.15122247e+03,
           -1.11486050e+03, -1.07965932e+03, -1.04555938e+03, -1.01254470e+03,
           -9.80569870e+02, -9.49601336e+02, -9.19615113e+02, -8.90571179e+02,
           -8.62447237e+02, -8.35206399e+02, -8.08834652e+02, -7.83290471e+02,
           -7.58554401e+02, -7.34600806e+02, -7.11399350e+02, -6.88932301e+02,
           -6.67176576e+02, -6.46105014e+02, -6.25701368e+02, -6.05944758e+02,
           -5.86806578e+02, -5.68276478e+02, -5.50327843e+02,  2.53186471e-01,
            2.53913651e-01,  2.56124847e-01,  2.59834740e-01,  2.65899565e-01,
            2.77267332e-01,  2.86076052e-01,  2.90492585e-01,  2.94612212e-01,
            2.92930782e-01,  2.87932093e-01,  2.80286964e-01,  2.70712747e-01,
            2.61781852e-01,  2.52714354e-01,  2.43170533e-01,  2.33049834e-01,
            2.22537941e-01,  2.11340616e-01,  1.98440910e-01,  1.83528317e-01,
            1.72015775e-01,  1.61997127e-01,  1.52515912e-01,  1.43107332e-01,
            1.32905147e-01,  1.24016206e-01,  1.15976598e-01,  1.08209864e-01,
            9.97256613e-02,  8.97635271e-02,  8.30301706e-02,  7.82142430e-02,
            7.47229113e-02,  7.26208586e-02,  7.12585183e-02,  7.07864149e-02,
            7.15771774e-02,  7.44207504e-02,  7.77157137e-02,  8.06904184e-02,
            8.46683767e-02,  9.03966082e-02,  9.40889980e-02,  9.73891893e-02,
            1.00207622e-01,  1.00218839e-01,  9.78121730e-02,  9.60047592e-02,
            9.38045171e-02,  8.99406449e-02,  8.61781999e-02,  8.43583440e-02,
            8.34921876e-02,  8.57064403e-02,  9.18507021e-02,  9.94173509e-02,
            1.06327860e-01,  1.15146788e-01,  1.21365031e-01,  1.24835047e-01,
            1.26582349e-01,  1.27246983e-01,  1.27458207e-01,  1.26841568e-01,
            1.25285962e-01,  1.22905863e-01,  1.20605297e-01,  1.18723417e-01,
            1.17178567e-01,  1.14265276e-01,  1.10552489e-01,  1.05937947e-01,
            1.01851419e-01,  9.80545873e-02,  9.43148767e-02,  9.05285107e-02,
            8.68562101e-02,  8.34239515e-02,  8.01314247e-02,  7.69427510e-02,
            7.38233657e-02,  7.08155682e-02,  6.78948916e-02,  6.50389667e-02,
            6.22413910e-02,  5.95149833e-02,  5.68251655e-02,  5.41877760e-02,
            5.15808835e-02,  4.88374406e-02,  4.62679411e-02,  4.39725709e-02,
            4.18147170e-02,  3.96960427e-02,  3.77285144e-02,  3.58914855e-02,
            3.41341371e-02,  3.24130708e-02,  3.06905640e-02,  2.91235108e-02,
            2.76571714e-02,  2.62382985e-02,  2.48691895e-02,  2.36161257e-02,
            2.24465060e-02,  2.13359072e-02,  2.02420970e-02,  1.92543249e-02,
            1.83583306e-02,  1.75411383e-02,  1.68100183e-02,  1.61300363e-02,
            1.54937050e-02,  1.49065273e-02,  1.43770643e-02,  1.39412397e-02,
            1.35498923e-02,  1.31358922e-02,  1.27249577e-02,  1.23243712e-02,
            1.19290090e-02,  1.15417083e-02,  1.12047984e-02,  1.09144622e-02,
            1.06666391e-02,  1.04935158e-02,  1.03228565e-02,  1.01599535e-02,
            1.00258965e-02,  9.93792679e-03,  9.94853367e-03,  9.89786943e-03,
            9.81455077e-03,  9.72356963e-03,  9.64246088e-03,  9.57704773e-03,
            9.48131246e-03,  9.36241425e-03,  9.23411968e-03,  9.09815927e-03,
            8.95223139e-03,  8.81117407e-03,  8.67381384e-03,  8.53880258e-03,
            8.40584937e-03,  8.27358395e-03,  8.14611890e-03,  8.02347666e-03,
            7.90474986e-03,  7.78956300e-03,  7.67831832e-03,  7.57165560e-03,
            7.46943707e-03,  7.37178394e-03,  7.27934774e-03,  7.19238735e-03,
            7.11157881e-03,  7.03786613e-03,  6.97306657e-03,  6.92112859e-03,
            6.89155525e-03,  6.88326317e-03,  6.87456924e-03,  6.82913051e-03,
            6.76381872e-03,  6.68727830e-03,  6.60322040e-03,  6.51406059e-03,
            6.42146924e-03,  6.33010735e-03,  6.23514399e-03,  6.13585053e-03,
            6.03084658e-03,  5.92128334e-03,  5.80842041e-03,  5.69191001e-03,
            5.57129332e-03,  5.44422825e-03,  5.29955003e-03,  5.14802566e-03,
            5.00703404e-03,  4.90959646e-03,  4.85623056e-03,  4.83089008e-03,
            4.79829096e-03,  4.74483670e-03,  4.68201497e-03,  4.61185352e-03,
            4.53790094e-03,  4.46241797e-03,  4.38846625e-03,  4.31562029e-03,
            4.24116991e-03,  4.15913323e-03,  4.07129458e-03,  3.97923215e-03,
            3.88678771e-03,  3.79380296e-03,  3.70182341e-03,  3.61738359e-03,
            3.53865124e-03,  3.46468459e-03,  3.39638469e-03,  3.33566876e-03,
            3.27892060e-03,  3.21478992e-03,  3.14646896e-03,  3.07661847e-03,
            3.00630971e-03,  2.93606570e-03,  2.86613085e-03,  2.79641085e-03,
            2.72725194e-03,  2.65898413e-03,  2.59204906e-03,  2.52602608e-03,
            2.46076321e-03,  2.39621591e-03,  2.33267301e-03,  2.27036122e-03,
            2.20953577e-03,  2.14956016e-03,  2.09033284e-03,  2.03185339e-03,
            1.97434407e-03,  1.91789276e-03,  1.86255174e-03,  1.80832964e-03,
            1.75523710e-03,  1.70334089e-03,  1.65262399e-03,  1.60308266e-03,
            1.55475875e-03,  1.50762093e-03,  1.46160043e-03,  1.41653948e-03,
            1.37244276e-03,  1.32934192e-03,  1.28725543e-03,  1.24616696e-03,
            1.20605393e-03,  1.16690486e-03,  1.12871389e-03,  1.09145858e-03,
            1.05513872e-03,  1.01962293e-03,  9.84865698e-04,  9.50854355e-04,
            9.17583281e-04,  8.85019337e-04,  8.53129687e-04,  8.21857964e-04,
            7.91167560e-04,  7.61065867e-04,  7.31402335e-04,  7.02039497e-04,
            6.72874692e-04,  6.43794240e-04,  6.14645023e-04,  5.85212729e-04,
            5.55162277e-04,  5.23987008e-04,  4.90844992e-04,  4.54286818e-04,
            4.10748663e-04,  3.48685016e-04,  5.46168839e-06,  5.44576737e-06,
            2.06660574e-04,  3.15530484e-04,  2.80954899e-04,  1.23865307e-04,
            1.24038152e-04,  3.46089827e-04,  3.82046147e-04,  3.98494594e-04,
            4.06330329e-04,  4.09048462e-04,  4.08298969e-04,  4.04971998e-04,
            3.99558621e-04,  3.92207149e-04,  3.82555765e-04,  3.68018667e-04,
            3.10161039e-04,  3.10150757e-04,  3.56802108e-04,  3.63300614e-04,
            3.61157967e-04,  3.56485884e-04,  3.50450647e-04,  3.43476824e-04,
            3.35694182e-04,  3.26846139e-04,  3.10692616e-04,  3.02405819e-04,
            3.02379341e-04,  3.12217964e-04,  3.06509876e-04,  2.99730490e-04,
            2.92289638e-04,  2.83518959e-04,  2.73214620e-04,  2.73190730e-04,
            2.77685648e-04,  2.72222071e-04,  2.66031311e-04,  2.59620249e-04,
            2.53128837e-04,  2.46625703e-04,  2.40156132e-04,  2.33743651e-04,
            2.27402757e-04,  2.21149437e-04,  2.14996695e-04,  2.08951927e-04,
            2.03018512e-04,  1.97200735e-04,  1.91506383e-04,  1.85934383e-04,
            1.80486338e-04,  1.75164218e-04,  1.69968940e-04,  1.64899553e-04,
            1.59956069e-04,  1.55137304e-04,  1.50442455e-04,  1.45870716e-04,
            1.41420106e-04,  1.37089038e-04,  1.32875812e-04,  1.28778165e-04,
            1.24794317e-04,  1.20922598e-04,  1.17159881e-04,  1.13504503e-04,
            1.09953993e-04,  1.06506341e-04,  1.03158959e-04,  9.99094322e-05,
            9.67560776e-05,  9.36955645e-05,  9.07263785e-05,  8.78454932e-05,
            8.50512522e-05,  8.23412462e-05,  7.97127145e-05,  7.71642148e-05,
            7.46927165e-05,  7.22970230e-05,  6.99741642e-05,  6.77226072e-05,
            6.55399830e-05,  6.34244534e-05,  6.13738917e-05,  5.93864951e-05,
            5.74602773e-05,  5.55933300e-05,  5.37839507e-05,  5.20302181e-05,
            5.03304370e-05,  4.86829982e-05,  4.70858346e-05,  4.55375014e-05,
            4.40363453e-05,  4.25808329e-05,  4.11689088e-05,  3.97993660e-05,
            3.84702886e-05,  3.71800866e-05,  3.59270737e-05,  3.47093600e-05,
            3.35273419e-05,  3.23759995e-05,  3.12473541e-05,  3.01426190e-05,
            2.90567285e-05,  2.79833157e-05,  2.69109686e-05,  2.58183838e-05,
            2.46516046e-05,  2.31783879e-05,  1.78751062e-05,  1.78748300e-05,
            2.22561046e-05,  2.23825476e-05,  2.20068615e-05,  2.14958364e-05,
            2.09223092e-05,  2.03088409e-05,  1.96583388e-05,  1.89437096e-05,
            1.78944267e-05,  1.62327304e-05,  1.62323191e-05,  1.77712898e-05,
            1.74098620e-05,  1.63429323e-05,  1.63428414e-05,  1.67106844e-05,
            1.66021431e-05,  1.62354386e-05,  1.58375973e-05,  1.54287048e-05,
            1.50173616e-05,  1.46073644e-05,  1.42013795e-05,  1.38008890e-05,
            1.34069389e-05,  1.30202104e-05,  1.26414887e-05,  1.22708908e-05,
            1.19087805e-05,  1.15552196e-05,  1.12101367e-05,  1.08738897e-05,
            1.05461919e-05,  1.02270791e-05,  9.91653359e-06,  9.61433940e-06,
            9.32045022e-06,  9.03463016e-06,  8.75690946e-06,  8.48703598e-06,
            8.22490131e-06,  7.97026650e-06,  7.72304223e-06,  7.48297051e-06,
            7.25001710e-06,  7.02394747e-06,  6.80458374e-06,  6.59174388e-06,
            6.38527571e-06,  6.18506449e-06,  5.99088351e-06,  5.80257491e-06,
            5.61998428e-06,  5.44299611e-06,  5.27142338e-06,  5.10514714e-06,
            4.94394796e-06,  4.78772292e-06,  4.63633805e-06,  4.48965636e-06,
            4.34751382e-06,  4.20984593e-06,  4.07696722e-06,  3.94557283e-06])




```python
#%% Plot the data.

fig, ax = plt.subplots()

x_min = 5e+1
ax.plot(energies[energies > x_min], deltas[energies > x_min], '-', label=r'$\delta_\mathrm{Pt}$')
ax.plot(energies[energies> x_min], betas[energies> x_min], '-', label=r'$\beta_\mathrm{Pt}$')

ax.axvline(energy_CuKa, color='C2')
ax.text(energy_CuKa, 0.2, r'$E=8.05$ keV $\equiv 1.54\,\r{A}$ ', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

# Set the ticks in log scale.
major_ticks = mpl.ticker.LogLocator(base=10.0)
minor_ticks = mpl.ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=10)

#ax.xaxis.set_major_locator(major_ticks)
#ax.xaxis.set_minor_locator(minor_ticks)
ax.yaxis.set_major_locator(major_ticks)
ax.yaxis.set_minor_locator(minor_ticks)

ax.set_xlabel(r'$E$ (in eV)')
ax.set_ylabel(r'$n(E)= 1 - \delta(E) + i \beta(E)$')
```




    Text(0, 0.5, '$n(E)= 1 - \\delta(E) + i \\beta(E)$')




    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_47_1.png)
    



```python
fig.savefig('dispersion_refractive_index_Pt.png')
fig.savefig('dispersion_refractive_index_Pt.pdf')
```


```python
from scipy.constants import epsilon_0
```


```python
epsilon_0
```




    8.8541878128e-12



## Direct plotting of refractive indices 

CXRO also provides database of the refractive indices.

Let us plot them.


```python
raw_data_n = 'raw_data_refractive_index_Pt.txt'
```


```python
#%% Extract the data
data_n = np.loadtxt(raw_data_n, skiprows=2, unpack=True)
print(data.shape)
```

    (3, 516)



```python
energies = data_n[0]
deltas = data_n[1]
betas = data_n[2]
```


```python
#%% Plotc the data.

fig, ax = plt.subplots()

ax.plot(energies, deltas, '-', label=r'$\delta_\mathrm{Pt}$')
ax.plot(energies, betas, '-', label=r'$\beta_\mathrm{Pt}$')

ax.axvline(energy_CuKa, color='C2')
ax.text(energy_CuKa, 0.2, r'$E=8.05$ keV $\equiv 1.54\,\r{A}$ ', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

# Set the ticks in log scale.
major_ticks = mpl.ticker.LogLocator(base=10.0)
minor_ticks = mpl.ticker.LogLocator(base=10.0, subs=np.arange(1, 10), numticks=10)

#ax.xaxis.set_major_locator(major_ticks)
#ax.xaxis.set_minor_locator(minor_ticks)
ax.yaxis.set_major_locator(major_ticks)
ax.yaxis.set_minor_locator(minor_ticks)

ax.set_xlabel(r'$E$ (in eV)')
ax.set_ylabel(r'$n(E)= 1 - \delta(E) + i \beta(E)$')
fig.subplots_adjust(left=0.25, top=0.9)
#ax.set_title(r'Material(Pt)')
```


    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_55_0.png)
    



```python
fig.savefig('refractive_index_Pt.png')
fig.savefig('refractive_index_Pt.pdf')
```

## Comparison of refractive indices of Co and Pt

Let us compare the refractive indices of different metals.


```python
raw_datas = [
                'raw_data_refractive_index_Ta.txt',
                'raw_data_refractive_index_Pt.txt',
                'raw_data_refractive_index_Co.txt',
                #'raw_data_refractive_index_Cu.txt',
                'raw_data_refractive_index_Fe.txt',
                'raw_data_refractive_index_B.txt'
            ]
labels = [
            r'$\delta_\mathrm{Ta}$',
            r'$\delta_\mathrm{Pt}$',
            r'$\delta_\mathrm{Co}$',
            #r'$\delta_\mathrm{Cu}$',
            r'$\delta_\mathrm{Fe}$',
            r'$\delta_\mathrm{B}$'
         ]
```


```python
#%% Extract the data
data_n_s = [np.loadtxt(raw_data, skiprows=2, unpack=True) for raw_data in raw_datas]
```


```python
energies = [data_n[0] for data_n in data_n_s]
deltas = [data_n[1] for data_n in data_n_s]
```


```python
#%% Plot the deltas.

fig, ax = plt.subplots()

xmin = 1e+2

for energy, delta, label in zip(energies, deltas, labels):
    ax.plot(energy[energy > xmin], delta[energy > xmin], '-', label=label)

#ax.plot(energies[2], deltas[2], '-', label=labels[2])

ax.axvline(energy_CuKa, color='C2')
ax.text(energy_CuKa, 0.2, r'$E=8.05$ keV', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E$ (in eV)')
ax.set_ylabel(r'$\delta(E)$')
fig.subplots_adjust(left=0.25, top=0.9)
#ax.set_title(r'Material(Pt)')
```


    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_61_0.png)
    



```python
fig.savefig('deltas_Ta_Pt_Co_Fe_B.png')
fig.savefig('deltas_Ta_Pt_Co_Fe_B.pdf')
```

### Refractive index of synthetic $Co_3FeB$

In a multi atom material, the refractive index is the weighted contributions from each of the atomic species.

Since the chemical formula of CoFeB alloy is $Co_{60}Fe_{20}B_{20}$ (at%), we have

$$
f_\mathrm{CoFeB} =  \frac{60 f_\mathrm{Co} + 20 f_\mathrm{Fe} + 20 f_\mathrm{B}}{100}, 
$$

and

$$
M_\mathrm{CoFeB} = \frac{60 M_\mathrm{Co} + 20 M_\mathrm{Fe} + 20 M_\mathrm{B}}{100}.
$$


```python
raw_datas = [
                'raw_data_refractive_index_Cu.txt',
                'raw_data_refractive_index_Co.txt',
                'raw_data_refractive_index_CoFeB.txt'            ]
labels = [
            r'$\delta_\mathrm{Cu}$',
            r'$\delta_\mathrm{Co}$',
            r'$\delta_\mathrm{CoFeB}$'
         ]
```


```python
#%% Extract the data
data_n_s = [np.loadtxt(raw_data, skiprows=2, unpack=True) for raw_data in raw_datas]
```


```python
energies = [data_n[0] for data_n in data_n_s]
deltas = [data_n[1] for data_n in data_n_s]
```


```python
#%% Plot the deltas.

fig, ax = plt.subplots()

xmin = 5e+3

for energy, delta, label in zip(energies, deltas, labels):
    ax.plot(energy[energy > xmin], delta[energy > xmin], '-', label=label)

#ax.plot(energies[2], deltas[2], '-', label=labels[2])

ax.axvline(energy_CuKa, color='C2')
ax.text(energy_CuKa, 0.2, r'$E=8.05$ keV', transform=ax.get_xaxis_transform(),
        rotation='vertical')

ax.legend()
#ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(right=1e+4)
ax.set_xlabel(r'$E$ (in eV)')
ax.set_ylabel(r'$\delta(E)$')
fig.subplots_adjust(left=0.25, top=0.9)
#ax.set_title(r'Material(Pt)')
```


    
![Jupyter Notebook Plot](/assets/notebooks/dispersion_in_refractive_index_files/dispersion_in_refractive_index_67_0.png)
    



```python
(60*24.6213 + 20*24.8575 + 20*5.00961)/100
```




    20.746202000000004




```python
(60*3.56198 + 20*3.20844 + 20 * 0.0040748)/100
```




    2.7796909600000004




```python
(60* 59 + 20* 56 + 20 * 11)/100
```




    48.8



## Comparison of $\delta$ and $\rho$ of the material

Given the $\delta$ and $f_1$, let us estimate the density $\rho$ of the material.

Since

$$
\delta_\mathrm{Pt} =  A_\mathrm{Pt} \lambda^2 f_1^\mathrm{Pt},
$$

we have

$$
A_\mathrm{Pt} = \frac{\delta_\mathrm{Pt}}{\lambda^2 f_1^\mathrm{Pt}}.
$$


```python
f_pt = 73.6633
d_pt = 3.8e-5

a_pt = d_pt/ lambda_xrr**2 / f_pt
```


```python
a_pt
```




    21751588784184.492



Since 

$$
A_\mathrm{Pt}  = \frac{N_A r_e}{2\pi} \frac{\rho_\mathrm{Pt} }{M_\mathrm{Pt} }, 
$$

we have

$$
\rho_\mathrm{Pt} =  \frac{2\pi}{N_A r_e} A_\mathrm{Pt} M_\mathrm{Pt} 
$$


```python
rho_pt_estimate = 2 * pi / (Avogadro* r_e) * a_pt * m_pt
```


```python
rho_pt_estimate
```




    15704451.3870913



The dispersion of atomic scattering factor for atomic species $s$ extracted from the experimentally observed dispersion in refractive index is given by

$$
f_1 = \frac{\delta}{A_{s} \lambda^2}.
$$


```python

```
