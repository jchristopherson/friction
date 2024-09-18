# friction
A library containing routines for calculating the frictional response of contacting bodies.

## Status
[![CMake](https://github.com/jchristopherson/friction/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/friction/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/friction/workflows/fpm/badge.svg)](https://github.com/jchristopherson/friction/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/friction/).

## Available Models
- Coulomb Model
```math
F = \text{sgn} \left( v \right) \mu_{c} N
```
- Lu-Gre Model
```math
F = \sigma_{0} z + \sigma_{1} \frac{dz}{dt} + \sigma_{2} v
```
```math
\frac{dz}{dt} = v - \frac{\left| v \right| z}{g(v)}
```
```math
g(v) = a_{1} + \frac{a_2}{1 + s^{\alpha}}
```
```math
a_{1} = \frac{\mu_c N}{\sigma_{0}}, a_{2} = \frac{\mu_s N - \mu_c N}{\sigma_{0}}, s = \frac{\left| v \right|}{v_s}
```
- Maxwell Model
```math
F = k \delta
```
```math
\delta_{i+1} = \text{sgn} \left( x_{i+1} - x_{i} + \delta_{i} \right) \min \left( \left| x_{i+1} - x_{i} + \delta_{i} \right|, \Delta \right)
```
```math
\Delta = \frac{N \mu_c}{k}
```
- Generalized Maxwell Slip Model
```math
F = \sum_{i=1}^{n} \left( k_i z_i + b_i \frac{dz_i}{dt} \right) + b_v v
```
```math
\begin{equation}
\frac{dz_i}{dt} = 
\begin{cases}
v & \text{if $|z_i| \le g(v)$} \\
\text{sgn} \left( v \right) \nu_i C \left( 1 - \frac{z_i}{\nu_i g(v)} \right) & \text{otherwise}
\end{cases}
\end{equation}
```
```math
g(v) = a_{1} + \frac{a_2}{1 + s^{\alpha}}
```
```math
a_{1} = \frac{\mu_c N}{\sigma_{0}}, a_{2} = \frac{\mu_s N - \mu_c N}{\sigma_{0}}, s = \frac{\left| v \right|}{v_s}
```
```math
\sum_{i=1}^{n} \nu_i = 1
```
- Stribeck Model
```math
F = \text{sgn} \left( v \right) \left( \mu_c N + N \left( mu_s - mu_c \right) \exp(-|\frac{v}{v_s}|^2) \right) + b_v v
```
- Modified Stribeck Model
```math
F = k \delta + b_v v
```
```math
\delta_{i+1} = \text{sgn} \left( x_{i+1} - x_{i} + \delta_{i} \right) \min \left( \left| x_{i+1} - x_{i} + \delta_{i} \right|, g(v) \right)
```
```math
g(v) = a_{1} + a_2 \exp(-|\frac{v}{v_s}|^2)
```
```math
a_{1} = \frac{\mu_c N}{k}, a_{2} = \frac{\mu_s N - \mu_c N}{k}
```

## References:
1. Al-Bender, Farid & Lampaert, Vincent & Swevers, Jan. (2004). Modeling of dry sliding friction dynamics: From heuristic models to physically motivated models and back. Chaos (Woodbury, N.Y.). 14. 446-60. 10.1063/1.1741752. 
2. Rizos, Demosthenis & Fassois, Spilios. (2009). Friction Identification Based Upon the LuGre and Maxwell Slip Models. Control Systems Technology, IEEE Transactions on. 17. 153 - 160. 10.1109/TCST.2008.921809. 
3. Al-Bender, Farid & Lampaert, Vincent & Swevers, Jan. (2005). The generalized Maxwell-Slip model: A novel model for friction simulation and compensation. Automatic Control, IEEE Transactions on. 50. 1883 - 1887. 10.1109/TAC.2005.858676. 
4. Tjahjowidodo, Tegoeh & Al-Bender, Farid & Brussel, H.. (2005). Friction identification and compensation in a DC motor. IFAC Proceedings Volumes. 32. 10.3182/20050703-6-CZ-1902.00093. 
5. Lampaert, Vincent & Al-Bender, Farid & Swevers, Jan. (2003). A generalized Maxwell-slip friction model appropriate for control purposes. 4. 1170- 1177 vol.4. 10.1109/PHYCON.2003.1237071. 
6. Al-Bender, Farid & Swevers, Jan. (2009). Characterization of friction force dynamics. Control Systems, IEEE. 28. 64 - 81. 10.1109/MCS.2008.929279. 
7. Al-Bender, Farid. (2010). Fundamentals of friction modeling. Proceedings - ASPE Spring Topical Meeting on Control of Precision Systems, ASPE 2010. 48. 
8. Al-Bender, Farid & Moerlooze, K.. (2011). Characterization and modeling of friction and wear: an overview. International Journal Sustainable Construction & Design. 2. 19-28. 10.21825/scad.v2i1.20431. 