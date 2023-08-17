# friction
A library containing routines for calculating the frictional response of contacting bodies.

## Work-In-Progress
This is a work in progress.  This library will eventually be a compilation of the friction modeling efforts I have undertaken as part of my everyday work.

## Available Models
- Coulomb Model
```math
F = \sign{ \left( v \right)} \mu_{c} N
```
- Lu-Gre Model
```math
F = \sigma_{0} z + \sigma_{1} e^{-\left( \frac{v}{v_s} \right)^{2} } \frac{dz}{dt} + \sigma_{2} v
```
```math
\frac{dz}{dt} = v - \frac{\left| v \right| z}{g(v)}
```
```math
g(v) = a_{1} + \frac{a_2}{1 + s^{\alpha}}
```
```math
a_{1} = \frac{\mu_c N}{\sigma_{0}}, a_{2} = \frac{\mu_s N + \mu_c N}{\sigma_{0}}, s = \frac{\left| v \right|}{v_s}
```

## References:
1. Al-Bender, Farid & Lampaert, Vincent & Swevers, Jan. (2004). Modeling of dry sliding friction dynamics: From heuristic models to physically motivated models and back. Chaos (Woodbury, N.Y.). 14. 446-60. 10.1063/1.1741752. 
2. Rizos, Demosthenis & Fassois, Spilios. (2009). Friction Identification Based Upon the LuGre and Maxwell Slip Models. Control Systems Technology, IEEE Transactions on. 17. 153 - 160. 10.1109/TCST.2008.921809. 