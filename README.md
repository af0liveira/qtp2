# Quantum Transport Properties 2 (QTP2)

This program calculates transport properties for a particle moving through a 1D potential energy barrier. 

The implmentation is based on the equations published by Poltavsky et al. (2018). 

The transmission coefficients for particles with energy $E$ moving through a potential barrier $U(z)$ are calculated as

$$
T(E) = \exp\left(-\frac{2}{\hbar} \sqrt{2m} \int_{z_{1}(E)}^{z_{2}(E)} \sqrt{U(z) - E} \ \mathrm{d}z \right)
$$

where $z_i(E)$ are the coordinates where $U(z_i) = E$. 
Not that $T(E) = 1$ if $E > U_\mathrm{max}$. 

The particle flux is calculated as

$$
u_\mathrm{total}(\beta) = u_\mathrm{class}(\beta) + u_\mathrm{tunnel}(\beta) 
$$

where $\beta$ is the reverse temperature defined as $\beta = 1/(k_\mathrm{B}T)$
and the classical and tunneling contributions are given as

$$
\begin{gather}
u_\mathrm{class}(\beta) = \frac{e^{-\beta U_\mathrm{max}}}{\sqrt{2\pi m \beta}} \\
u_\mathrm{tunnel}(\beta) = 
    \sqrt{\frac{\beta}{2\pi m}} \int_0^\infty T(E) e^{\beta U(z)} \mathrm{d}E
\end{gather}
$$

The average transmission probabilities $k$ are then computed as

$$
k_\alpha =  u_\alpha \sqrt{2\pi m \beta}
$$

where $\alpha$ indicates _classical_, _tunneling_, or _total_ probabilities. 

Finally, the Arrhenius equation 

$$
\ln k = \ln A - \beta E_\mathrm{a}
$$

is used to calculate the activation energies $E_\mathrm{a}$ and pre-exponential factors $A$. 
The activation energies are obtained from the first derivative of the Arrhenius equation above; once $E_\mathrm{a}$ is known, $A$ can be obtained. 

## Requirements

This project was developed on Python 3.10.6 using the following libraries

- numpy 1.23.3
- scipy 1.9.1
- pandas 1.4.2

## Setup

As long as the Python environment with the packages listed under _Requirements_, not special setup should be required. 
The program should work by simply excuting the `qtp2` command from the `bin` folder in the program's repository. 

A list of available options can be obtained by running
```sh
${REPO}/bin/qtp2 -h
``` 
where `${REPO}` is the folder containing the repository. 

## Running the program

The program runs from the command line by simply excuting the `qtp2` command within the `bin` folder in the repository. 
A full list of available options can be obtained by running
```sh
bin/qtp2 -h
```

### Input data

The program requires the following input information:

- mass of the particle traveling throught the potential energy barrier 
- the temperature (or range of temperatures) for the calculations 
- the potential energy barrier $U(z)$ as a XY file with $z$ values as X and $U(z)$ values as Y 
  
### Output data

The output contains three blocks of information

1. transmission coefficients associated to the potential energy barrier heights
2. the particle fluxes and average transmission probabilities (rate constants) for each temperature value
3. the Arrhenius activation energies and exponential prefactors -- note that at least two temperature values are required for this, though

* * *

[Poltavsky et al.(2018)] 
I Poltavsky et al. 
"Quantum tunneling of thermal protons through pristine graphene." 
_Journal of Chemical Physics_ 148(20), 204707 
(doi: 10.1063/1.5024317)

