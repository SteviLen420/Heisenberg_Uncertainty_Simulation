## Results and Discussion

This chapter presents the results of a numerical simulation of Gaussian wave packets, with a particular focus on the quantitative and qualitative verification of the Heisenberg uncertainty principle. The analysis covers the investigation of initial states, the behavior of the uncertainty product, and the time evolution of the system.

### Presentation of Results

In the simulation, we examined the behavior of quantum mechanical wave packets in one dimension. The initial states were described by Gaussian functions, the standard deviation ($σ_x$) of which was systematically varied.

**1. The Relationship Between Position and Momentum Uncertainty**

The `heisenberg_scan.csv` dataset and the `uncertainty_product.png` plot generated from it clearly demonstrate the inverse relationship between position uncertainty ($Δx$) and momentum uncertainty ($Δp$). As we increased the standard deviation of the initial wave packet in position space, i.e., the value of $Δx$, the measured standard deviation in momentum representation, $Δp$, decreased accordingly. This behavior stems from the fundamental property of the Fourier transform, which connects the position and momentum space representations.

![The uncertainty product as a function of the initial position's standard deviation.](uploaded:uncertainty_product.png)
*Figure 1: The value of the product $Δx·Δp$ as a function of the initial position's standard deviation, $σ_x$. The product is constant and approximately 0.5, which, in atomic units (where $ħ=1$), corresponds to the theoretical minimum of $ħ/2$.*

The most significant result is that the product of these two quantities, $Δx·Δp$, remained constant across the investigated range. Based on Figure 1 and the `heisenberg_scan.csv` data, the value of the product consistently hovers around $0.5$. In an atomic unit system ($ħ=1$), this corresponds precisely to the minimum value allowed by the Heisenberg relation, $Δx·Δp ≥ ħ/2$.

**2. Wave Packet Density Distributions**

To visualize the structure of the wave packets, Figures 2 and 3 show the probability density distributions in position and momentum space for a representative state with an initial standard deviation of $σ_x = 3.0$. Both distributions, as theoretically expected, have a Gaussian shape. A wider distribution in position space (larger $Δx$) results in a narrower distribution in momentum space (smaller $Δp$), visually confirming the inverse proportionality inherent in the uncertainty principle.

<center>
![Position probability density distribution.](uploaded:position_density_sigma3.0.png) 
</center>

*Figure 2: The probability density of the wave packet in position space ($|\psi(x)|^2$) for an initial state with a standard deviation of $σ_x = 3.0$.*

<center>
![Momentum probability density distribution.](uploaded:momentum_density_sigma3.0.png) 
</center>

*Figure 3: The probability density of the wave packet in momentum space ($|\phi(p)|^2$) for an initial state with a standard deviation of $σ_x = 3.0$.*

**3. Time Evolution of the Wave Packet**

The simulation was also extended to investigate the free time evolution of the wave packet (Figure 4). During this process, the position uncertainty, $Δx(t)$, monotonically increases over time. This phenomenon is known as "wave packet spreading" and arises because the different plane wave components making up the wave packet propagate at different velocities. It is important to note that since no external force acts on the particle, its momentum distribution—and thus its momentum uncertainty $Δp$—remains constant in time.

![Time evolution of the position uncertainty.](uploaded:free_spreading_sigma_x_t.png)
*Figure 4: The change in position uncertainty $Δx(t)$ over time for a freely propagating Gaussian wave packet. The wave packet spreads out in time, resulting in an increase in $Δx$.*

---

### Discussion

The presented numerical results confirm and illustrate fundamental concepts of quantum mechanics from several perspectives.

**Numerical Verification of the Heisenberg Uncertainty Principle**
The simulation clearly and quantitatively validates the Heisenberg uncertainty principle ($Δx·Δp ≥ ħ/2$). The data shows that an unavoidable, inverse relationship exists between the uncertainties of these two physical quantities. Their product never falls below the theoretical limit of $ħ/2$, which is an inherent property of quantum systems.

**The Gaussian Wave Packet as a Minimum Uncertainty State**
Our results highlight the special role that Gaussian wave packets play in quantum mechanics. The fact that their uncertainty product $Δx·Δp$ assumes the minimum possible value, $ħ/2$, means that these states are **minimum uncertainty wave packets**. In other words, a Gaussian wave packet describes the "most classical-like" state possible, where a particle's position and momentum are simultaneously defined with the highest possible precision.

**Time Evolution of Uncertainty**
The study of time evolution shows that although the position uncertainty ($Δx(t)$) increases during free evolution, the momentum uncertainty ($Δp$) remains constant. Consequently, the uncertainty product, $Δx(t)·Δp$, also increases over time. This is in perfect agreement with the Heisenberg relation, as the product continues to satisfy the inequality $Δx(t)·Δp ≥ ħ/2$; it simply moves away from the minimum value as time progresses.

**Educational and Illustrative Value**
Finally, this simulation possesses outstanding educational and demonstrative value. The underlying Python code (which forms the basis of the simulation) is easy to run and reproduce, allowing students and researchers to interactively explore one of the most important and least intuitive principles of quantum mechanics. The visual results (plots and density distributions) effectively aid in understanding these concepts, bridging the gap between abstract mathematical formalism and physical reality. ⚛️🔬
