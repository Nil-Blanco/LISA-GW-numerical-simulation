# Mathematical Framework: Gravitational Wave Timing Residuals

## 1. Natural Units and Geometrical Constants
[cite_start]To simplify the computational models, all calculations are performed in natural units where the speed of light is $c=1$, converting all spatial dimensions to seconds[cite: 5].

* [cite_start]**LISA Arm Length ($L$):** $2.5 \times 10^6 \text{ km} \Rightarrow 8.34 \text{ s}$[cite: 6].
* [cite_start]**GW Frequency ($f$):** $100 \text{ mHz} \Rightarrow \omega = 0.628 \text{ s}^{-1}$[cite: 7, 8].
* [cite_start]**Hubble Constant ($H_0$):** $70 \text{ km/s/Mpc} \Rightarrow 2.27 \times 10^{-18} \text{ s}^{-1}$[cite: 9].
* **Source Distance ($Z_A$):** $1 \text{ Gpc} \Rightarrow 1.03 \times 10^{17} \text{ s}$[cite: 9].

## 2. Geometrical Model and Phase Derivation
[cite_start]Considering a triangle formed by the source ($S$) and two LISA satellites ($A$ and $B$) separated by distance $L$, where the distance from $A$ to the source is $Z_A$ and the incidence angle is $\alpha$[cite: 15, 16]. [cite_start]Assuming $Z_A \gg L$, the position along the LISA arm is approximated as[cite: 19]:

[cite_start]$$\vec{R}(x) \approx Z_A + xL\cos(\alpha)$$ [cite: 21]

[cite_start]The integral to calculate the timing residual $\tau_{GW}$ takes the form[cite: 23]:

[cite_start]$$\tau_{GW}(T_A) = -\frac{L\epsilon}{2}\sin^2(\alpha)\int_{-1}^{0}\frac{\epsilon_{ij}^{\prime}}{R(x)}(1+HT(x))\cos[\Theta(x)]dx$$ [cite: 23]

[cite_start]Through algebraic substitution, the phase $\Theta(x)$ becomes a quadratic function $\Theta(x) = Ax^2 + Bx + C$[cite: 28], where the coefficients are:
* $A = \frac{\omega HL^2}{2}(\cos\alpha - 2)\cos\alpha$ [cite: 30]
* [cite_start]$B = L\omega(-HZ_A + (-HT_A + HZ_A - 1)\cos\alpha + 1)$ [cite: 31]

## 3. Analytical Approach: Stationary Phase Approximation (SPA)
To analytically find the optimal angle that maximizes the response, we apply the Stationary Phase Approximation. [cite_start]The function reaches a maximum when the cosine argument is evaluated at $B=0$[cite: 36, 37]. 

[cite_start]Solving for $\alpha$, and noting that $Z_A \approx T_A$, we obtain the optimal angle expression[cite: 38, 41, 42]:

[cite_start]$$\alpha_{optim} = 2\arcsin\left(\sqrt{\frac{HZ_A}{2}}\right)$$ [cite: 42]

[cite_start]Using the defined cosmological values ($H = 2.3 \times 10^{-18} \text{ s}^{-1}$ and $Z_A = 1.0296 \times 10^{17} \text{ s}$), the analytical optimal angle is **$\alpha \approx 40.25^\circ$ ($0.702 \text{ rad}$)**[cite: 44, 46, 48].

## 4. Computational Method and Discrepancy Analysis
[cite_start]To validate the analytical model, the integral was resolved numerically[cite: 54]:

[cite_start]$$\tau_{GW}(T_A) = K \sin^2(\alpha) \int_{-1}^{0} g(x, \alpha) \cos(f(x, \alpha)) dx$$ [cite: 61]

[cite_start]An algorithm was developed to evaluate this integral using the trapezoidal rule for equidistant $\alpha$ values between $0^\circ$ and $90^\circ$[cite: 66, 67].

**Results & Conclusions:**
[cite_start]While the analytical SPA yielded $40.25^\circ$, the computational numerical integration located the absolute maximum at **$\alpha \approx 48.25^\circ$**[cite: 72]. 
This discrepancy is theoretically consistent. [cite_start]The analytical approach maximizes only the internal phase ($\Theta(x)$), whereas the numerical method accounts for the global $\sin^2(\alpha)$ pre-factor[cite: 73, 74]. [cite_start]Because $\sin^2(\alpha)$ naturally maximizes at $90^\circ$, it actively "pushes" the overall optimal angle higher than the purely phase-derived value[cite: 74].