# Mathematical Framework: Gravitational Wave Timing Residuals

## 1. Natural Units and Geometrical Constants
To simplify the computational models, all calculations are performed in natural units where the speed of light is $c=1$, converting all spatial dimensions to seconds.

* **LISA Arm Length ($L$):** $2.5 \times 10^6 \text{ km} \Rightarrow 8.34 \text{ s}$.
* **GW Frequency ($f$):** $100 \text{ mHz} \Rightarrow \omega = 0.628 \text{ s}^{-1}$.
* **Hubble Constant ($H_0$):** $70 \text{ km/s/Mpc} \Rightarrow 2.27 \times 10^{-18} \text{ s}^{-1}$.
* **Source Distance ($Z_A$):** $1 \text{ Gpc} \Rightarrow 1.03 \times 10^{17} \text{ s}$.

## 2. Geometrical Model and Phase Derivation
Considering a triangle formed by the source ($S$) and two LISA satellites ($A$ and $B$) separated by distance $L$, where the distance from $A$ to the source is $Z_A$ and the incidence angle is $\alpha$. Assuming $Z_A \gg L$, the position along the LISA arm is approximated as:

$$
\vec{R}(x) \approx Z_A + xL\cos(\alpha)
$$

The integral to calculate the timing residual $\tau_{GW}$ takes the form:

$$
\tau_{GW}(T_A) = -\frac{L\epsilon}{2}\sin^2(\alpha)\int_{-1}^{0}\frac{\epsilon_{ij}^{\prime}}{R(x)}(1+HT(x))\cos[\Theta(x)]dx
$$

Through algebraic substitution, the phase $\Theta(x)$ becomes a quadratic function $\Theta(x) = Ax^2 + Bx + C$, where the coefficients are:
* $A = \frac{\omega HL^2}{2}(\cos\alpha - 2)\cos\alpha$
* $B = L\omega(-HZ_A + (-HT_A + HZ_A - 1)\cos\alpha + 1)$

## 3. Analytical Approach: Stationary Phase Approximation (SPA)
To analytically find the optimal angle that maximizes the response, we apply the Stationary Phase Approximation. The function reaches a maximum when the cosine argument is evaluated at $B=0$. 

Solving for $\alpha$, and noting that $Z_A \approx T_A$, we obtain the optimal angle expression:

$$
\alpha_{optim} = 2\arcsin\left(\sqrt{\frac{HZ_A}{2}}\right)
$$

Using the defined cosmological values ($H = 2.3 \times 10^{-18} \text{ s}^{-1}$ and $Z_A = 1.0296 \times 10^{17} \text{ s}$), the analytical optimal angle is **$\alpha \approx 40.25^\circ$ ($0.702 \text{ rad}$)**.

## 4. Computational Method and Discrepancy Analysis
To validate the analytical model, the integral was resolved numerically:

$$
\tau_{GW}(T_A) = K \sin^2(\alpha) \int_{-1}^{0} g(x, \alpha) \cos(f(x, \alpha)) dx
$$

An algorithm was developed to evaluate this integral using the trapezoidal rule for equidistant $\alpha$ values between $0^\circ$ and $90^\circ$.

**Results & Conclusions:**
While the analytical SPA yielded $40.25^\circ$, the computational numerical integration located the absolute maximum at **$\alpha \approx 48.25^\circ$**. 
This discrepancy is theoretically consistent. The analytical approach maximizes only the internal phase ($\Theta(x)$), whereas the numerical method accounts for the global $\sin^2(\alpha)$ pre-factor. Because $\sin^2(\alpha)$ naturally maximizes at $90^\circ$, it actively "pushes" the overall optimal angle higher than the purely phase-derived value.
