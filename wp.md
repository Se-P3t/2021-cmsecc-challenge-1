$$
a_i = 2^k y_i + z_i \\
a_{i+j} = \sum_{l=0}^{n-1} q_{j,l} a_{i+l} \bmod m
$$


$$
\sum_{l=0}^{n-1} q_{j,l} z_l - z_j \equiv 2^k \left(y_j - \sum_{l=0}^{n-1} q_{j,l} y_l \right) \equiv c_j, (n \leq j < d)
$$

$$
\left(\begin{matrix}
-1 \\ z_0 \\ \vdots \\ z_{n-1} \\ k_n \\ \vdots \\ k_{d-1}
\end{matrix}\right)^T \cdot
\left(\begin{matrix}
1 & & & &c_n &\cdots &c_{d-1} \\
&1 & & &q_{n,0} &\cdots &q_{d-1,0} \\
& &\ddots & &\vdots &\cdots &\vdots \\
& & &1 &q_{n,n-1} &\cdots &q_{d-1,n-1} \\
& & & &m & & \\
& & & & &\ddots & \\
& & & & & &m \\
\end{matrix}\right) =
\left(-1, z_0, \ldots, z_{n-1}, z_n, \ldots, z_{d-1} \right)
$$



scale: $\left(m, m/z, \ldots, m/z, m/z, \ldots, m/z\right)$

dim: $d+1$

$||v||_2 \leq \sqrt{d+1} m $

$$
det(L) = \frac{m^{2d-n+1}}{z^{d+1}}
$$

$$
||b_1||_2 \leq 2^{d/4} z^{-1} m^{1+\frac{d-n}{d+1}}
$$





$Y_i = (y_i, y_{i+1}, \ldots, y_{i+t-1}), i = 0, \ldots, {r-1}, t < r$

there exists $\vec{\eta}$ s.t. $\sum_{i=0}^{r-1} \eta_i Y_i = 0$ (SIS)

define $A_i = (a_i, a_{i+1}, \ldots, a_{i+t-1})$, let $U = \sum_{i=0}^{r-1} \eta_i A_i$, we have

$$
U = \sum_{i=0}^{r-1} \eta_i A_i - 2^k \sum_{i=0}^{r-1} \eta_i Y_i = \sum_{i=0}^{r-1} \eta_i Z_i,
$$

where $Z_i = (z_i, z_{i+1}, \ldots, z_{i+t-1})$. (bounded by $2^k$)

$$
\left(\begin{matrix}
a_i \\ \vdots \\ a_{i+n-1} \\ k_{i+n} \\ \vdots \\ k_{i+t-1}
\end{matrix}\right)^T \cdot
\left(\begin{matrix}
&1 & & &q_{n,0} &\cdots &q_{t-1,0} \\
& &\ddots & &\vdots &\cdots &\vdots \\
& & &1 &q_{n,n-1} &\cdots &q_{t-1,n-1} \\
& & & &m & & \\
& & & & &\ddots & \\
& & & & & &m \\
\end{matrix}\right) =
\left(a_i, \ldots, a_{i+n-1}, a_{i+n}, \ldots, a_{i+t-1} \right) = 
A_i
$$

$A_i \in L \longrightarrow U \in L$

if $||U||_2 < \lambda_1(L)$, then $U = \vec{0} = \sum_{i=0}^{r-1} \eta_i A_i$

$||U|| < \sqrt{r} 2^k ||\vec\eta||$, $\lambda_1(L) \approx \sqrt{\frac{t}{2\pi e}}det(L)^{1/t} = O(m^{1-n/t})$





Let $g_i = \eta_i + \sum_{j=n}^{r-1} \eta_j q_{j,i}$, we can get

$$
\left(\begin{matrix}
a_0 &a_1 &\cdots &a_{n-1} \\
a_1 &a_2 &\cdots &a_n \\
\vdots &\vdots &\vdots &\vdots \\
a_{t-1} &a_t &\cdots &a_{t+n-2}
\end{matrix}\right) \cdot
\left(\begin{matrix}
g_0 \\ g_1 \\ \vdots \\ g_{n-1}
\end{matrix}\right) \equiv \vec{0} \pmod{m}
$$

Thus $g_i \bmod m = 0$


$$
L(g_i) = 
\left(\begin{matrix}
m &0 &0 &\cdots &0 \\
-q_{n,i} &1 &0 &\cdots &0 \\
-q_{n+1, i} &0 &1 &\cdots &0 \\
\vdots &\vdots &\vdots &\ddots &\vdots \\
-q_{r-1,i} &0 &0 &\cdots &1
\end{matrix}\right)
$$

$\vec\eta(i) = (\eta_i, \eta_n, \eta_{n+1}, \ldots, \eta_{r-1}) \in L(g_i)$













$$
\begin{aligned}
\left(\begin{matrix} 1 \\ q_{n,i} \\ \vdots \\ q_{r-1,i} \\ -u^{(0)} \\ \vdots \\ -u_{(d-1)} \end{matrix}\right)^T 
&\cdot \left(\begin{matrix}
1 & & & &\eta_i^{(0)} &\cdots &\eta_i^{(d-1)} \\
 &1 & & &\eta_n^{(0)} &\cdots &\eta_n^{(d-1)} \\
 & &\ddots & &\vdots &\cdots &\vdots \\
 & & &1 &\eta_{r-1}^{(0)} &\cdots &\eta_{r-1}^{(d-1)} \\
 & & & &m & & \\
 & & & & &\ddots & \\
 & & & & & &m\\
\end{matrix}\right) \\
&= \left(\begin{matrix} 1 &q_{n,i} &\cdots &q_{r-1,i} &0 &\ldots &0 \end{matrix}\right)
\end{aligned}
$$

**failed**: wanted vector too large      (could used as $\eta'$ ? )










