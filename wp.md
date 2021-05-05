# challenge 1

## category 1

## category 2

### [recover_coefficients__kernel.sage](./recover_coefficients__kernel.sage)

after we reduced the kernel lattice, we can get $d$ pairs of $\eta_i$ such that $\sum _i^{r-1} \eta_i Z_{j+i} = 0$ ($0 \leq j < t$)

combine them in one matrix, given that
$$
\begin{aligned}
\left(\begin{matrix}
z_0 \\ \vdots \\ z_{r-1} \\ \vdots \\ z_{N-1}
\end{matrix}\right)^T
&\cdot \left(\begin{matrix}
\eta_0^{(0)} &\cdots &\eta_0^{(d-1)} & & & & \\
\vdots &\cdots &\vdots &\ddots & & & \\
\eta_{r-1}^{(0)} &\cdots &\eta_{r-1}^{(d-1)} &\ddots &\eta_0^{(0)} &\cdots &\eta_0^{(d-1)} \\
 & & &\ddots &\vdots &\cdots &\vdots \\
 & & & &\eta_{r-1}^{(0)} &\cdots &\eta_{r-1}^{(d-1)}
\end{matrix}\right) \\
&= \left(\begin{matrix} 0 &\cdots &0 &\cdots &0 &\cdots &0 \end{matrix}\right)
\end{aligned}
$$

where $N = r+t-1$

namely, the $N$ dimensional vector $\vec z$ belongs to the left kernel of the matrix $ETA$

noticed that we also have $\vec y$ in the null space, so we might find $\vec z$ using LLL if $|\vec z| \ll |\vec y|$

#### coefficients

$A_{i+1} \equiv A_i \cdot Q \pmod m$

just solve the equation over $\mathbb{Z}/m\mathbb{Z}$

#### note

if we have enough pairs of $\eta_i$, the rank of $ker(ETA)$, denoted as $B$, will reduce to $2$

since both $\vec y$ and $\vec z$ belong to the linear span, we have $gcd(B[0, j], B[1, j])$ divides $gcd(\vec y[j], \vec z_[j])$

to leak more information of $\vec z$, we can substitute these variable; for example, let $\vec y' := 2 \vec y + 1$, then $\vec z' = z' - 2^{zbits-1}$, relatively, we need to adjust the parameter $r$ and $t$

it is easy to see that this method will quickly enlarge the dimension of the lattice, resulting in low efficiency

OPEN PROBLEM: what can we learn further from the kernel when $|\vec z|$ is big



## category 3





# misc

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

there exists $\vec{\eta}$ s.t. $\sum_{i=0}^{r-1} \eta_i Y_i = 0$ (kernel lattice)

define $A_i = (a_i, a_{i+1}, \ldots, a_{i+t-1})$, let $U = \sum_{i=0}^{r-1} \eta_i A_i$, we have

$$
U = \sum_{i=0}^{r-1} \eta_i A_i - 2^k \sum_{i=0}^{r-1} \eta_i Y_i = \sum_{i=0}^{r-1} \eta_i Z_i,
$$

where $Z_i = (z_i, z_{i+1}, \ldots, z_{i+t-1})$. (bounded by $2^k$)

$$
\begin{aligned}
\left(\begin{matrix}
a_i \\ \vdots \\ a_{i+n-1} \\ k_{i+n} \\ \vdots \\ k_{i+t-1}
\end{matrix}\right)^T
&\cdot \left(\begin{matrix}
&1 & & &q_{n,0} &\cdots &q_{t-1,0} \\
& &\ddots & &\vdots &\cdots &\vdots \\
& & &1 &q_{n,n-1} &\cdots &q_{t-1,n-1} \\
& & & &m & & \\
& & & & &\ddots & \\
& & & & & &m \\
\end{matrix}\right) \\
&= \left(\begin{matrix} a_i &\ldots &a_{i+n-1} &a_{i+n} &\ldots &a_{i+t-1} \end{matrix}\right) \\
&= A_i
\end{aligned}
$$

$A_i \in L \longrightarrow U \in L$

if $||U||_2 < \lambda_1(L)$, then $U = \vec{0} = \sum_{i=0}^{r-1} \eta_i A_i$

$||U|| < \sqrt{rt} 2^k ||\vec\eta||$, $\lambda_1(L) \approx \sqrt{\frac{t}{2\pi e}}det(L)^{1/t} = O(m^{1-n/t})$



improve: set $Y_i' := 2Y_i + 1$, then $U = \sum_i \eta_i A_i = 2^{k-1} \sum_i \eta_i Y_i' + \sum_i \eta_i (Z_i - 2^{k-1})$. since $||Z_i-2^{k-1}|| < \sqrt{rt} 2^{k-1}$



**problem**: in category 2, we're given exactly $r+t-1$ $y_i$'s, which is impossible to reach full rank: in the kernel lattice, we can get at most $r-t$ rows with leading zeros, even if they all satisfy $U_i = \sum_i \eta_i A_i = 0$, we still need more to make $L(g_i)^* = L(g_i)$. Moreover, we did not use $m$.

**problem**: the bound for $||\eta||$ is too ambiguous.

we need NEW idea



$U = 0 \longrightarrow \sum \eta A = 0$

$vector(z_0, \ldots, z_{d-1})$ is left kernel of ${ETA}^T$











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









