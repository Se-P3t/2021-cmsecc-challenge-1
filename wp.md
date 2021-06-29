for the details on the results, refer to [challenge.md](./challenge.md)

# challenge 1

## category 1

### reading list

- [**The General Sieve Kernel and New Records in Lattice Reduction**](https://eprint.iacr.org/2019/089)
- [**On Bounded Distance Decoding with Predicate: Breaking the "Lattice Barrier" for the Hidden Number Problem**](https://eprint.iacr.org/2020/1540) <https://github.com/malb/bdd-predicate>
- [**Advanced Lattice Sieving on GPUs, with Tensor Cores**](https://eprint.iacr.org/2021/141) <https://github.com/WvanWoerden/G6K-GPU-Tensor>
- [**Guessing Bits: Improved Lattice Attacks on (EC)DSA**](https://eprint.iacr.org/2021/455)
- [**Lattice Enumeration on GPUs for fplll**](https://eprint.iacr.org/2021/430) <https://github.com/FeanorTheElf/fplll-CUDA-enumeration>
- [**Dual lattice attacks for closest vector problems (with preprocessing)**](https://eprint.iacr.org/2021/557)

### [recover_initial_state.py](./recover_initial_state.py)

TODO

### [recover_initial_state__embedding.py](./recover_initial_state__embedding.py)

$$
a_i = 2^k y_i + z_i \\
a_{i+j} = \sum_{l=0}^{n-1} q_{j,l} a_{i+l} \bmod m
$$

$$
\sum_{l=0}^{n-1} q_{j,l} z_l - z_j \equiv 2^k \left(y_j - \sum_{l=0}^{n-1} q_{j,l} y_l \right) \equiv c_j, (n \leq j < d)
$$

$$
\begin{aligned}
\begin{pmatrix}
-1 \\ z_0 \\ \vdots \\ z_{n-1} \\ k_n \\ \vdots \\ k_{d-1}
\end{pmatrix}^T
&\cdot \begin{pmatrix}
1 & & & &c_n &\cdots &c_{d-1} \\
&1 & & &q_{n,0} &\cdots &q_{d-1,0} \\
& &\ddots & &\vdots &\cdots &\vdots \\
& & &1 &q_{n,n-1} &\cdots &q_{d-1,n-1} \\
& & & &m & & \\
& & & & &\ddots & \\
& & & & & &m \\
\end{pmatrix} \\
&= \begin{pmatrix} -1 &z_0 &\ldots &z_{n-1} &z_n &\ldots &z_{d-1}\end{pmatrix}
\end{aligned}
$$

bound of $z_i$: $z = 2^k > z_i \geq 0$

scale: $\left(m, m/z, \ldots, m/z, m/z, \ldots, m/z\right)$

dim: $d+1$

norm of target vector:

$$
||v||_2 = \sqrt{m^2 + \sum_{i=0}^{d-1} (z_i \cdot \frac{m}{z})^2} < \sqrt{d+1} m
$$

determinant of basis matrix:

$$
det(L) = \frac{m^{2d-n+1}}{z^d}
$$

expected $\lambda_1$ according to the Gaussian heuristic:

$$
gh(L) \approx \sqrt{\frac{dim}{2\pi e}} det(L)^{1/dim} = \sqrt{\frac{d+1}{2\pi e}} z^{-d/(d+1)} m^{1+\frac{d-n}{d+1}}
$$

norm of the first LLL-reduced ($\delta = 3/4$) vector:

$$
||b_1||_2 \leq 2^{(dim-1)/4} det(L)^{1/dim} = 2^{d/4} z^{-d/(d+1)} m^{1+\frac{d-n}{d+1}}
$$

expected norm of the shortest vector found by BKZ-$\beta$:

$$
||b_1||_2 \approx \delta_\beta^{dim-1} \cdot det(L)^{1/dim},
$$
​where $\delta_\beta = \sqrt{\beta/(2\pi e)}^{1/(\beta-1)}$.

firstly, the target vector must be the shortest one: $||v||_2 < gh(L)$

for category 1, we have

$$
d > \frac{498.0471}{29.9529-zbits}
$$

| zbits (mbits=31) | 2  | 3  | 14 | 21 | 24 | 25  | 26  | 27  | 28  | 29  |
| :--------------: | -- | -- | -- | -- | -- | --- | --- | --- | --- | --- |
| lower bound of d | 18 | 19 | 32 | 56 | 84 | 101 | 126 | 169 | 256 | 523 |

#### result

solve level 1-6 in minutes

## category 2

### [recover_coefficients.py](./recover_coefficients.py)

TODO: coding

recover the coefficients with given `ETA`

in this category, we're given exactly $r+t-1$ $y_i$'s, which is impossible to reach full rank: in the kernel lattice, we can get at most $r-t$ rows with leading zeros, even if they all satisfy $U_i = \sum_i \eta_i A_i = 0$, we still need more to make $L(g_i)^* = L(g_i)$.

Thus this script is not to be used in this category.

### [recover_coefficients__resultant.sage](./recover_coefficients__resultant.sage)

$g_i$ are polynomials, with degree at most $r-n$, sharing the root $(c_0, \ldots, c_{n-1})$ over $\mathbb{Z}/m\mathbb{Z}$

try to solve with resultant

#### result

we need at least two vectors such that $U_i = 0$

however, the degree increases rapidly, and I could only solve the first two level on my laptop

### [recover_coefficients__kernel.sage](./recover_coefficients__kernel.sage)

after we reduced the kernel lattice, we can get $d$ pairs of $\eta_i$ such that $\sum_i^{r-1} \eta_i Z_{j+i} = 0$ ($0 \leq j < t$)

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

namely, the $N$ dimensional vector $\vec z$ belongs to the left kernel of the matrix $ETA$, denoted as $B$

if we have enough pairs of $\eta_i$, the rank of $B$ will reduce to $2$

noticed that we also have $\vec y$ in the null space, generally, $L(B)$ is exactly lattice generated by $[\vec y, \vec z]$ ,so we might find $\vec z$ using LLL if $|\vec z| \ll |\vec y|$

it is easy to get the coefficients after we find the initial state:

$A_{i+1} \equiv A_i \cdot Q \pmod m$

just solve the equation over $\mathbb{Z}/m\mathbb{Z}$

#### note

since both $\vec y$ and $\vec z$ belong to the linear span, we have $gcd(B[0, j], B[1, j])$ divides $gcd(\vec y[j], \vec z_[j])$

to leak more information of $\vec z$, we can substitute these variable; for example, let $\vec y' := 2 \vec y + 1$, then $\vec z' = z' - 2^{zbits-1}$, relatively, we need to adjust the parameter $r$ and $t$

it is easy to see that this method will quickly enlarge the dimension of the lattice, resulting in low efficiency

#### problem

what can we learn further from the kernel when $|\vec z|$ is big

- see the next section
- ...

#### improvement

when $|\vec y| \ll |\vec z|$, the reduced basis matrix is likely to be $[\vec y, \pm\vec z + k_1 \vec y]$

similarly, we could get $[2\vec y +1, (\vec z - 2^{zbits-1})+k_2(2\vec y +1)]$ if $|\vec y'| \ll |\vec z'|$

FIND the left kernel of $[\vec y, \vec z +k_1\vec y, \vec z + k_2(2\vec y +1), \vec 1]$

for the sign of each vector, we *simply iterate all possible situations*

moreover, this method cannot be applied to category 1, since $||\vec y||$ is very small ($2$ ~ $3$ bit), if we only use the info. of $y_i$, the parameter $r$, $t$ will become super large

#### problem 2

we did not use `m`, can we do better

#### result

solve all levels of challenges with fewer BKZ calls (1 or 2 times)

## category 3

### [recover_modulus.py](./recover_modulus.py)

TODO: code

#### improvement

set $Y_i' := 2Y_i + 1$, then $U = \sum_i \eta_i A_i = 2^{k-1} \sum_i \eta_i Y_i' + \sum_i \eta_i (Z_i - 2^{k-1})$. since $||Z_i-2^{k-1}|| < \sqrt{rt} 2^{k-1}$

TODO: compare the performance

#### problem

the bound for $||\eta||$ is too ambiguous, we need to choose $r$ and $t$ manually

$$
\begin{aligned}
|u_i| = \sum_{j=0}^{r-1} \eta_j z_{i+j} &= \sum_{j=0}^{r-1} \eta_j (z_{i+j}-2^{zbits-1}) + 2^{zbits-1} \sum_{j=0}^{r-1} \eta_j \\
&< 2^{zbits-1}(\sqrt{r} ||\vec \eta|| + \sum_{j=0}^{r-1} \eta_j)
\end{aligned}
$$

### [recover_modulus__kernel.sage](./recover_modulus__kernel.sage)

we restored the initial state by the method mentioned in the second category;

since $A_{i+1} \equiv A_i \cdot Q \pmod m$, we have $A_{i+1} \cdot A_i^{-1} \equiv A_{j+1} \cdot A_j^{-1} \pmod m$

so we can compute GCD several times and finally get the modulus

then the rest are the same
