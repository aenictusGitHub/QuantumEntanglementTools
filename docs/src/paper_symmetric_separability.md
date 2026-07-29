# Symmetric separability and SAPPT witnesses

This executable case study reproduces selected calculations from
[J. Louvet, E. Serrano-Ensástiga, T. Bastin, and J. Martin, *Phys. Rev. A*
**111**, 042418 (2025)](https://doi.org/10.1103/PhysRevA.111.042418).
It illustrates three logically different conclusions:

- an explicit convex decomposition certifies separability;
- a negative entanglement-witness expectation certifies entanglement; and
- the paper's spectral theorem establishes symmetric absolute positive partial
  transposition (SAPPT).

These conclusions are deliberately kept separate. In particular, PPT and
SAPPT do not by themselves imply separability for the five-qubit states below.

## The one-parameter family

The symmetric subspace of ``N`` qubits has the ordered Dicke basis

```math
\left\{
|D_N^{(0)}\rangle,\ldots,|D_N^{(N)}\rangle
\right\},
```

where ``|D_N^{(k)}\rangle`` is the normalized equal superposition of
computational-basis vectors with ``k`` excitations. The paper considers

```math
\rho(p)
=
p\rho_0+(1-p)|\psi_0\rangle\langle\psi_0|,
\qquad
\rho_0=\frac{I_{N+1}}{N+1},
\qquad 0\leq p\leq1.
```

Here ``I_{N+1}`` is the identity **inside the symmetric subspace**, not the
``2^N``-dimensional identity on the full qubit Hilbert space. The parameter
``p`` is the weight of this symmetric maximally mixed state.

For any normalized symmetric ``|\psi_0\rangle``, the spectrum is

```math
\mathrm{spec}\,\rho(p)
=
\left(
1-\frac{Np}{N+1},
\underbrace{\frac{p}{N+1},\ldots,\frac{p}{N+1}}_{N\ \mathrm{times}}
\right).
```

For the bipartition ``k|N-k``, the paper obtains the unitary-orbit lower
bound

```math
\min_{U\in SU(N+1)}
\lambda_{\min}\!\left[
    \left(U\rho(p)U^\dagger\right)^{T_A}
\right]
\geq
\frac{p}{(N+1)\binom Nk}-\frac{1-p}{2}.
```

The strictest condition is at ``k=\lfloor N/2\rfloor``. The resulting spectral
SAPPT threshold is

```math
p_{\min}
=
\left[
1+
\frac{2}{(N+1)\binom{N}{\lfloor N/2\rfloor}}
\right]^{-1}.
```

For five qubits,

```math
p_{\min}=\frac{30}{31}\approx0.9677419355.
```

The theorem concerns the whole symmetric-unitary orbit: every symmetric state
with this spectrum is SAPPT when ``p\geq p_{\min}``. A PPT computation on one
representative would not, by itself, prove that absolute statement.

## Same spectrum, different separability

Set ``N=5`` and ``p=p_{\min}``. Consider first

```math
\rho_{\mathrm{sep}}
=
p\frac{I_6}{6}
+(1-p)|D_5^{(0)}\rangle\langle D_5^{(0)}|.
```

Since ``|D_5^{(0)}\rangle=|0\rangle^{\otimes5}``, its second term is a product
projector. The first term also has an explicit product decomposition. To make
that statement executable, the tutorial uses three Gauss--Legendre nodes

```math
z_i\in\left\{-\sqrt{\frac35},0,\sqrt{\frac35}\right\},
\qquad
\omega_i\in\left\{\frac5{18},\frac49,\frac5{18}\right\},
```

and six phases ``\phi_j=2\pi j/6`` for ``j=0,\ldots,5``. Define the normalized
single-qubit vectors

```math
|\varphi_{ij}\rangle
=
\sqrt{\frac{1+z_i}{2}}\,|0\rangle
+
\sqrt{\frac{1-z_i}{2}}\,e^{i\phi_j}|1\rangle.
```

The corresponding symmetric products have Dicke coefficients

```math
\langle D_5^{(k)}|\varphi_{ij}\rangle^{\otimes5}
=
\sqrt{\binom5k}
\left(\frac{1+z_i}{2}\right)^{(5-k)/2}
\left(\frac{1-z_i}{2}\right)^{k/2}
e^{ik\phi_j}.
```

The six phases cancel every off-diagonal Dicke-basis entry. Three-point
Gauss--Legendre quadrature integrates the remaining degree-five diagonal
polynomials exactly, giving

```math
\frac{I_6}{6}
=
\sum_{i=1}^{3}\sum_{j=0}^{5}
\frac{\omega_i}{6}
|\varphi_{ij}\rangle^{\otimes5}
\langle\varphi_{ij}|^{\otimes5}.
```

Thus ``\rho_{\mathrm{sep}}`` has a concrete 19-term decomposition:

- 18 spin-coherent product projectors with weights
  ``p\omega_i/6``; and
- ``|0\rangle^{\otimes5}\langle0|^{\otimes5}`` with weight ``1-p``.

All weights are nonnegative and sum to one. The tutorial constructs these 19
terms and checks their normalization and reconstruction residual directly.

Now replace the distinguished product vector by

```math
|\mathrm{GHZ}_5^+\rangle
=
\frac{|D_5^{(0)}\rangle+|D_5^{(5)}\rangle}{\sqrt2}
```

to obtain

```math
\rho_{\mathrm{GHZ}}(p)
=
p\frac{I_6}{6}
+(1-p)|\mathrm{GHZ}_5^+\rangle\langle\mathrm{GHZ}_5^+|.
```

The two density matrices have the same spectrum and are connected by a
unitary acting within the symmetric subspace. Nevertheless, the witness below
certifies that ``\rho_{\mathrm{GHZ}}(p_{\min})`` is entangled. Consequently,
the explicitly separable ``\rho_{\mathrm{sep}}`` is SAPPT but not symmetric
absolutely separable (SAS): another state on its allowed unitary orbit is
entangled.

## Reconstructing the published five-qubit witness

In the ordered Dicke basis, the rounded witness published in the paper is

```math
W_5=
\begin{pmatrix}
a&0&0&0&0&c\\
0&b&0&0&0&0\\
0&0&1&0&0&0\\
0&0&0&1&0&0\\
0&0&0&0&b&0\\
c&0&0&0&0&a
\end{pmatrix},
\qquad
\begin{aligned}
a&=0.0366656,\\
b&=-0.134595,\\
c&=-9.31947.
\end{aligned}
```

This is a **symmetric** entanglement witness. Every separable symmetric state
is a convex mixture of vectors ``|\varphi\rangle^{\otimes5}``, so witness
validity reduces to checking

```math
\langle\varphi|^{\otimes5}W_5|\varphi\rangle^{\otimes5}\geq0
```

for every normalized single-qubit ``|\varphi\rangle``.

The matrix construction itself is direct Julia:

```julia
using LinearAlgebra

a, b, c = 0.0366656, -0.134595, -9.31947
W5 = Matrix(Diagonal(ComplexF64[a, b, 1, 1, b, a]))
W5[1, end] = c
W5[end, 1] = c
```

The tutorial's `published_symmetric_witness(5)` performs this construction
and also accepts a GHZ `phase` keyword.

To verify block positivity, write

```math
|\varphi\rangle
=x|0\rangle+y e^{i\phi}|1\rangle,
\qquad
x=\cos\frac{\theta}{2},
\quad
y=\sin\frac{\theta}{2}.
```

Then

```math
\begin{aligned}
F_5(\theta,\phi)
={}&a(x^{10}+y^{10})
+5b(x^8y^2+x^2y^8)\\
&+10(x^6y^4+x^4y^6)
+2c x^5y^5\cos(5\phi).
\end{aligned}
```

Because ``c<0``, the phase minimum has ``\cos(5\phi)=1``. With
``t=xy\in[0,1/2]``, the remaining expression is

```math
f(t)
=
a+5(b-a)t^2+(5a-15b+10)t^4+2ct^5.
```

The tutorial treats the displayed decimal coefficients as exact rationals and
proves positivity without a floating grid or root finder. Set ``u=2t``. Exact
polynomial division gives

```math
f(u/2)-f(1/2)=(1-u)q(u),
```

where

```math
q(u)=
\frac{
542429+542429u-2882783u^2-2882783u^3+9319470u^4
}{16000000}.
```

The code converts ``q`` to the Bernstein basis and raises its degree from four
to eight. Its exact degree-eight Bernstein coefficients are

```math
\left(
\frac{542429}{16000000},
\frac{4881861}{128000000},
\frac{2012779}{56000000},
\frac{96373}{4000000},
\frac{246253}{32000000},
\frac{155429}{896000000},
\frac{5213133}{224000000},
\frac{6851851}{64000000},
\frac{2319381}{8000000}
\right).
```

Every coefficient is positive, and every Bernstein basis polynomial is
nonnegative on ``u\in[0,1]``. Hence ``q(u)>0`` throughout the interval and the
global minimum occurs exactly at ``t=1/2``, corresponding to
``(\theta,\phi)=(\pi/2,0)``:

```math
\min_{\theta,\phi}F_5(\theta,\phi)
=
\frac{a+5b+10+c}{16}
=\frac{221103}{80000000}
=0.0027637875>0.
```

This is an exact positivity certificate for the rounded decimal matrix as
printed. It does not recover any unprinted higher-precision SDP coefficients.

For ``|\mathrm{GHZ}_5^+\rangle``, define

```math
m=\mathrm{Tr}\left(W_5\frac{I_6}{6}\right)
=\frac{a+b+1}{3},
\qquad
g=\langle\mathrm{GHZ}_5^+|W_5|\mathrm{GHZ}_5^+\rangle
=a+c.
```

The witness expectation is therefore the affine function

```math
\mathrm{Tr}\!\left[W_5\rho_{\mathrm{GHZ}}(p)\right]
=pm+(1-p)g.
```

It is negative below the rounded-coefficient endpoint

```math
p_{W_5}
=-\frac{g}{m-g}
\approx0.9686241593.
```

At ``p=p_{\min}``, the value is approximately ``-0.0084547871``. A zero
expectation is inconclusive, so executable checks use a point strictly inside
the negative interval rather than treating the rounded endpoint as a
certificate.

## A strict SAPPT bound-entangled point

Choose

```math
p=\frac{121}{125}=0.968.
```

This satisfies

```math
\frac{30}{31}<\frac{121}{125}<p_{W_5}.
```

The paper's theorem therefore makes the state SAPPT, while the reconstructed
witness gives

```math
\mathrm{Tr}\!\left[
W_5\rho_{\mathrm{GHZ}}\!\left(\frac{121}{125}\right)
\right]
\approx-0.0059816272<0.
```

The tutorial also evaluates the representative state's partial transposes
within the supported products of the two subsystems' symmetric sectors. It
finds restricted-sector minimum eigenvalues approximately ``0.0162667`` for
``1|4`` and ``0.000133333`` for ``2|3``. The full computational-space
partial transposes also contain zero modes orthogonal to those supported
sectors. The positive restricted values are useful consistency checks; the
SAPPT conclusion itself comes from the spectral theorem, not from checking
only these two matrices. The state is therefore a strict, non-boundary SAPPT
bound-entangled example.

## Below the SAPPT threshold: a decomposable NPT witness

For the same GHZ family below ``p_{\min}``, the ``2|3`` partial transpose has
the eigenvalue

```math
\lambda_{\min}
=
\frac{31p-30}{60}.
```

At ``p=29/30`` this is exactly

```math
\lambda_{\min}=-\frac1{1800}.
```

An associated normalized eigenvector is

```math
|\eta\rangle
=
\frac{
|D_2^{(0)}\rangle|D_3^{(3)}\rangle
-
|D_2^{(2)}\rangle|D_3^{(0)}\rangle
}{\sqrt2}.
```

Set ``Q=|\eta\rangle\langle\eta|\succeq0`` and construct

```math
W_{\mathrm{NPT}}=Q^{T_A}.
```

In the full five-qubit computational basis, that construction is:

```julia
using QuantumEntanglementTools

eta = (
    tensor_product(dicke_state(2, 0), dicke_state(3, 3)) -
    tensor_product(dicke_state(2, 2), dicke_state(3, 0))
) / sqrt(2)
Q = eta * eta'
W_npt = partial_transpose(Q, ntuple(_ -> 2, 5); systems=(1, 2))
```

This is a decomposable witness: for every state ``\sigma`` that is PPT across
the ``2|3`` split,

```math
\mathrm{Tr}(W_{\mathrm{NPT}}\sigma)
=
\mathrm{Tr}(Q\sigma^{T_A})\geq0.
```

For the target state,

```math
\mathrm{Tr}\!\left[
W_{\mathrm{NPT}}\rho_{\mathrm{GHZ}}\!\left(\frac{29}{30}\right)
\right]
=
\langle\eta|\rho_{\mathrm{GHZ}}^{T_A}|\eta\rangle
=-\frac1{1800}.
```

Unlike ``W_5``, which is represented as a ``6\times6`` symmetric-subspace
witness, the tutorial constructs ``W_{\mathrm{NPT}}`` in the full
``2^5\times2^5`` computational space with [`dicke_state`](@ref),
[`tensor_product`](@ref), and [`partial_transpose`](@ref).

## Matching the GHZ phase convention

The journal PDF writes the GHZ vector with a minus sign,

```math
|\mathrm{GHZ}_5^-\rangle
=
\frac{|D_5^{(0)}\rangle-|D_5^{(5)}\rangle}{\sqrt2},
```

while its displayed negative corner ``c`` is phase-matched to
``|\mathrm{GHZ}_5^+\rangle``. Taken together without adjustment, the printed
minus vector and negative corner do not reproduce the paper's claimed
negative expectation. Indeed,

```math
\langle\mathrm{GHZ}_5^-|W_5|\mathrm{GHZ}_5^-\rangle=a-c>0,
```

whereas the plus convention gives ``a+c<0``.

For the mixed state at ``p_{\min}``, combining the written minus vector with
the unadjusted negative-corner matrix gives approximately ``+0.59280134``.
The executable regression checks that this value is positive before applying
the phase correction.

The tutorial uses ``|\mathrm{GHZ}_5^+\rangle`` with the displayed negative
corner. It also demonstrates the equivalent phase-matched construction. For

```math
|\mathrm{GHZ}_5(\delta)\rangle
=
\frac{|D_5^{(0)}\rangle+e^{i\delta}|D_5^{(5)}\rangle}{\sqrt2},
```

conjugating the witness by the same endpoint phase changes its upper-right
corner to ``ce^{-i\delta}`` and its lower-left corner to
``ce^{i\delta}``. At ``\delta=\pi``, this flips the real corner sign and
restores exactly the same expectation for the minus-phase state. This is a
phase convention, not a different entanglement result.

## Published rounded witnesses for five, seven, and nine qubits

The tutorial helper constructs all three reported witnesses as

```math
W_N
=
\mathrm{Diag}(d_0,\ldots,d_N)
+z\left(
|D_N^{(0)}\rangle\langle D_N^{(N)}|
+|D_N^{(N)}\rangle\langle D_N^{(0)}|
\right).
```

The rounded coefficients are:

| ``N`` | Dicke-basis diagonal ``(d_0,\ldots,d_N)`` | corner ``z`` | extension split used in the paper |
|---:|---|---:|---:|
| 5 | ``(0.0366656, -0.134595, 1, 1, -0.134595, 0.0366656)`` | ``-9.31947`` | ``1|4`` |
| 7 | ``(0.00197514, 0.0643064, -0.189017, 1, 1, -0.189017, 0.0643064, 0.00197514)`` | ``-31.2405`` | ``1|6`` |
| 9 | ``(0.00235791, -0.013747, 0.0621661, -0.1636915, 1, 1, -0.1636915, 0.0621661, -0.013747, 0.00235791)`` | ``-114.305`` | ``4|5`` |

Recomputing expectations and detection endpoints from those rounded numbers
gives:

| ``N`` | ``p_{\min}`` | ``\mathrm{Tr}[W_N\rho(p_{\min})]`` | rounded-coefficient endpoint |
|---:|---:|---:|---:|
| 5 | ``30/31`` | ``-0.0084547871`` | ``0.9686241593`` |
| 7 | ``140/141`` | ``-0.0037891203`` | ``0.9930282522`` |
| 9 | ``630/631`` | ``-0.0040092993`` | ``0.9984502358`` |

The paper reports minimum symmetric-product expectations of approximately
``0.00276`` for ``W_5``, ``0.001975`` for ``W_7``, and ``0.0002234`` for
``W_9``. The extra digits in the table above are arithmetic consequences of
the published rounded coefficients, not recovered higher-precision witness
data.

## What the executable example certifies

The provenance of each conclusion matters:

- The 19-term positive decomposition is an explicit separability certificate
  for the ``|D_5^{(0)}\rangle`` representative.
- The product-state minimization and negative ``W_5`` expectation reproduce an
  explicit symmetric entanglement-witness certificate.
- ``W_{\mathrm{NPT}}=Q^{T_A}`` is constructed directly and certifies the
  below-threshold NPT example.
- The SAPPT conclusion uses the paper's analytic spectral theorem; the package
  does not currently expose a general SAPPT predicate.
- The published ``W_5``, ``W_7``, and ``W_9`` arose from duals of failed
  two-copy PPT symmetric-extension problems. The tutorial reconstructs and
  checks the reported matrices; it does **not** rerun those witness-generation
  semidefinite programs.
- The paper reports, using a truncated-moment semidefinite method with
  ``10^{-5}`` parameter resolution, that the five-qubit GHZ family is
  entangled for ``p_{\min}\leq p< p_{\mathrm{ent}}=0.96953`` and separable
  for ``p_{\mathrm{ent}}\leq p\leq1``. The tutorial does not rerun that
  optimization, and ``0.96953`` is not presented as a package-owned exact
  separability certificate.

Outside the interval detected by a particular witness, its expectation is
merely inconclusive. It must not be inverted into a separability claim.

## Run the complete example

The repository tutorial performs all constructions above, checks their
invariants, and prints values from the live calculation:

```@example tutorial-symmetric-sappt-witnesses
using QuantumEntanglementTools

path = joinpath(
    pkgdir(QuantumEntanglementTools),
    "tutorials",
    "symmetric_sappt_witnesses.jl",
)
include(path);
TutorialSymmetricSAPPTWitnesses.run()
```
