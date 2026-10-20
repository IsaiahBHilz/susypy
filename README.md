<h1 align="center">
<img src="https://raw.githubusercontent.com/IsaiahBHilz/susypy/main/branding/logo/logo.png" width="150">
<br> SusyPy
</h1><br>

**SusyPy** (short for **Su**per**sy**mmetry **Py**thon) is a symbolic algebra system, built on top of [Cadabra2](https://cadabra.science/), for performing intricate tensor arithmetic in supersymmetry. It is described in

> S. J. Gates, Jr., I. B. Hilsenrath, and S. Hilsenrath, *Modern Tensor-Spinor Symbolic Algebra Algorithms and Computing Non-Closure Geometry & Holoraumy in 11D, N = 1 Supergravity*, [arXiv:2212.00614](https://arxiv.org/abs/2212.00614).

## **Features:**
* LaTeX-like inputs for tensor expressions.
* Complete compatibility with Cadabra2's existing library, `Ex`/`ExNode` objects, and tree structure.
* Complete simplification of tensor expressions in a canonical way.
* Spinor-index canonicalization.
* Two- and four-index Fierz expansions.
* Symbolic Fourier and inverse Fourier transforms.
* Evaluation of tensor traces.
* Handling of Feynman propagators.
* Supersymmetry multiplet solvers.
* SUSY-invariant action solver.
* Calculation of a multiplet's holoraumy.


# Installation:

SusyPy requires **Python 3.9 or newer** and **Cadabra2 2.6 or newer** (with its Python module).

### Linux

Cadabra2 publishes Linux wheels on PyPI, so a single command installs everything:

```bash
pip install susypy
```

### macOS and Windows

PyPI has no Cadabra2 wheels for macOS or Windows, so install Cadabra2 first and then SusyPy:

1. Install Cadabra2 (2.6 or newer) and make sure its Python module imports in the Python you will use:
   ```bash
   python -c "import cadabra2"
   ```
   See the [Cadabra2 download page](https://cadabra.science/download.html) and its build-from-source instructions for [macOS](https://github.com/kpeeters/cadabra2#macos) and [Windows](https://github.com/kpeeters/cadabra2#windows). You only need the Python module, so you can configure the build with `-DENABLE_FRONTEND=OFF -DENABLE_PY_JUPYTER=OFF`.
2. Install SusyPy into that same Python:
   ```bash
   pip install susypy
   ```

> NOTE (macOS with a conda/micromamba Python): conda's Python is statically linked, and Cadabra2 2.6's build links its Python module against `libpython`, which makes `import cadabra2` crash with a segmentation fault. Fixing it means changing Cadabra2's `core/CMakeLists.txt` to build the module with `add_library(cadabra2 SHARED ...)` plus `target_link_libraries(cadabra2 PRIVATE pybind11::module)` instead of `pybind11_add_module(cadabra2 SHARED ...)`, and changing `core/packages/CMakeLists.txt` to use `pybind11_add_module(... MODULE ...)`.

Once installed, to use SusyPy, you simply need to write in any Python file

```python
import susypy
```

The optional text-message notifier `susypy.tools.send_msg` needs Twilio, which you can install with `pip install "susypy[notify]"`.


# Basic Use:

This section follows §5.1 of [the paper](https://arxiv.org/abs/2212.00614v1). The outputs shown are from SusyPy 1.1.0 with Cadabra2 2.6.0. (Cadabra2 2.6 prints index brackets as `(Γ^{a})_{α β}`, while older versions printed `\indexbracket(Γ^{a})_{α β}`.)

## **A Quick Tour of Cadabra2 in Python:**

SusyPy is a Python module built on Cadabra2. Cadabra2 has its own notebook language, but all of its functions can also be used from plain Python, which is what we use here. Import both packages with

```python
>>> import cadabra2 as cdb
>>> import susypy as susy
```

Input to Cadabra2 is a LaTeX-like string inside an `Ex()` wrapper. Sums are written with `+`, products by juxtaposition, Greek letters by their LaTeX names (e.g. `\Psi`), covariant indices with `_`, and contravariant indices with `^`. `Ex` objects can themselves be added and multiplied, and Cadabra2 prints them in Unicode:

```python
>>> ex1 = cdb.Ex(r'(A_{a} B_{b} + E_{a b}) \Psi^{c}')
>>> ex2 = cdb.Ex(r'G_{a} \Theta_{b}^{c}')
>>> ex = ex1 + ex2
>>> ex
(A_{a} B_{b} + E_{a b}) Ψ^{c} + G_{a} ϴ_{b}^{c}
```

Substitution rules are written as a comma-separated list of `->` rules in an `Ex()` wrapper, with the extra argument `False` to mark it as a set of rules rather than an ordinary expression. `substitute()` applies them. The indices in a rule adapt to the expression: new dummy indices are introduced as needed, and a rule may use index names that do not appear in the expression at all, since only the relative index structure matters.

```python
>>> sub = cdb.Ex(r'A_{a} -> G_{a}^{d} H_{d}, \Theta_{e}^{f} -> \Psi_{e}^{f}', False)
>>> cdb.substitute(ex, sub)
(G_{a}^{d} H_{d} B_{b} + E_{a b}) Ψ^{c} + G_{a} Ψ_{b}^{c}
```

Like Cadabra2's functions, almost all SusyPy functions modify the expression you pass in (as well as returning it), so you do not need to store the result in a new variable.

## **Setting up the Environment:**

Before doing any arithmetic, create an environment with `susy_env()`. It defines the tensor-spinors SusyPy needs and fixes the properties of the Clifford algebra representation you are working in. Its parameters are:

* `D`: the dimension.
* `lorentz_indices`: the Lorentz indices you intend to use (at least `D` of them, so that Levi-Civita tensors can be built).
* `spinor_indices`: the spinor indices you intend to use.
* `desired_syms`: a length-2 list (t₀, t₁) fixing the spinor-index symmetries of the spinor metric and the 1-index gamma matrix. It only matters in even dimensions, since the symmetries are fixed in odd dimensions.
* `rep`: the sign (±1) of the gamma-matrix reduction formula, which picks the Clifford algebra representation. It only matters in odd dimensions.

> :warning: The result of `susy_env()` must be stored in a variable named `__cdbkernel__`, which is how Cadabra2 finds the active kernel.

Calling `susy_env()` with no arguments gives 11D with the conventions of Appendix A of the paper (the environment used for the 11D supergravity calculations). The line below instead creates a 4D environment with Lorentz indices a, b, c, d and spinor indices α, β, γ, η, in which (t₀, t₁) = (1, 1) makes the spinor metric and the 1-index gamma matrix antisymmetric:

```python
>>> __cdbkernel__ = susy.susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma', r'\eta'], desired_syms=[1,1], rep=-1)
```

The environment defines the following objects:

| SusyPy Name | Cadabra Class | Use |
| --- | --- | --- |
| `\Gamma` | `GammaMatrix` | Dirac gamma matrix |
| `\delta` | `KroneckerDelta` | Kronecker delta |
| `\epsilon` | `EpsilonTensor` | Levi-Civita symbol |
| `D` | `Derivative` | Supercovariant derivative |
| `\partial` | `PartialDerivative` | Partial derivative |
| `I` | `ImaginaryI` | Imaginary unit |
| `Tr` | `Trace` | Trace function |
| `C` | N/A | Spinor metric |
| `\Gamma'` | N/A | Highest-rank element γ<sub>∗</sub> of the Clifford algebra |

To check that the environment has the properties we asked for, we use SusyPy's simplification function `evaluate()`; Cadabra2's `canonicalise()`, which uses symmetries to put indices in a standard lexicographic order; and Cadabra2's `rename_dummies()`, which renames dummy indices to the earliest available names. We expect the following:

* The trace (γ<sub>a</sub>)<sub>α</sub><sup>β</sup>(γ<sup>b</sup>)<sub>β</sub><sup>α</sup> = 2<sup>⌊D/2⌋</sup>δ<sub>a</sub><sup>b</sup> should give 4δ<sub>a</sub><sup>b</sup> in 4D.
* Swapping the spinor indices of the (antisymmetric) 1-index gamma matrix should give a minus sign.
* The dummy pairs d, γ should be renamed to a, α, the first indices in each list.

```python
>>> ex1 = cdb.Ex(r'(\Gamma_{a})_{\alpha}^{\beta} (\Gamma^{b})_{\beta}^{\alpha}')
>>> susy.evaluate(ex1)
4δ_{a}^{b}
>>> ex2 = cdb.Ex(r'(\Gamma^{a})_{\beta \alpha}')
>>> cdb.canonicalise(ex2)
-(Γ^{a})_{α β}
>>> ex3 = cdb.Ex(r'(A_{d} B^{d})_{\gamma}^{\gamma}')
>>> cdb.rename_dummies(ex3)
(A_{a} B^{a})_{α}^{α}
```

To change the environment, call `susy_env()` again with new parameters. Each call starts from a fresh Cadabra2 kernel, so properties from the previous environment (including any you declared yourself) are discarded; assign the result to `__cdbkernel__` again. Below we switch to D = 11 with Greek Lorentz indices, Latin spinor indices, and the "positive" representation of the 11D Clifford algebra. `desired_syms` is not needed in odd dimensions; here it is (1, −1), so the spinor metric is antisymmetric and the 1-index gamma matrix is symmetric.

```python
>>> __cdbkernel__ = susy.susy_env(D = 11, lorentz_indices=[r'\alpha', r'\beta', r'\gamma', r'\zeta', r'\eta', r'\theta', r'\iota', r'\kappa', r'\lambda', r'\mu', r'\nu'], spinor_indices=['a', 'b', 'c', 'd'], rep=1)
```

Now the trace should be 32δ<sub>α</sub><sup>β</sup>, the 1-index gamma matrix should no longer change sign when its indices are swapped, and the dummy indices should follow the new index lists. The positive representation also means a six-index gamma matrix reduces as

γ<sup>αβηγθζ</sup> = −(1/120) ϵ<sup>αβηγθζικλμν</sup> γ<sub>ικλμν</sub>.

```python
>>> ex1 = cdb.Ex(r'(\Gamma_{\alpha})_{a}^{b} (\Gamma^{\beta})_{b}^{a}')
>>> susy.evaluate(ex1)
32δ_{α}^{β}
>>> ex2 = cdb.Ex(r'(\Gamma^{\alpha})_{b a}')
>>> cdb.canonicalise(ex2)
(Γ^{α})_{a b}
>>> ex3 = cdb.Ex(r'(A_{\gamma} B^{\gamma})_{d}^{d}')
>>> cdb.rename_dummies(ex3)
(A_{α} B^{α})_{a}^{a}
>>> ex4 = cdb.Ex(r'\Gamma^{\alpha \beta \eta \gamma \theta \zeta}')
>>> susy.evaluate(ex4)
 -  1/120 ε^{α β η γ θ ζ ι κ λ μ ν} Γ_{ι κ λ μ ν}
```

## **Simplifying Tensor Expressions:**

At its core, SusyPy is a tool for explicit spinor-index arithmetic and canonicalization following the NW-SE convention. Spinor indices are distinguished from vector indices using Cadabra2's *index brackets*: a set of parentheses carrying indices, around the object they act on. So a tensor-spinor is written as its Lorentz tensor in parentheses, with the spinor indices on the parentheses. For example, Ψ<sub>a η</sub>, with Lorentz index a and spinor index η, is

```python
>>> __cdbkernel__ = susy.susy_env()
>>> cdb.Ex(r'''(\Psi_{a})_{\eta}''')
(Ψ_{a})_{η}
```

> NOTE: The mixed-index spinor metric is often written as δ<sub>α</sub><sup>β</sup>, but in SusyPy Kronecker deltas are pure Lorentz tensors, so you must use `C` for the spinor metric. `C` is the only tensor whose indices are read as spinor indices without an index bracket.

The function that does the heavy lifting with spinor indices is `spinor_combine()`. It uses symmetries, the spinor metric, and a careful ordering of tensor-spinor chains to combine index brackets as much as possible, following the NW-SE convention. Below, it contracts out the spinor metric and combines the two gamma matrices as an ordinary matrix product:

```python
>>> ex = cdb.Ex(r'(\Gamma_{a})_{\alpha \beta} (\Gamma^{b})^{\beta \gamma} C_{\gamma \eta}')
>>> susy.spinor_combine(ex)
-(Γ_{a} Γ^{b})_{α η}
```

The general simplification function is `evaluate()`. It:

1. applies `spinor_combine()` to combine index brackets as much as possible and put chains of index brackets in lexicographic order (canonicalizing the spinor-index structure as far as symmetries allow);
2. multiplies all gamma matrices using Cadabra2's Lorentz-tensor simplification;
3. simplifies the traces implicit in products of tensor-spinors, all derivatives, and all Kronecker-delta contractions;
4. pulls constants and pure Lorentz tensors (e.g. Kronecker deltas and Levi-Civita tensors) out of derivatives, Fourier transforms, and index brackets;
5. canonicalizes the Lorentz indices inside index brackets; and
6. flattens the expression into a sum of products of single-term factors.

The result is a standard, simplified form, so different expressions can be reliably compared. For example, finishing the product of gamma matrices above:

```python
>>> ex = cdb.Ex(r'(\Gamma_{a})_{\alpha \beta} (\Gamma^{b})^{\beta \gamma} C_{\gamma \eta}')
>>> susy.evaluate(ex)
-(Γ_{a}^{b})_{α η}-C_{α η} δ_{a}^{b}
```

The next example pulls coefficients and constants out of index brackets, contracts a spinor metric and a Kronecker delta, combines index brackets per the NW-SE convention, multiplies the gamma matrices, and canonicalizes the result:

```python
>>> ex = cdb.Ex(r'(\Gamma_{a})_{\alpha \beta} (I \Gamma^{b})^{\beta \gamma} C_{\gamma \eta} (5 \delta_{b}{}^{a} \Psi_{d})_{\zeta}')
>>> susy.evaluate(ex)
-55I C_{α η} (Ψ_{d})_{ζ}
```

## **Compatibility with Cadabra2:**

Almost all SusyPy functions take Cadabra2 `Ex` objects as input and return `Ex` objects, so you can freely pass expressions back and forth between Cadabra2 and SusyPy functions. SusyPy also works with Cadabra2's usual ways of defining new tensors and their properties.

In the 4D code below, three tensors A<sub>a</sub>, λ<sub>γ</sub>, and d (a pseudo-scalar) are declared with Cadabra2's `Depends` to have nonzero supercovariant and partial derivatives, and λ<sub>γ</sub> is declared to anticommute with the supercovariant derivative, since both are fermionic. Properties of a tensor with spinor indices must be given both to the underlying Lorentz tensor and to the tensor inside an index bracket, which Cadabra2 treats as a separate object; in property declarations you must write `\indexbracket` explicitly instead of using parentheses.

```python
>>> from cadabra2 import Ex
>>> __cdbkernel__ = susy.susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])
>>> cdb.Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
>>> cdb.Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
>>> cdb.Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
>>> cdb.AntiCommuting(Ex(r'''\lambda, D{#}'''))
>>> cdb.AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
>>> cdb.Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))
```

Now we combine SusyPy and Cadabra2 functions on a single expression: the anticommutator of supercovariant derivatives acting on λ<sub>γ</sub>. We apply SusyPy's `susy_expand()`, `evaluate()`, and `fierz_expand_2index()`, followed by Cadabra2's `factor_in()` to collect terms with respect to the constants u, v, w, x, y. (This computes the closure on the fermion of the off-shell 4D vector multiplet; see §5.9 of the paper.)

```python
>>> susy_rule = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + w C_{\alpha \beta} \partial^{a}(A_{a}) + x (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> y (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta})'''
>>> basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
>>> indices = [r'_{\alpha}', r'_{\beta}']
>>> ex = cdb.Ex(r'D_{\alpha}(D_{\beta}((\lambda)_{\gamma})) + D_{\beta}(D_{\alpha}((\lambda)_{\gamma}))')
>>> susy.susy_expand(ex, susy_rule)
>>> susy.evaluate(ex, to_perform_subs=False)
>>> susy.fierz_expand_2index(ex, basis, indices)
>>> cdb.factor_in(ex, cdb.Ex('u, v, w, x, y', False))
( -  3/2 u v -  1/2 u w +  1/2 x y) (Γ^{a})_{α β} \partial_{a}((λ)_{γ}) + ( 1/2 u v -  1/2 u w +  1/2 x y) (Γ^{a})_{α β} (Γ_{a}^{b})_{γ}^{\alpha'} \partial_{b}((λ)_{\alpha'}) + ( 1/2 u v +  1/2 u w +  1/2 x y) (Γ_{a}^{b})_{α β} (Γ^{a})_{γ}^{\alpha'} \partial_{b}((λ)_{\alpha'}) + ( 1/4 u v +  1/4 u w +  1/4 x y) I ε_{a b c}^{d} (Γ^{a b})_{α β} (Γ_{d} \Gamma')_{γ}^{\alpha'} \partial^{c}((λ)_{\alpha'})
```

(In the last term, the rank-3 gamma matrix (Γ<sub>ab</sub><sup>c</sup>) of the paper's output has been rewritten with the 4D duality relation in terms of ϵ and Γ<sub>d</sub>Γ′.)

## **Fierz Expansions:**

To get the Fierz expansion of an expression, pass the expression, the basis, and the spinor indices of the factors in the expansion to `fierz_expand()`. For example, the following code generates (A.24) of [Gates2019](https://doi.org/10.1007/JHEP07(2019)063) in the default 11D environment:

```python
>>> from cadabra2 import Ex
>>> from susypy import fierz_expand, susy_env
>>> __cdbkernel__ = susy_env()
>>> ex = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} + \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
>>> basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
>>> commuted = [r'_{\alpha}', r'_{\beta}']
>>> uncommuted = [r'_{\eta}', r'_{\gamma}']
>>> fierz_expand(ex, basis, commuted, uncommuted)
 -  1/16 (Γ^{a})_{α β} (Γ_{a})^{η γ} +  1/32 (Γ^{a b})_{α β} (Γ_{a b})^{η γ} -  1/1920 (Γ^{a b c d e})_{α β} (Γ_{a b c d e})^{η γ}
```

## **More Examples:**

The [`paper/`](https://github.com/IsaiahBHilz/susypy/tree/main/paper) folder has the full code for every calculation in the paper (the files with `code_blocks` in their names), including the multiplet solvers, the SUSY-invariant action solver, propagators, and the holoraumy of 11D, N = 1 supergravity.


# License:

SusyPy is released under the GNU General Public License, version 3 or (at your option) any later version; see [`COPYING.txt`](https://github.com/IsaiahBHilz/susypy/blob/main/COPYING.txt).
