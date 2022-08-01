<h1 align="center">
<img src="/branding/logo/logo.png" width="150"> 
<br> SusyPy
</h1><br>

**SusyPy** (short for **Su**per**sy**mmetry **Py**thon) is a symbolic algebra system, built on top of Cadabra2, for performing intricate tensor arithmetic in supersymmetry.

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

To install, you will need the SusyPy wheel file, which has the name `susypy-1.0.0-py3-none-any.whl`. Now, if you are using pip, go to your terminal and run the command

```bash
pip install susypy-1.0.0-py3-none-any.whl
```

Once installed, to use SusyPy, you simply need to write in any Python file

```python
import susypy
```

> :warning: WARNING: SusyPy will only run on Python 3.8 or newer, so make sure to update your version of Python.

> TIP: To make sure that you are installing SusyPy into the right version of Python, use `pip3.8` instead of `pip` in the installation command.

# User Guide:
## **Getting Started:**

Input into SusyPy is given as LaTeX strings into `Ex()` wrappers. All expressions, whether vectorial or spinorial, are evaluated with the `evaluate()` function. Below is a simple example of evaluating a product of gamma matrices and the spinor metric. Notice that the algorithm enforces NW-SE convention.

```python
>>> ex = Ex(r'(\Gamma_{a})_{\alpha \beta} (\Gamma^{b})^{\beta \gamma} C_{\gamma \eta}')
>>> evaluate(ex)
'-\indexbracket(Γ_{a}^{b})_{α η}-C_{α η} δ_{a}^{b}'
```

In order to generate the Fierz expansion of an expression, one can input the expression, the basis, and the desired spinor indices of the factors in the expansion into `fierz_expand()`. For example, the following code generates (A.24) in [Gates2019](https://doi.org/10.1007/JHEP07(2019)063).

```python
>>> ex = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} + \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
>>> basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
>>> commuted = [r'_{\alpha}', r'_{\beta}']
>>> uncommuted = [r'_{\eta}', r'_{\gamma}']
>>> fierz_expand(ex, basis, commuted, uncommuted)
```

(To do: use this subsection to demonstrate the other algorithms.)