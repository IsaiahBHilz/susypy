**SusyPy** (short for **Su**per**sy**mmetry **Py**thon) offers a symbolic algebra system for supersymmetry calculations, particularly for 10-11D supergravity. Given the basis for a multiplet (by default that of the 11D supergravity multiplet), the algorithms of the library can complete intricate tensor arithmetic.

## **Features:**
* The program extends ordinary arithmetic in Python to gamma matrices and other `SUSYTensor`s.
* The algorithm simplifies expressions completely, giving results linear in the basis elements.
* Commutators can be treated in arithmetic without expansion, giving more succinct and useful solutions.

# Installation:

To install, you will need the SusyPy wheel file, which has the name `susypy-0.9.2-py3-none-any.whl`. Now, if you are using pip, go to your terminal and run the command

```bash
pip install susypy-0.9.2-py3-none-any.whl
```

Once installed, to use SusyPy, you simply need to write in any Python file

```python
import susypy
```

.. warning::
	SusyPy will only run on Python 3.9 or newer, so make sure to update your version of Python.

.. tip::
	To make sure that you are installing SusyPy into the right version of Python, use `pip3.9` instead of `pip` in the installation command.

# User Guide:
## **Defining Indices and Gammas:**
Here are a few simple examples to demonstrate SusyPy's tools. First, the Lorentz and Dirac index types must be imported and tensor indices created. The example below is in 11D, so the 11D Lorentz and Dirac index types must be imported from the basis module. (See more information on the basis [here](basis/index.html).)

```python
from susypy import tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b = tensor_indices('a, b', Lorentz11d)
alpha, beta = tensor_indices('alpha, beta', Dirac11d)
```

The `Gamma` class can be used to create a gamma matrix. For example, below is the code for \((\gamma^a)^{\alpha\beta}\).

```python
from susypy import Gamma

u = Gamma([a], [alpha, beta])
```

By default, tensor indices are contravariant. They can be negated to give covariant indices. For example, below is the code for \((\gamma_a)_\alpha{}^\beta\).

```python
v = Gamma([-a], [-alpha, beta])
```

## **Addition:**

Below are created two (identical) gamma matrices, each with two contravariant Lorentz indices, one covariant spinor index, and one contravariant spinor index.

```python
u = Gamma([a, b], [-alpha, beta])
v = Gamma([a, b], [-alpha, beta])
```

These gamma matrices can be added.

```python
>>> w = u+v
>>> w
2.00000000000000*gamma{a, b, -alpha, beta}
```

If the indices do not match, adding the gamma matrices will raise an error.

```python
>>> c = tensor_indices('c', Lorentz11d)
>>> u = Gamma([a, b], [-alpha, beta])
>>> v = Gamma([a, c], [-alpha, beta])
>>> w = u+v
Exception: Cannot add terms with different indices.
```

Full code:

```python
from susypy import Gamma, tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b = tensor_indices('a, b', Lorentz11d)
alpha, beta = tensor_indices('alpha, beta', Dirac11d)

u = Gamma([a, b], [-alpha, beta])
v = Gamma([a, b], [-alpha, beta])
w = u+v
```

## **Multiplication:**

Gamma matrices can also be multiplied. The result may involve commutators.

```python
>>> a, b, c = tensor_indices('a, b, c', Lorentz11d)
>>> alpha, beta, gamma = tensor_indices('alpha, beta, gamma', Dirac11d)

>>> u = Gamma([a], [-alpha, gamma])
>>> v = Gamma([b, c], [-gamma, beta])
>>> w = u*v
>>> w
gamma{a, b, c, -alpha, beta} + comm(eta{a, b, -alpha, zeta}*gamma{c, -zeta, beta}, [b, c])
```

In the result above, `comm(eta{a, b, -alpha, zeta}*gamma{c, -zeta, beta}, [b, c])` means \((\eta^{a[b})_\alpha{}^\zeta(\gamma^{c]})_\zeta{}^\beta = (\eta^{a[b}\gamma^{c]})_\alpha{}^\beta\). In other words, `eta{a, b, -alpha, zeta}*gamma{c, -zeta, beta}` is commuted with respect to the indices `b, c`. If northwest-southeast convention is not satisfied, multiplying the gamma matrices will raise an error.

```python
>>> u = Gamma([a], [-alpha, gamma])
>>> v = Gamma([b, c], [gamma, beta])
>>> w = u*v
Exception: Northwest-southeast convention not satisfied.
```

Full code:

```python
from susypy import Gamma, tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b, c = tensor_indices('a, b, c', Lorentz11d)
alpha, beta, gamma = tensor_indices('alpha, beta, gamma', Dirac11d)

u = Gamma([a], [-alpha, gamma])
v = Gamma([b, c], [-gamma, beta])
w = u*v
```

## **Other SUSYTensors:**
In addition to gamma matrices, calculations may involve metrics and Levi-Civita symbols, or even commutator inputs. For example, below is the code for \((\eta^{ab})_\alpha{}^\beta\).

```python
from susypy import Gamma, tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b = tensor_indices('a, b', Lorentz11d)
alpha, beta = tensor_indices('alpha, beta', Dirac11d)

eta = SUSYTensor([a, b], [-alpha, beta], 'eta', L=Lorentz11d.metric(a, b))
```

In the last line, the "Lorentz part" of the `SUSYTensor` is given as `L=Lorentz.metric(a,b)`, the Lorentz metric tensor with contravariant indices `a, b`. Notice that the `SUSYTensor` must be given a name, in this case, `'eta'`. Below is the code for \((\epsilon^{abcdefghijk})_\alpha{}^\beta\).

```python
from susypy import Gamma, tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b, c, d, e, f, g, h, i, j, k = tensor_indices('a, b, c, d, e, f, g, h, i, j, k', Lorentz11d)
alpha, beta = tensor_indices('alpha, beta', Dirac11d)

epsilon = SUSYTensor([a, b, c, d, e, f, g, h, i, j, k], [-alpha, beta], 'epsilon', L=Lorentz11d.epsilon(a, b))
```

Finally, commutators can be defined using the `Commutator` class. Below is the code for \((\eta^{a[b})_\alpha{}^\zeta(\gamma^{c]})_\zeta{}^\beta\)

```python
from susypy import Gamma, tensor_indices
from susypy.basis import Lorentz11d, Dirac11d

a, b, c = tensor_indices('a, b, c', Lorentz11d)
alpha, beta, zeta = tensor_indices('alpha, beta, zeta', Dirac11d)

eta = SUSYTensor([a, b], [-alpha, zeta], 'eta', L=Lorentz11d.metric(a, b))
gamma1 = Gamma([c], [-zeta, beta])
comm = Commutator([a, b, c], [-alpha, beta], [b, c], eta*gamma1)
```

In the last line, `a, b, c` are the Lorentz indices of the commuted expression, `-alpha, beta` are the spinor indices of the commuted expression, `b, c` are the Lorentz indices with respect to which the expression is commuted, and `eta*gamma1` is the commuted expression. (Compare with the result from the multiplication example.)

## **Setting New Defaults:**
By default, the gamma matrices are written in terms of the 11D basis, and internal index manipulations use the Lorentz11d and Dirac11d index types. To set new defaults, use `set_defaults()`. Below is the code to set the 11D basis as default.

```python
from susypy.basis import Lorentz11d, Dirac11d, basis11d

set_defaults(Lorentz11d, Dirac11d, basis11d)
```

The defaults can also be altered for a particular class or list of classes. For example, the code below sets the default for the `Gamma` class.

```python
from susypy import Gamma
from susypy.basis import Lorentz11d, Dirac11d, basis11d

set_defaults(Lorentz11d, Dirac11d, basis11d, object_types=[Gamma])
```

Now that some demonstrations have been given, the key classes and functions are
expounded below.