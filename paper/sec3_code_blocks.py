import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric, factor_in
from susypy import evaluate, susy_env, susy_expand, make_action_susy_inv

__cdbkernel__ = susy_env()

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''h{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''c{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''\indexbracket{\zeta{#}}{#}'''), Ex(r'''\partial{#}'''))
Symmetric(Ex(r'''h{#}'''))
AntiSymmetric(Ex(r'''A{#}'''))
AntiCommuting(Ex(r'''\Psi, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\Psi{#}}{#}, D{#}'''))

susy = r'D_{\alpha}(h_{a b}) -> u ((\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}), D_{\alpha}((\Psi_{b})^{\beta}) -> 2 v \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + x (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) + y (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> 2 z (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - 2 z (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + 2 z (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}, c_{a b c} -> \partial_{a}(h_{b c}) - \partial_{b}(h_{a c})'

# SECTION 3
## Section 3.1
### Sec. 3.1: Code Block #1
L = Ex(r'D_{\gamma}(l (-1/4 c^{a b c} c_{a b c} + 1/2 c^{a b c} c_{c a b} + c^{a}_{b}^{b} c_{a c}^{c}) + m (1/12) (\Psi_{a})^{\alpha} (\Gamma^{a b c})_{\alpha \beta} \partial_{b}((\Psi_c)^{\beta}) + n (3/4) (\partial_{a}(A_{b c d}) - \partial_{b}(A_{a c d}) + \partial_{c}(A_{a b d}) - \partial_{d}(A_{a b c})) (\partial^{a}(A^{b c d}) - \partial^{b}(A^{a c d}) + \partial^{c}(A^{a b d}) - \partial^{d}(A^{a b c})))')

### Sec. 3.1: Code Block #2
sol = make_action_susy_inv(L, susy, ['u', 'v', 'x', 'y', 'z', 'l', 'm', 'n'], [r'\Psi'])
print('\n', sol)

## Section 3.2
### Sec. 3.2: Code Block #1
susy_expand(L, susy)
evaluate(L)
factor_in(L, Ex('u, v, x, y, z, l, m, n', False))
print('\n', L)