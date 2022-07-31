import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric, factor_in
from susypy import evaluate, susy_env, susy_expand

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

# SECTION 6
## Section 6.1
### Sec. 6.1: Code Block #1
susy = r'D_{\alpha}(h_{a b}) -> (\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}, D_{\alpha}((\Psi_{b})^{\beta}) -> -1/4 I \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + 1/48 I (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) - 1/144 I (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}'

### Sec. 6.1: Code Block #2
ex = Ex(r'D_{\alpha}(D_{\beta}(h_{a b})) - D_{\beta}(D_{\alpha}(h_{a b}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)

### Sec. 6.1: Code Block #3
ex = Ex(r'D_{\alpha}(D_{\beta}(A_{a b c})) - D_{\beta}(D_{\alpha}(A_{a b c}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)

## Section 6.2
### Sec. 6.2: Code Block #1
ex = Ex(r'D_{\alpha}(D_{\beta}(\indexbracket(\Psi_{a})^{\gamma})) - D_{\beta}(D_{\alpha}(\indexbracket(\Psi_{a})^{\gamma}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)