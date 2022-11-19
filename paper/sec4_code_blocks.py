import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric, factor_in
from susypy import evaluate, fierz_expand_2index, holoraumy, susy_env, susy_expand
from susypy.tools import pretty_print

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

# SECTION 4
## Section 4.0
### Sec. 4.0: Code Block #1
fields = [Ex('h_{a b}'), Ex('A_{a b c}'), Ex(r'(\Psi_{a})^{\gamma}')]
susy = r'D_{\alpha}(h_{a b}) -> u ((\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}), D_{\alpha}((\Psi_{b})^{\beta}) -> 2 v \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + x (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) + y (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> 2 z (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - 2 z (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + 2 z (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}, c_{a b c} -> \partial_{a}(h_{b c}) - \partial_{b}(h_{a c})'
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
subs = Ex('u -> 1, v -> (-1/8) I, x -> (1/48) I, y -> (-1/144) I, z -> 1/2, l -> 1/4, m -> 24 I, n -> -1/36', False)
indices = [r'_{\alpha}', r'_{\beta}']
hol_data = holoraumy(fields, susy, basis, subs, indices)
print('\n')
pretty_print(hol_data)

## Section 4.1
### Sec. 4.1: Eq. 6.1.1
susy = r'D_{\alpha}(h_{a b}) -> (\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}, D_{\alpha}((\Psi_{b})^{\beta}) -> -1/4 I \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + 1/48 I (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) - 1/144 I (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}'
ex = Ex(r'D_{\alpha}(D_{\beta}(h_{a b})) - D_{\beta}(D_{\alpha}(h_{a b}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)

### Sec. 4.1: Eq. 6.1.3
susy = r'D_{\alpha}(h_{a b}) -> (\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}, D_{\alpha}((\Psi_{b})^{\beta}) -> -1/4 I \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + 1/48 I (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) - 1/144 I (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}'
ex = Ex(r'D_{\alpha}(D_{\beta}(A_{a b c})) - D_{\beta}(D_{\alpha}(A_{a b c}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)

## Section 4.2
### Sec. 4.2: Code Block #1
susy = r'D_{\alpha}(h_{a b}) -> (\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}, D_{\alpha}((\Psi_{b})^{\beta}) -> -1/4 I \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + 1/48 I (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) - 1/144 I (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}'
ex = Ex(r'D_{\alpha}(D_{\beta}(\indexbracket(\Psi_{a})^{\gamma})) - D_{\beta}(D_{\alpha}(\indexbracket(\Psi_{a})^{\gamma}))')
susy_expand(ex, susy)
evaluate(ex)
print('\n', ex)

### Sec. 4.2: Code Block #2
fierz_expand_2index(ex, basis, indices)
print('\n', ex)