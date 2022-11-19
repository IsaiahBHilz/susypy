import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric, factor_in
from susypy import evaluate, fierz_expand_2index, susy_env, rarita_schwinger_prop, susy_expand, susy_solve, susy_solve_propagator
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

# SECTION 2
## Section 2.1
### Sec. 2.1: Code Block #1
susy = r'D_{\alpha}(h_{a b}) -> u ((\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}), D_{\alpha}((\Psi_{b})^{\beta}) -> 2 v \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + x (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) + y (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> 2 z (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - 2 z (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + 2 z (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}, c_{a b c} -> \partial_{a}(h_{b c}) - \partial_{b}(h_{a c})'

### Sec. 2.1: Code Block #2
bosons = [Ex('h_{a b}'), Ex('A_{a b c}')]
fermions = [Ex(r'(\Psi_{a})^{\gamma}')]
gauge_transs = [Ex(r'\partial_{a}((\zeta)^{\gamma})')]
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
consts = ['u', 'v', 'x', 'y', 'z']
indices = [r'_{\alpha}', r'_{\beta}']
susy_sol, susy_data = susy_solve(bosons, fermions, gauge_transs, susy, basis, consts, indices)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

### Sec. 2.1: Code Block #3
fermion_propagators = [rarita_schwinger_prop()]
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices)
print('\n', sol)

## Section 2.2
### Sec. 2.2: Code Block #1
ex = Ex(r'D_{\alpha}(D_{\beta}(h_{a b})) + D_{\beta}(D_{\alpha}(h_{a b}))')
susy_expand(ex, susy)
evaluate(ex)
factor_in(ex, Ex('u, v, x, y, z', False))
print('\n', ex)

## Section 2.3
### Sec. 2.3: Code Block #1
ex = Ex(r'D_{\alpha}(D_{\beta}(A_{a b c})) + D_{\beta}(D_{\alpha}(A_{a b c}))')
susy_expand(ex, susy)
evaluate(ex)
factor_in(ex, Ex('u, v, x, y, z', False))
print('\n', ex)

## Section 2.4
### Sec. 2.4: Code Block #1
ex = Ex(r'D_{\alpha}(D_{\beta}((\Psi_{a})^{\gamma})) + D_{\beta}(D_{\alpha}((\Psi_{a})^{\gamma}))')
susy_expand(ex, susy)
evaluate(ex)
factor_in(ex, Ex('u, v, x, y, z', False))
print('\n', ex)

### Sec. 2.4: Code Block #2
ex = Ex(r'D_{\alpha}(D_{\beta}((\Psi_{a})^{\gamma})) + D_{\beta}(D_{\alpha}((\Psi_{a})^{\gamma}))')
susy_expand(ex, susy)
evaluate(ex)
fierz_expand_2index(ex, [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'])
factor_in(ex, Ex('u, v, x, y, z', False))
print('\n', ex)