import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, Depends, Ex
from cadabra2 import canonicalise, eliminate_kronecker, collect_terms
from susypy import evaluate, evaluate_traces, fierz_expand, fierz_expand_2index, fourier, inverse_fourier, make_action_susy_inv, susy_env, susy_expand, susy_solve
from susypy.tools import pretty_print

# SECTION 2
## Section 2.2
__cdbkernel__ = susy_env()

### Sec. 2.2: Code Block #1
ex = Ex(r'(\Gamma_{a b} \Gamma^{c d})_{\alpha}^{\alpha}')
evaluate_traces(ex)
print('\n', ex)

### Sec. 2.2: Code Block #2
ex = Ex(r'(-1/32) (\Gamma_{a} \Gamma_{b} \Gamma_{c} \Gamma_{d} \Gamma_{e} \Gamma_{f} \Gamma_{g} \Gamma_{h} \Gamma_{i} \Gamma_{j} \Gamma_{k})^{\alpha}_{\alpha}')
evaluate_traces(ex)
eliminate_kronecker(ex)
canonicalise(ex)
collect_terms(ex)
print('\n', ex)


## Section 2.5
__cdbkernel__ = susy_env()

### Sec. 2.5: Code Block #1
ex = Ex(r'\delta_{\alpha}^{\eta} (\Gamma_{a})^{\rho \gamma} - \delta_{\alpha}^{\rho} (\Gamma_{a})^{\eta \gamma}')
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
fierz_expand(ex, basis, [r'_{\alpha}', r'_{\gamma}'], [r'_{\eta}', r'_{\rho}'])
print('\n', ex)

### Sec. 2.5: Code Block #2
ex = Ex(r'(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}')
fierz_expand(ex, basis, [r'_{\alpha}', r'_{\gamma}'], [r'_{\gamma}', r'_{\eta}'])
print('\n', ex)

## Section 2.6
### Sec. 2.6: Code Block #1
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\lambda, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

susy = r'''D_{\alpha}(A_{a}) -> (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> -I (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> I (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
ex = Ex(r'D_{\alpha}(D_{\beta}((\lambda)_{\gamma})) + D_{\beta}(D_{\alpha}((\lambda)_{\gamma}))')
susy_expand(ex, susy)
print('\n', ex)
fierz_expand_2index(ex, basis, [r'_{\alpha}', r'_{\beta}'])
print('\n', ex)

### Sec. 2.6: Code Block #2
__cdbkernel__ = susy_env()
ex = Ex(r'''(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}''')
fierz_expand_2index(ex, basis, [r'_{\alpha}', r'_{\beta}'])
print('\n', ex)

## Section 2.7
### Sec. 2.7: Code Block #1
__cdbkernel__ = susy_env()
Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''C'''), Ex(r'''D{#}, \partial{#}'''))

ex = Ex(r'''\partial_{a}(A B) \partial_{b c}(C) + \partial_{a b c}(A B C)''')
fourier(ex, [Ex('A'), Ex('B'), Ex('C')])
canonicalise(ex)
print('\n', ex)

### Sec. 2.7: Code Block #2
inverse_fourier(ex)
canonicalise(ex)
print('\n', ex)

## Section 2.9
### Sec. 2.9: Code Block #1
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))

fields = fields = [Ex('A'), Ex('B'), Ex(r'(\Psi)_{\gamma}')]
gauge_transs = [Ex(r'0'), Ex(r'0'), Ex(r'0')]
susy = r'''D_{\alpha}(A) -> u (\Psi)_{\alpha}, D_{\alpha}(B) -> v (\Gamma')_{\alpha}^{\beta} (\Psi)_{\beta}, D_{\alpha}((\Psi)_{\beta}) -> w (\Gamma^{a})_{\alpha \beta} \partial_{a}(A) + x (\Gamma' \Gamma^{a})_{\alpha \beta} \partial_{a}(B)'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v', 'w', 'x']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = susy_solve(fields, gauge_transs, susy, basis, consts, indices, comm_coef=2)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

### Sec. 2.9: Code Block #2
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\lambda, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

fields = [Ex('A_{a}'), Ex(r'd'), Ex(r'(\lambda)_{\gamma}')]
gauge_transs = [Ex(r'\partial_{a}(\zeta)'), Ex(r'0'), Ex(r'0')]
susy = r'''D_{\alpha}(A_{a}) -> p (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> q (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + t C_{\alpha \beta} \partial^{a}(A_{a}) + r (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> s (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['p', 'q', 'r', 's', 't']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = susy_solve(fields, gauge_transs, susy, basis, consts, indices, comm_coef=2)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

## Section 2.12
### Sec. 2.12: Code Block #1
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\lambda, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

L = Ex(r'D_{\gamma}(m F_{a b} F^{a b} + n (\Gamma^{a})^{\alpha \beta} (\lambda)_{\alpha} \partial_{a}((\lambda)_{\beta}) + o d d)')
susy = r'''D_{\alpha}(A_{a}) -> p (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> q (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + t C_{\alpha \beta} \partial^{a}(A_{a}) + r (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> s (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta}), F_{a b} -> \partial_{a}(A_{b}) - \partial_{b}(A_{a})'''
sol = make_action_susy_inv(L, susy, ['m', 'n', 'o', 'p', 'q', 'r', 's', 't'], [r'\lambda'])
print('\n', sol)