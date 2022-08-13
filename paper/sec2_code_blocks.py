import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric
from cadabra2 import canonicalise, collect_terms, eliminate_kronecker, split_gamma
from susypy import evaluate, evaluate_traces, fierz_expand, fierz_expand_2index, find_non_gauge_inv, fourier, inverse_fourier, make_action_susy_inv, rarita_schwinger_prop, susy_env, susy_expand, susy_solve, sym_proj_decomp, sym_prop, susy_solve_propagator
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

## Section 2.8
### Sec. 2.8: Code Block #1
__cdbkernel__ = susy_env(D=4)
AntiSymmetric(Ex(r'''B{#}'''))
Depends(Ex(r'''B{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''D{#}, \partial{#}'''))
field = Ex(r'B_{a b}')
gauge_trans = Ex(r'\partial_{a}(\zeta_{b}) - \partial_{b}(\zeta_{a})')
ex = Ex(r'(u + v) (\Gamma^{a})_{\alpha \beta} \partial_{a}(B_{b c}) - v (\Gamma^{a})_{\alpha \beta} \partial_{b}(B_{a c}) + (u + v) (\Gamma^{a})_{\alpha \beta} \partial_{c}(B_{a b})')
print('\n', find_non_gauge_inv(ex, field, [field], gauge_trans))

## Section 2.9
### Sec. 2.9: Code Block #1
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))

bosons = [Ex('A'), Ex('B')]
fermions = [Ex(r'(\Psi)_{\gamma}')]
gauge_transs = [Ex(r'0')]
susy = r'''D_{\alpha}(A) -> u (\Psi)_{\alpha}, D_{\alpha}(B) -> v (\Gamma')_{\alpha}^{\beta} (\Psi)_{\beta}, D_{\alpha}((\Psi)_{\beta}) -> w (\Gamma^{a})_{\alpha \beta} \partial_{a}(A) + x (\Gamma' \Gamma^{a})_{\alpha \beta} \partial_{a}(B)'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v', 'w', 'x']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = susy_solve(bosons, fermions, gauge_transs, susy, basis, consts, indices, comm_coef=2)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

### Sec. 2.9: Code Block #2
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])
Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))

bosons = [Ex('A_{a}')]
fermions = [Ex(r'(\lambda)_{\gamma}')]
gauge_transs = [Ex(r'0')]
susy = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = susy_solve(bosons, fermions, gauge_transs, susy, basis, consts, indices, comm_coef=2)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

### Sec. 2.9: Code Block #3
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\lambda, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

bosons = [Ex('A_{a}'), Ex(r'd')]
fermions = [Ex(r'(\lambda)_{\gamma}')]
gauge_transs = [Ex(r'0')]
susy = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + w C_{\alpha \beta} \partial^{a}(A_{a}) + x (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> y (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v', 'w', 'x', 'y']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = susy_solve(bosons, fermions, gauge_transs, susy, basis, consts, indices, comm_coef=2)
print('\n', susy_sol)
print('\n')
pretty_print(susy_data)

## Section 2.10
__cdbkernel__ = susy_env(D=4)

### Sec. 2.10: Code Block #1
subs = sym_proj_decomp(3)
print('\n', subs)

### Sec. 2.10: Code Block #2
subs = sym_proj_decomp(3/2)
print('\n', subs)

### Sec. 2.10: Code Block #3
prop = sym_prop(2)
print('\n', prop)

### Sec. 2.10: Code Block #4
prop = sym_prop(3/2)
print('\n', prop)

### Sec. 2.10: Code Block #5
ex = sym_prop(3/2)
split_gamma(ex, on_back=True)
split_gamma(ex, on_back=True)
evaluate(ex, to_join_gamma=False)
print('\n', ex)

### Sec. 2.10: Code Block #6
__cdbkernel__ = susy_env(D=11)
prop = rarita_schwinger_prop()
print('\n', prop)

## Section 2.11
### Sec. 2.11: Equation 2.11.5
__cdbkernel__ = susy_env(D=4)
print(evaluate(rarita_schwinger_prop()))

### Sec. 2.11: Code Block #1 (Axial-Vector)
Depends(Ex(r'''U{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))

bosons = [Ex('U_{a}')]
fermions = [Ex(r'(\lambda)_{\gamma}')]
fermion_propagators = [sym_prop(1/2)]
susy = r'''D_{\alpha}(U_{a}) -> u (\Gamma' \Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma' \Gamma^{a b})_{\alpha \beta} \partial_{a}(U_{b})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

### Sec. 2.11: Code Block #2 (Matter-Gravitino)
Depends(Ex(r'''B{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''\indexbracket{\zeta{#}}{#}'''), Ex(r'''\partial{#}'''))

bosons = [Ex('B_{a}')]
fermions = [Ex(r'(\Psi_{a})_{\gamma}')]
fermion_propagators = [rarita_schwinger_prop()]
susy = r'''D_{\alpha}(B_{a}) -> u (\Psi_{a})_{\alpha}, D_{\alpha}((\Psi_{a})_{\beta}) -> v (\Gamma_{a} \Gamma^{b c})_{\alpha \beta} \partial_{b}(B_{c})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

### Sec. 2.11: Code Block #3 (Supergravity)
Depends(Ex(r'''h{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''\indexbracket{\zeta{#}}{#}'''), Ex(r'''\partial{#}'''))
Symmetric(Ex(r'''h{#}'''))

bosons = [Ex('h_{a b}')]
fermions = [Ex(r'(\Psi_{a})_{\gamma}')]
fermion_propagators = [rarita_schwinger_prop()]
susy = r'''D_{\alpha}(h_{a b}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\Psi_{b})_{\beta} + u (\Gamma_{b})_{\alpha}^{\beta} (\Psi_{a})_{\beta}, D_{\alpha}((\Psi_{a})_{\beta}) -> v (\Gamma^{b c})_{\alpha \beta} \partial_{b}(h_{c a})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

## Section 2.12
### Sec. 2.12: Code Block #1
__cdbkernel__ = susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])

Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\lambda, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

L = Ex(r'D_{\gamma}(l F_{a b} F^{a b} + m (\Gamma^{a})^{\alpha \beta} (\lambda)_{\alpha} \partial_{a}((\lambda)_{\beta}) + n d d)')
susy = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + w C_{\alpha \beta} \partial^{a}(A_{a}) + x (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> y (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta}), F_{a b} -> \partial_{a}(A_{b}) - \partial_{b}(A_{a})'''
sol = make_action_susy_inv(L, susy, ['l', 'm', 'n', 'u', 'v', 'w', 'x', 'y'], [r'\lambda'])
print('\n', sol)