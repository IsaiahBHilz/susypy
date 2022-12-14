import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric
from cadabra2 import canonicalise, collect_terms, decompose_product, distribute, eliminate_kronecker, factor_in, split_gamma
from susypy import evaluate, evaluate_traces, fierz_expand, fierz_expand_2index, find_non_gauge_inv, fourier, holoraumy, indexbracket_hex, inverse_fourier, make_action_susy_inv, rarita_schwinger_prop, substitute, susy_env, susy_expand, susy_solve, sym_proj_decomp, sym_prop, susy_solve_propagator
from susypy.tools import pretty_print

# SECTION 5
## Section 5.1
### Sec. 5.1: Code Block #1
import cadabra2 as cdb
import susypy as susy

### Sec. 5.1: Code Block #2
ex1 = cdb.Ex(r'(A_{a} B_{b} + E_{a b}) \Psi^{c}') 
ex2 = cdb.Ex(r'G_{a} \Theta_{b}^{c}')
ex = ex1 + ex2
print('\n', ex)

### Sec. 5.1: Code Block #3
sub = cdb.Ex(r'A_{a} -> G_{a}^{d} H_{d}, \Theta_{e}^{f} -> \Psi_{e}^{f}', False)
cdb.substitute(ex, sub)
print('\n', ex)

### Sec. 5.1: Code Block #4
__cdbkernel__ = susy.susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma', r'\eta'], desired_syms=[1,1], rep=-1)

### Sec. 5.1: Code Block #5
ex1 = cdb.Ex(r'(\Gamma_{a})_{\alpha}^{\beta} (\Gamma^{b})_{\beta}^{\alpha}')
susy.evaluate(ex1)
print('\n', ex1)
ex2 = cdb.Ex(r'(\Gamma^{a})_{\beta \alpha}')
cdb.canonicalise(ex2)
print('\n', ex2)
ex3 = cdb.Ex(r'(A_{d} B^{d})_{\gamma}^{\gamma}')
cdb.rename_dummies(ex3)
print('\n', ex3)

### Sec. 5.1: Code Block #6
__cdbkernel__ = susy.susy_env(D = 11, lorentz_indices=[r'\alpha', r'\beta', r'\gamma', r'\zeta', r'\eta', r'\theta', r'\iota', r'\kappa', r'\lambda', r'\mu', r'\nu'], spinor_indices=['a', 'b', 'c', 'd'], rep=1)

### Sec. 5.1: Code Block #7
ex1 = cdb.Ex(r'(\Gamma_{\alpha})_{a}^{b} (\Gamma^{\beta})_{b}^{a}')
susy.evaluate(ex1)
print('\n', ex1)
ex2 = cdb.Ex(r'(\Gamma^{\alpha})_{b a}')
cdb.canonicalise(ex2)
print('\n', ex2)
ex3 = cdb.Ex(r'(A_{\gamma} B^{\gamma})_{d}^{d}')
cdb.rename_dummies(ex3)
print('\n', ex3)
ex4 = cdb.Ex(r'\Gamma^{\alpha \beta \eta \gamma \theta \zeta}')
susy.evaluate(ex4)
print('\n', ex4)

### Sec. 5.1: Code Block #8
__cdbkernel__ = susy.susy_env()
cdb.Ex(r'''(\Psi_{a})_{\eta}''')


### Sec. 5.1: Code Block #9
ex = cdb.Ex(r'(\Gamma_{a})_{\alpha \beta} (\Gamma^{b})^{\beta \gamma} C_{\gamma \eta}')
susy.spinor_combine(ex)
print('\n', ex)

### Sec. 5.1: Code Block #10
ex = Ex(r'(\Gamma_{a})_{\alpha \beta} (I \Gamma^{b})^{\beta \gamma} C_{\gamma \eta} (5 \delta_{b}{}^{a} \Psi_{d})_{\zeta}')
susy.evaluate(ex)
print('\n', ex)

### Sec. 5.1: Code Block #11
__cdbkernel__ = susy.susy_env(D = 4, lorentz_indices=['a', 'b', 'c', 'd'], spinor_indices=[r'\alpha', r'\beta', r'\gamma'])
cdb.Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
cdb.Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
cdb.Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
cdb.AntiCommuting(Ex(r'''\lambda, D{#}'''))
cdb.AntiCommuting(Ex(r'''\indexbracket{\lambda{#}}{#}, D{#}'''))
cdb.Depends(Ex(r'''d'''), Ex(r'''D{#}, \partial{#}'''))

susy_rule = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b}) + w C_{\alpha \beta} \partial^{a}(A_{a}) + x (\Gamma')_{\alpha \beta} d, D_{\alpha}(d) -> y (\Gamma' \Gamma^{a})_{\alpha}^{\beta} \partial_{a}((\lambda)_{\beta})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v', 'w', 'x', 'y']
indices = [r'_{\alpha}', r'_{\beta}']

ex = cdb.Ex(r'D_{\alpha}(D_{\beta}((\lambda)_{\gamma})) + D_{\beta}(D_{\alpha}((\lambda)_{\gamma}))')
susy.susy_expand(ex, susy_rule)
susy.evaluate(ex, to_perform_subs=False)
susy.fierz_expand_2index(ex, basis, indices, to_perform_subs=False)
print('\n', ex)
cdb.factor_in(ex, Ex('u, v, w, x, y', False))
print('\n', ex)

# SECTION 5
## Section 5.2
__cdbkernel__ = susy_env()

### Sec. 5.2: Code Block #1
ex = Ex(r'(\Gamma_{a b} \Gamma^{c d})_{\alpha}^{\alpha}')
evaluate_traces(ex)
print('\n', ex)

### Sec. 5.2: Code Block #2
ex = Ex(r'(-1/32) (\Gamma_{a} \Gamma_{b} \Gamma_{c} \Gamma_{d} \Gamma_{e} \Gamma_{f} \Gamma_{g} \Gamma_{h} \Gamma_{i} \Gamma_{j} \Gamma_{k})^{\alpha}_{\alpha}')
evaluate_traces(ex)
eliminate_kronecker(ex)
canonicalise(ex)
collect_terms(ex)
print('\n', ex)

## Section 5.5
__cdbkernel__ = susy_env()

### Sec. 5.5: Code Block #1
ex = Ex(r'\delta_{\alpha}^{\eta} (\Gamma_{a})^{\rho \gamma} - \delta_{\alpha}^{\rho} (\Gamma_{a})^{\eta \gamma}')
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
fierz_expand(ex, basis, [r'_{\alpha}', r'_{\gamma}'], [r'_{\eta}', r'_{\rho}'])
print('\n', ex)

### Sec. 5.5: Code Block #2
ex = Ex(r'(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}')
fierz_expand(ex, basis, [r'_{\alpha}', r'_{\gamma}'], [r'_{\gamma}', r'_{\eta}'])
print('\n', ex)

## Section 5.6
### Sec. 5.6: Code Block #1
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

### Sec. 5.6: Code Block #2
__cdbkernel__ = susy_env()
ex = Ex(r'''(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}''')
fierz_expand_2index(ex, basis, [r'_{\alpha}', r'_{\beta}'])
print('\n', ex)

## Section 5.7
### Sec. 5.7: Code Block #1
__cdbkernel__ = susy_env()
Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''C'''), Ex(r'''D{#}, \partial{#}'''))

ex = Ex(r'''\partial_{a}(A B) \partial_{b c}(C) + \partial_{a b c}(A B C)''')
fourier(ex, [Ex('A'), Ex('B'), Ex('C')])
canonicalise(ex)
print('\n', ex)

### Sec. 5.7: Code Block #2
inverse_fourier(ex)
canonicalise(ex)
print('\n', ex)

## Section 5.8
### Sec. 5.8: Code Block #1
__cdbkernel__ = susy_env(D=4)
AntiSymmetric(Ex(r'''B{#}'''))
Depends(Ex(r'''B{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''D{#}, \partial{#}'''))
field = Ex(r'B_{a b}')
gauge_trans = Ex(r'\partial_{a}(\zeta_{b}) - \partial_{b}(\zeta_{a})')
ex = Ex(r'(u + v) (\Gamma^{a})_{\alpha \beta} \partial_{a}(B_{b c}) - v (\Gamma^{a})_{\alpha \beta} \partial_{b}(B_{a c}) + (u + v) (\Gamma^{a})_{\alpha \beta} \partial_{c}(B_{a b})')
print('\n', find_non_gauge_inv(ex, field, [field], gauge_trans))

## Section 5.9
### Sec. 5.9: Code Block #1
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

### Sec. 5.9: Code Block #2
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

### Sec. 5.9: Code Block #3
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

## Section 5.10
__cdbkernel__ = susy_env(D=4)

### Sec. 5.10: Code Block #1
subs = sym_proj_decomp(3)
print('\n', subs)

### Sec. 5.10: Code Block #2
subs = sym_proj_decomp(3/2)
print('\n', subs)

### Sec. 5.10: Code Block #3
prop = sym_prop(2)
print('\n', prop)

### Sec. 5.10: Code Block #4
prop = sym_prop(3/2)
print('\n', prop)

### Sec. 5.10: Code Block #5
ex = sym_prop(3/2)
split_gamma(ex, on_back=True)
split_gamma(ex, on_back=True)
evaluate(ex, to_join_gamma=False)
print('\n', ex)

### Sec. 5.10: Code Block #6
__cdbkernel__ = susy_env(D=11)
prop = rarita_schwinger_prop()
print('\n', prop)

## Section 5.11
### Sec. 5.11: Equation 2.11.5
__cdbkernel__ = susy_env(D=4)
print(evaluate(rarita_schwinger_prop()))

### Sec. 5.11: Code Block #1 (Axial-Vector)
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

### Sec. 5.11: Code Block #2 (Matter-Gravitino)
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

### Sec. 5.11: Code Block #3 (Supergravity)
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

## Section 5.12
### Sec. 5.12: Code Block #1
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

## Section 5.12
### Sec. 5.12: Code Block #1
__cdbkernel__ = susy_env(D = 4)
Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
AntiCommuting(Ex(r'''\Psi, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\Psi{#}}{#}, D{#}'''))

fields = [Ex('A'), Ex('B'), Ex(r'(\Psi)_{\gamma}')]
susy = r'''D_{\alpha}(A) -> u (\Psi)_{\alpha}, D_{\alpha}(B) -> v (\Gamma')_{\alpha}^{\beta} (\Psi)_{\beta}, D_{\alpha}((\Psi)_{\beta}) -> w (\Gamma^{a})_{\alpha \beta} \partial_{a}(A) + x (\Gamma' \Gamma^{a})_{\alpha \beta} \partial_{a}(B)'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
subs = Ex('u -> 1, v -> I, w -> I, x -> -1', False)
indices = [r'_{\alpha}', r'_{\beta}']
T = holoraumy(fields, susy, basis, subs, indices)
print('\n')
pretty_print(T)

### Sec. 5.12: Code Block #2
__cdbkernel__ = susy_env(D = 4)
Depends(Ex(r'''h{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta{#}'''), Ex(r'''\partial{#}'''))
Depends(Ex(r'''\indexbracket{\zeta{#}}{#}'''), Ex(r'''\partial{#}'''))
Symmetric(Ex(r'''h{#}'''))
AntiCommuting(Ex(r'''\Psi, D{#}'''))
AntiCommuting(Ex(r'''\indexbracket{\Psi{#}}{#}, D{#}'''))

fields = [Ex('h_{a b}'), Ex(r'(\Psi_{a})_{\gamma}')]
susy = r'''D_{\alpha}(h_{a b}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\Psi_{b})_{\beta} + u (\Gamma_{b})_{\alpha}^{\beta} (\Psi_{a})_{\beta}, D_{\alpha}((\Psi_{a})_{\beta}) -> v (\Gamma^{b c})_{\alpha \beta} \partial_{b}(h_{c a})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
subs = Ex('u -> 1/2, v -> -I', False)
indices = [r'_{\alpha}', r'_{\beta}']
holo_data = holoraumy(fields, susy, basis, subs, indices)
print('\n')
pretty_print(holo_data)

### Sec. 5.12: Code Block #3
ex1 = Ex(r'D_{\alpha}(D_{\beta}((\Psi_{a})_{\gamma})) - D_{\beta}(D_{\alpha}((\Psi_{a})_{\gamma}))')
susy_expand(ex1, susy)
substitute(ex1, subs)
distribute(ex1)
evaluate(ex1, to_perform_subs=False)
fierz_expand_2index(ex1, basis, indices)
ex2a = Ex(r'''2 I (-(3/4) \delta_{a b} (\Gamma' \Gamma_{c})_{\alpha \beta} (\Gamma')_{\gamma}^{\eta} + (3/4) \delta_{a c} (\Gamma' \Gamma_{b})_{\alpha \beta} (\Gamma')_{\gamma}^{\eta} - (1/8) I \epsilon_{a b c d} (\Gamma' \Gamma^{d})_{\alpha \beta} C_{\gamma}^{\eta} + (1/8) (\Gamma' \Gamma_{b})_{\alpha \beta} (\Gamma' \Gamma_{c a})_{\gamma}^{\eta} - (1/8) (\Gamma' \Gamma_{c})_{\alpha \beta} (\Gamma' \Gamma_{b a})_{\gamma}^{\eta}) \partial^{b}((\Psi^{c})_{\eta})''')
ex2b = Ex(r'''I (-C_{\alpha \beta} (\Gamma_{a})_{\gamma}^{\eta} - (\Gamma')_{\alpha \beta} (\Gamma' \Gamma_{a})_{\gamma}^{\eta} + (1/2) (\Gamma' \Gamma^{d})_{\alpha \beta} (\Gamma' \Gamma_{d} \Gamma_{a})_{\gamma}^{\eta} - (1/4) (\Gamma' \Gamma_{a})_{\alpha \beta} (\Gamma')_{\gamma}^{\eta}) (\Gamma^{e f})_{\eta}^{\rho} \partial_{e}((\Psi_{f})_{\rho}) + (1/4) (5 C_{\alpha \beta} C_{\gamma}^{\eta} + 5 (\Gamma')_{\alpha \beta} (\Gamma')_{\gamma}^{\eta} - 2 (\Gamma' \Gamma^{d})_{\alpha \beta} (\Gamma' \Gamma_{d})_{\gamma}^{\eta}) \epsilon_{a}^{e f g} (\Gamma' \Gamma_{e})_{\eta}^{\rho} \partial_{f}((\Psi_{g})_{\rho})''')
ex2c = Ex(r'''(3/4) I (C_{\alpha \beta} (\Gamma^{b})_{\gamma}^{\eta} + (\Gamma')_{\alpha \beta} (\Gamma' \Gamma^{b})_{\gamma}^{\eta} + (\Gamma' \Gamma^{b})_{\alpha \beta} (\Gamma')_{\gamma}^{\eta} - (1/3) (\Gamma' \Gamma_{c})_{\alpha \beta} (\Gamma' \Gamma^{c b})_{\gamma}^{\eta}) \partial_{a}((\Psi_{b})_{\eta})''')
ex = ex1 + ex2a + ex2b + ex2c
evaluate(ex)
evaluate(ex)
print('\n', ex)
indexbracket_hex(ex, to_decompose_product=True)
evaluate(ex)
print('\n', ex)

### Sec. 5.12: Code Block #4
__cdbkernel__ = susy_env(D = 4)
Depends(Ex(r'''E{#}'''), Ex(r'''\partial{#}'''))
AntiSymmetric(Ex(r'''B{#}'''))
ex = Ex(r'''-  1/4 \epsilon_{a}^{b c d} A^{e} B_{b c} \partial_{e}(E_{d}) +  1/4 \epsilon_{a}^{b c d} A^{e} B_{b c} \partial_{d}(E_{e}) +  1/2 \epsilon_{a}^{b c d} A_{e} B_{b}^{e} \partial_{c}(E_{d}) +  1/4 \epsilon^{b c d e} A_{a} B_{b c} \partial_{d}(E_{e})''')
decompose_product(ex)
canonicalise(ex)
collect_terms(ex)
print('\n', ex)