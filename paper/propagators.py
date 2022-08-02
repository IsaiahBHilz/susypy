import sys
sys.path.append('../')

from cadabra2 import Depends, Ex, Symmetric, AntiSymmetric
from susypy import susy_env, susy_solve_propagator, rarita_schwinger_prop

__cdbkernel__ = susy_env(D=4)

# Chiral
Depends(Ex(r'''A'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''B'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\Psi{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\Psi{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))

bosons = [Ex('A'), Ex('B')]
fermions = [Ex(r'(\Psi)_{\gamma}')]
fermion_propagators = [Ex(r'I (k1^{a} k1_{a})**(-1) (\Gamma^{b})_{\alpha \beta} k1_{b}')]
susy = r'''D_{\alpha}(A) -> u (\Psi)_{\alpha}, D_{\alpha}(B) -> v (\Gamma')_{\alpha}^{\beta} (\Psi)_{\beta}, D_{\alpha}((\Psi)_{\beta}) -> w (\Gamma^{a})_{\alpha \beta} \partial_{a}(A) + x (\Gamma' \Gamma^{a})_{\alpha \beta} \partial_{a}(B)'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v', 'w', 'x']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

# Vector
Depends(Ex(r'''A{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))

bosons = [Ex('A_{a}')]
fermions = [Ex(r'(\lambda)_{\gamma}')]
fermion_propagators = [Ex(r'I (k1^{a} k1_{a})**(-1) (\Gamma^{b})_{\alpha \beta} k1_{b}')]
susy = r'''D_{\alpha}(A_{a}) -> u (\Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma^{a b})_{\alpha \beta} \partial_{a}(A_{b})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

# Axial-Vector
Depends(Ex(r'''U{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\lambda{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\indexbracket{\lambda{#}}{#}'''), Ex(r'''D{#}, \partial{#}'''))
Depends(Ex(r'''\zeta'''), Ex(r'''\partial{#}'''))

bosons = [Ex('U_{a}')]
fermions = [Ex(r'(\lambda)_{\gamma}')]
fermion_propagators = [Ex(r'I (k1^{a} k1_{a})**(-1) (\Gamma^{b})_{\alpha \beta} k1_{b}')]
susy = r'''D_{\alpha}(U_{a}) -> u (\Gamma' \Gamma_{a})_{\alpha}^{\beta} (\lambda)_{\beta}, D_{\alpha}((\lambda)_{\beta}) -> v (\Gamma' \Gamma^{a b})_{\alpha \beta} \partial_{a}(U_{b})'''
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'''(\Gamma')_{\alpha \beta}'''), Ex(r'''(\Gamma' \Gamma^{a})_{\alpha \beta}''')]
consts = ['u', 'v']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

# Matter-Gravitino
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
#consts = ['v', 'w', 'x', 'y']
indices = [r'_{\alpha}', r'_{\beta}']
sol = susy_solve_propagator(bosons, fermions, fermion_propagators, susy, basis, consts, indices, comm_coef=2)
print('\n', sol)

# Supergravity
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