import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric, distribute, collect_terms, sort_product, asym

from susypy import susy_env, susy_solve, holoraumy, make_action_susy_inv, susy_expand, evaluate, substitute, spinor_combine, indexbracket_hex, rarita_schwinger_prop
from susypy.tools import logger, timer, pretty_print


#__cdbkernel__ = susy_env(spinor_indices = [r'\alpha', r'\beta', r'\gamma', r'\eta', r'\rho', r'\sigma', r'\zeta'])

__cdbkernel__ = susy_env(D=4)

print(evaluate(rarita_schwinger_prop()))

"""
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

bosons = [Ex('h_{a b}'), Ex('A_{a b c}')]
fermions = [Ex(r'(\Psi_{a})^{\gamma}')]
gauge_transs = [Ex(r'\partial_{a}((\zeta)^{\gamma})')]
susy = r'D_{\alpha}(h_{a b}) -> u ((\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}), D_{\alpha}((\Psi_{b})^{\beta}) -> 2 v \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + x (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) + y (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> 2 z (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - 2 z (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + 2 z (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}, c_{a b c} -> \partial_{a}(h_{b c}) - \partial_{b}(h_{a c})'
L = Ex(r'D_{\gamma}(l (-1/4 c^{a b c} c_{a b c} + 1/2 c^{a b c} c_{c a b} + c^{a}_{b}^{b} c_{a c}^{c}) + m (1/12) (\Psi_{a})^{\alpha} (\Gamma^{a b c})_{\alpha \beta} \partial_{b}((\Psi_c)^{\beta}) + n (1/48) (\partial_{a}(A_{b c d})-\partial_{a}(A_{b d c})-\partial_{a}(A_{c b d}) + \partial_{a}(A_{c d b}) + \partial_{a}(A_{d b c})-\partial_{a}(A_{d c b})-\partial_{b}(A_{a c d}) + \partial_{b}(A_{a d c}) + \partial_{b}(A_{c a d})-\partial_{b}(A_{c d a})-\partial_{b}(A_{d a c}) + \partial_{b}(A_{d c a}) + \partial_{c}(A_{a b d})-\partial_{c}(A_{a d b})-\partial_{c}(A_{b a d}) + \partial_{c}(A_{b d a}) + \partial_{c}(A_{d a b})-\partial_{c}(A_{d b a})-\partial_{d}(A_{a b c}) + \partial_{d}(A_{a c b}) + \partial_{d}(A_{b a c})-\partial_{d}(A_{b c a})-\partial_{d}(A_{c a b}) + \partial_{d}(A_{c b a})) (\partial^{a}(A^{b c d})-\partial^{a}(A^{b d c})-\partial^{a}(A^{c b d}) + \partial^{a}(A^{c d b}) + \partial^{a}(A^{d b c})-\partial^{a}(A^{d c b})-\partial^{b}(A^{a c d}) + \partial^{b}(A^{a d c}) + \partial^{b}(A^{c a d})-\partial^{b}(A^{c d a})-\partial^{b}(A^{d a c}) + \partial^{b}(A^{d c a}) + \partial^{c}(A^{a b d})-\partial^{c}(A^{a d b})-\partial^{c}(A^{b a d}) + \partial^{c}(A^{b d a}) + \partial^{c}(A^{d a b})-\partial^{c}(A^{d b a})-\partial^{d}(A^{a b c}) + \partial^{d}(A^{a c b}) + \partial^{d}(A^{b a c})-\partial^{d}(A^{b c a})-\partial^{d}(A^{c a b}) + \partial^{d}(A^{c b a})))')
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
consts = ['u', 'v', 'x', 'y', 'z', 'l', 'm', 'n']
indices = [r'_{\alpha}', r'_{\beta}']

ex = Ex(r'\partial_{a}^{b}(A_{c d e}) \Gamma_{b}^{c d e}')
print(evaluate(Ex('120')*asym(ex, Ex('_{a}, _{c}, _{d}, _{e}, _{b}', False))))
"""

"""
susy_sol, susy_data = timer(susy_solve)(bosons, fermions, gauge_transs, susy, basis, consts, indices)
ex_main = susy_data[r'\Psi']['desired_terms'] + susy_data[r'\Psi']['gauge_terms'] + susy_data[r'\Psi']['undesired_terms'] + susy_data[r'\Psi']['lorentz_proper_terms']
print(ex_main)

#print(evaluate(Ex(r'C_{\eta \rho} (\Gamma^{a b})_{\alpha \beta} C^{\gamma \beta} (\Gamma_{c d})^{\alpha}_{\gamma} (\Gamma^{c d})^{\rho}_{\sigma}')))

#ex1 = Ex(r'(-11/4 u v + 63 x z) (\Gamma^{b})_{\alpha \beta} \partial_{b}((\Psi_{a})^{\gamma}) + (-1/2 u v - 27 x z - 63 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a})_{\eta}^{\gamma} \partial_{b}((\Psi_{c})^{\eta}) + (-9/4 u v - 27 x z) (\Gamma_{b}^{c})_{\alpha \beta} (\Gamma^{b})_{\eta}^{\gamma} \partial_{c}((\Psi_{a})^{\eta}) + (1/4 u v + 3/2 x z - 9/2 y z) (\Gamma_{a b c}^{d e})_{\alpha \beta} (\Gamma^{b c})_{\eta}^{\gamma} \partial_{d}((\Psi_{e})^{\eta}) + (1/12 u v + 3/2 x z + 3/2 y z) (\Gamma_{b c d}^{e f})_{\alpha \beta} (\Gamma_{a}^{b c d})_{\eta}^{\gamma} \partial_{e}((\Psi_{f})^{\eta}) + (-1/32 u v - 3/8 x z) (\Gamma_{b c d e}^{f})_{\alpha \beta} (\Gamma^{b c d e})_{\eta}^{\gamma} \partial_{f}((\Psi_{a})^{\eta}) + (-1/8 u v - 27 x z - 15 y z) (\Gamma_{a})_{\alpha \beta} R^{\gamma} + (1/8 u v - 6 x z - 15 y z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{a b})^{\gamma \zeta} R_{\zeta} + (-1/8 u v + 15 x z + 6 y z) (\Gamma_{a}^{b})_{\alpha \beta} (\Gamma_{b})^{\gamma \zeta} R_{\zeta} + (1/16 u v + 3 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a b c})^{\gamma \zeta} R_{\zeta} + (-1/192 u v + 1/8 x z + 1/8 y z) (\Gamma_{a}^{b c d e})_{\alpha \beta} (\Gamma_{b c d e})^{\gamma \zeta} R_{\zeta} + (1/115200 u v + 1/4800 y z) \epsilon_{a b c d e f g h i j k} (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma^{g h i j k}^{\gamma \zeta}) R_{\zeta} + (1/4 u v - 21 x z - 42 y z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{a})^{\gamma \zeta} (E_{b})_{\zeta} + (-5/4 u v + 21 x z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{b})^{\gamma \zeta} (E_{a})_{\zeta} + (1/4 u v - 42 x z - 21 y z) (\Gamma_{a}^{c})_{\alpha \beta} (E_{c})^{\gamma} + (1/4 u v + 6 x z + 21 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a b})^{\gamma \zeta} (E_{c})_{\zeta} + (1/2 u v + 3 x z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{b c})^{\gamma \zeta} (E_{a})_{\zeta} + (-1/24 u v + 1/2 x z + y z) (\Gamma_{a}^{b c d e})_{\alpha \beta} (\Gamma_{b c d})^{\gamma \zeta} (E_{e})_{\zeta} + (1/96 u v + 1/8 x z + 1/4 y z) (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma_{a b c d e})^{\gamma \zeta} (E_{f})_{\zeta} + (-1/480 u v - 1/40 x z) (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma_{b c d e f})^{\gamma \zeta} (E_{a})_{\zeta} + (3/2 u v - 63 x z) (\Gamma^{b})_{\alpha \beta} \partial_{a}((\Psi_{b})^{\gamma}) - u v (\Gamma^{b})_{\alpha \beta} (\Gamma_{b}^{c})_{\eta}^{\gamma} \partial_{a}((\Psi_{c})^{\eta}) + (5/4 u v + 27 x z) (\Gamma_{b}^{c})_{\alpha \beta} (\Gamma^{b})_{\eta}^{\gamma} \partial_{a}((\Psi_{c})^{\eta}) - (3/8 u v) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{b c}^{d})_{\eta}^{\gamma} \partial_{a}((\Psi_{d})^{\eta}) + (1/48 u v + 3/8 x z) (\Gamma_{b c d e}^{f})_{\alpha \beta} (\Gamma^{b c d e})_{\eta}^{\gamma} \partial_{a}((\Psi_{f})^{\eta})')

ex2 = Ex(r'(-11/4 u v + 63 x z) (\Gamma^{b})_{\alpha \beta} \partial_{b}((\Psi_{a})^{\gamma}) + (-1/2 u v - 27 x z - 63 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a})_{\eta}^{\gamma} \partial_{b}((\Psi_{c})^{\eta}) + (-9/4 u v - 27 x z) (\Gamma_{b}^{c})_{\alpha \beta} (\Gamma^{b})_{\eta}^{\gamma} \partial_{c}((\Psi_{a})^{\eta}) + (1/4 u v + 3/2 x z - 9/2 y z) (\Gamma_{a b c}^{d e})_{\alpha \beta} (\Gamma^{b c})_{\eta}^{\gamma} \partial_{d}((\Psi_{e})^{\eta}) + (1/12 u v + 3/2 x z + 3/2 y z) (\Gamma_{b c d}^{e f})_{\alpha \beta} (\Gamma_{a}^{b c d})_{\eta}^{\gamma} \partial_{e}((\Psi_{f})^{\eta}) + (-1/32 u v - 3/8 x z) (\Gamma_{b c d e}^{f})_{\alpha \beta} (\Gamma^{b c d e})_{\eta}^{\gamma} \partial_{f}((\Psi_{a})^{\eta}) + (-1/8 u v - 27 x z - 15 y z) (\Gamma_{a})_{\alpha \beta} R^{\gamma} + (1/8 u v - 6 x z - 15 y z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{a b})^{\gamma \zeta} R_{\zeta} + (-1/8 u v + 15 x z + 6 y z) (\Gamma_{a}^{b})_{\alpha \beta} (\Gamma_{b})^{\gamma \zeta} R_{\zeta} + (1/16 u v + 3 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a b c})^{\gamma \zeta} R_{\zeta} + (-1/192 u v + 1/8 x z + 1/8 y z) (\Gamma_{a}^{b c d e})_{\alpha \beta} (\Gamma_{b c d e})^{\gamma \zeta} R_{\zeta} + (1/115200 u v + 1/4800 y z) \epsilon_{a b c d e f g h i j k} (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma^{g h i j k}^{\gamma \zeta}) R_{\zeta} + (1/4 u v - 21 x z - 42 y z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{a})^{\gamma \zeta} (E_{b})_{\zeta} + (-5/4 u v + 21 x z) (\Gamma^{b})_{\alpha \beta} (\Gamma_{b})^{\gamma \zeta} (E_{a})_{\zeta} + (-1/4 u v + 42 x z + 21 y z) (\Gamma_{a}^{c})_{\alpha \beta} (E_{c})^{\gamma} + (1/4 u v + 6 x z + 21 y z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{a b})^{\gamma \zeta} (E_{c})_{\zeta} + (1/2 u v + 3 x z) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{b c})^{\gamma \zeta} (E_{a})_{\zeta} + (-1/24 u v + 1/2 x z + y z) (\Gamma_{a}^{b c d e})_{\alpha \beta} (\Gamma_{b c d})^{\gamma \zeta} (E_{e})_{\zeta} + (1/96 u v + 1/8 x z + 1/4 y z) (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma_{a b c d e})^{\gamma \zeta} (E_{f})_{\zeta} + (-1/480 u v - 1/40 x z) (\Gamma^{b c d e f})_{\alpha \beta} (\Gamma_{b c d e f})^{\gamma \zeta} (E_{a})_{\zeta} + (3/2 u v - 63 x z) (\Gamma^{b})_{\alpha \beta} \partial_{a}((\Psi_{b})^{\gamma}) - u v (\Gamma^{b})_{\alpha \beta} (\Gamma_{b}^{c})_{\eta}^{\gamma} \partial_{a}((\Psi_{c})^{\eta}) + (5/4 u v + 27 x z) (\Gamma_{b}^{c})_{\alpha \beta} (\Gamma^{b})_{\eta}^{\gamma} \partial_{a}((\Psi_{c})^{\eta}) - (3/8 u v) (\Gamma^{b c})_{\alpha \beta} (\Gamma_{b c}^{d})_{\eta}^{\gamma} \partial_{a}((\Psi_{d})^{\eta}) + (1/48 u v + 3/8 x z) (\Gamma_{b c d e}^{f})_{\alpha \beta} (\Gamma^{b c d e})_{\eta}^{\gamma} \partial_{a}((\Psi_{f})^{\eta})')

subs = Ex(r'R_{\alpha} -> (\Gamma^{b c})_{\alpha \eta} \partial_{b}((\Psi_{c})^{\eta}), R^{\alpha} -> (\Gamma^{b c})^{\alpha}_{\eta} \partial_{b}((\Psi_{c})^{\eta}), (E_{c})_{\alpha} -> (\Gamma^{b})_{\alpha \eta} (\partial_{b}((\Psi_{c})^{\eta}) - \partial_{c}((\Psi_{b})^{\eta})), (E_{c})^{\alpha} -> (\Gamma^{b})^{\alpha}_{\eta} (\partial_{b}((\Psi_{c})^{\eta}) - \partial_{c}((\Psi_{b})^{\eta}))', False)

#substitute(ex1, subs)
substitute(ex2, subs)

#diff1 = ex1-ex_main

diff2 = ex2-ex_main

distribute(diff2)
distribute(diff2)
distribute(diff2)
distribute(diff2)
#canonicalise(diff2)
spinor_combine(diff2)
sort_product(diff2)
#indexbracket_hex(diff2)
collect_terms(diff2)
print('\n', diff2)

print('\n', timer(evaluate)(diff2, to_join_gamma=False, to_perform_subs=False))
#print('\n', timer(evaluate)(ex2-ex_main))
"""