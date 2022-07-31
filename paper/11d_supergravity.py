import sys
sys.path.append('../')

from cadabra2 import AntiCommuting, AntiSymmetric, Depends, Ex, Symmetric

from susypy import susy_env, susy_solve, holoraumy, make_action_susy_inv, susy_expand, evaluate
from susypy.tools import logger, timer, pretty_print

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

fields = [Ex('h_{a b}'), Ex('A_{a b c}'), Ex(r'(\Psi_{a})^{\gamma}')]
gauge_transs = [Ex(r'\partial_{a}(\zeta_{b}) + \partial_{b}(\zeta_{a})'), Ex(r'\partial_{a}(\zeta_{b c}) - \partial_{a}(\zeta_{c b}) + \partial_{c}(\zeta_{a b}) - \partial_{c}(\zeta_{b a}) + \partial_{b}(\zeta_{c a}) - \partial_{b}(\zeta_{a c})'), Ex(r'\partial_{a}((\zeta)^{\gamma})')]
susy = r'D_{\alpha}(h_{a b}) -> u ((\Gamma_{a})_{\alpha \beta} (\Psi_{b})^{\beta} + (\Gamma_{b})_{\alpha \beta} (\Psi_{a})^{\beta}), D_{\alpha}((\Psi_{b})^{\beta}) -> 2 v \partial_{e}(h_{b d}) (\Gamma^{d e})_{\alpha}^{\beta} + x (\Gamma_{b} \Gamma^{c d e f})_{\alpha}^{\beta} \partial_{c}(A_{d e f}) + y (\Gamma^{c d e f} \Gamma_{b})_{\alpha}^{\beta} \partial_{c}(A_{d e f}), D_{\alpha}(A_{b c d}) -> 2 z (\Gamma_{b c})_{\alpha \beta} (\Psi_{d})^{\beta} - 2 z (\Gamma_{b d})_{\alpha \beta} (\Psi_{c})^{\beta} + 2 z (\Gamma_{c d})_{\alpha \beta} (\Psi_{b})^{\beta}, c_{a b c} -> \partial_{a}(h_{b c}) - \partial_{b}(h_{a c})'
L = Ex(r'D_{\gamma}(l (-1/4 c^{a b c} c_{a b c} + 1/2 c^{a b c} c_{c a b} + c^{a}_{b}^{b} c_{a c}^{c}) + m (1/12) (\Psi_{a})^{\alpha} (\Gamma^{a b c})_{\alpha \beta} \partial_{b}((\Psi_c)^{\beta}) + n (1/48) (\partial_{a}(A_{b c d})-\partial_{a}(A_{b d c})-\partial_{a}(A_{c b d}) + \partial_{a}(A_{c d b}) + \partial_{a}(A_{d b c})-\partial_{a}(A_{d c b})-\partial_{b}(A_{a c d}) + \partial_{b}(A_{a d c}) + \partial_{b}(A_{c a d})-\partial_{b}(A_{c d a})-\partial_{b}(A_{d a c}) + \partial_{b}(A_{d c a}) + \partial_{c}(A_{a b d})-\partial_{c}(A_{a d b})-\partial_{c}(A_{b a d}) + \partial_{c}(A_{b d a}) + \partial_{c}(A_{d a b})-\partial_{c}(A_{d b a})-\partial_{d}(A_{a b c}) + \partial_{d}(A_{a c b}) + \partial_{d}(A_{b a c})-\partial_{d}(A_{b c a})-\partial_{d}(A_{c a b}) + \partial_{d}(A_{c b a})) (\partial^{a}(A^{b c d})-\partial^{a}(A^{b d c})-\partial^{a}(A^{c b d}) + \partial^{a}(A^{c d b}) + \partial^{a}(A^{d b c})-\partial^{a}(A^{d c b})-\partial^{b}(A^{a c d}) + \partial^{b}(A^{a d c}) + \partial^{b}(A^{c a d})-\partial^{b}(A^{c d a})-\partial^{b}(A^{d a c}) + \partial^{b}(A^{d c a}) + \partial^{c}(A^{a b d})-\partial^{c}(A^{a d b})-\partial^{c}(A^{b a d}) + \partial^{c}(A^{b d a}) + \partial^{c}(A^{d a b})-\partial^{c}(A^{d b a})-\partial^{d}(A^{a b c}) + \partial^{d}(A^{a c b}) + \partial^{d}(A^{b a c})-\partial^{d}(A^{b c a})-\partial^{d}(A^{c a b}) + \partial^{d}(A^{c b a})))')
basis = [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')]
consts = ['u', 'v', 'x', 'y', 'z', 'l', 'm', 'n']
indices = [r'_{\alpha}', r'_{\beta}']

susy_sol, susy_data = logger('susy_solve_test')(timer(susy_solve))(fields, gauge_transs, susy, basis, consts, indices)
pretty_print(susy_data)

act_sol = timer(make_action_susy_inv)(S, susy, consts, [r'\Psi'])
print('\n', act_sol)

subs = Ex('u -> 1, v -> (-1/8) I, x -> 1/48 I, y -> (-1/144) I, z -> 1/2', False)
holo_data = logger('holo_test')(timer(holoraumy))(fields, susy, basis, subs, indices)
pretty_print(holo_data)