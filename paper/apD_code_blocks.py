import sys
sys.path.append('../')

from cadabra2 import Ex, collect_terms
from susypy import indexbracket_hex, susy_env

__cdbkernel__ = susy_env()

# Appendix D
## Appendix D.0
## Ap. D.0: Code Block #1
ex = Ex(r'(\Gamma^{a b})_{\alpha \beta} (\Theta_{b a})_{\gamma \eta} + (\Gamma^{a b})_{\alpha \beta} (\Theta_{a b})_{\gamma \eta}')
indexbracket_hex(ex)
collect_terms(ex)
print('\n', ex)