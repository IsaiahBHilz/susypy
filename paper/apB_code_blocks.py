import sys
sys.path.append('../')

from cadabra2 import Ex
from susypy import fierz_expand, susy_env

__cdbkernel__ = susy_env()

# APPENDIX B
## Appendix B.1
### Ap. B.1: Code Block #1:
fierz1 = fierz_expand(Ex(r'(\Gamma^{d e})_{\eta \alpha} (\Gamma_{d e a b c})_{\beta}^{\gamma} + (\Gamma^{d e})_{\eta \beta} (\Gamma_{d e a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz1)

### Ap. B.1: Code Block #2:
fierz2 = fierz_expand(Ex(r'(\Gamma_{a b})_{\eta \alpha} (\Gamma^{a b c})_{\beta}^{\gamma} + (\Gamma_{a b})_{\eta \beta} (\Gamma^{a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz2)

### Ap. B.1: Code Block #3:
fierz3 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{a c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{a c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz3)

### Ap. B.1: Code Block #4:
fierz4 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{b c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz4)

### Ap. B.1: Code Block #5:
fierz5 = fierz_expand(Ex(r'(\Gamma_{d a})_{\eta \alpha} (\Gamma^{d b c})_{\beta}^{\gamma} + (\Gamma_{d a})_{\eta \beta} (\Gamma^{d b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz5)

## Appendix B.2
### Ap. B.2: Code Block #1:
fierz6 = fierz_expand(Ex(r'(\Gamma^{d e})_{\eta \alpha} (\Gamma_{d e a b c})_{\beta}^{\gamma} - (\Gamma^{d e})_{\eta \beta} (\Gamma_{d e a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz6)

## Appendix B.7
### Ap. B.2: Code Block #2:
fierz7 = fierz_expand(Ex(r'(\Gamma_{a b})_{\eta \alpha} (\Gamma^{a b c})_{\beta}^{\gamma} - (\Gamma_{a b})_{\eta \beta} (\Gamma^{a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz7)

## Appendix B.8
### Ap. B.2: Code Block #3:
fierz8 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{a c})_{\beta}^{\gamma} - (\Gamma_{a})_{\eta \beta} (\Gamma^{a c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz8)

## Appendix B.9
### Ap. B.2: Code Block #4:
fierz9 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{b c})_{\beta}^{\gamma} - (\Gamma_{a})_{\eta \beta} (\Gamma^{b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz9)

## Appendix B.10
### Ap. B.2: Code Block #5:
fierz10 = fierz_expand(Ex(r'(\Gamma_{d a})_{\eta \alpha} (\Gamma^{d b c})_{\beta}^{\gamma} - (\Gamma_{d a})_{\eta \beta} (\Gamma^{d b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz10)