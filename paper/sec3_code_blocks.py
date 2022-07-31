import sys
sys.path.append('../')

from cadabra2 import Ex
from susypy import fierz_expand, susy_env

__cdbkernel__ = susy_env()

# SECTION 3
## Section 3.1
### Sec. 3.1: Code Block #1:
fierz1 = fierz_expand(Ex(r'(\Gamma^{d e})_{\eta \alpha} (\Gamma_{d e a b c})_{\beta}^{\gamma} + (\Gamma^{d e})_{\eta \beta} (\Gamma_{d e a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz1)

## Section 3.2
### Sec. 3.2: Code Block #1:
fierz2 = fierz_expand(Ex(r'(\Gamma_{a b})_{\eta \alpha} (\Gamma^{a b c})_{\beta}^{\gamma} + (\Gamma_{a b})_{\eta \beta} (\Gamma^{a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz2)

## Section 3.3
### Sec. 3.3: Code Block #1:
fierz3 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{a c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{a c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz3)

## Section 3.4
### Sec. 3.4: Code Block #1:
fierz4 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{b c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz4)

## Section 3.5
### Sec. 3.5: Code Block #1:
fierz5 = fierz_expand(Ex(r'(\Gamma_{d a})_{\eta \alpha} (\Gamma^{d b c})_{\beta}^{\gamma} + (\Gamma_{d a})_{\eta \beta} (\Gamma^{d b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz5)

## Section 3.6
### Sec. 3.6: Code Block #1:
fierz6 = fierz_expand(Ex(r'(\Gamma^{d e})_{\eta \alpha} (\Gamma_{d e a b c})_{\beta}^{\gamma} - (\Gamma^{d e})_{\eta \beta} (\Gamma_{d e a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz6)

## Section 3.7
### Sec. 3.7: Code Block #1:
fierz7 = fierz_expand(Ex(r'(\Gamma_{a b})_{\eta \alpha} (\Gamma^{a b c})_{\beta}^{\gamma} - (\Gamma_{a b})_{\eta \beta} (\Gamma^{a b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz7)

## Section 3.8
### Sec. 3.8: Code Block #1:
fierz8 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{a c})_{\beta}^{\gamma} - (\Gamma_{a})_{\eta \beta} (\Gamma^{a c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz8)

## Section 3.9
### Sec. 3.9: Code Block #1:
fierz9 = fierz_expand(Ex(r'(\Gamma_{a})_{\eta \alpha} (\Gamma^{b c})_{\beta}^{\gamma} - (\Gamma_{a})_{\eta \beta} (\Gamma^{b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz9)

## Section 3.10
### Sec. 3.10: Code Block #1:
fierz10 = fierz_expand(Ex(r'(\Gamma_{d a})_{\eta \alpha} (\Gamma^{d b c})_{\beta}^{\gamma} - (\Gamma_{d a})_{\eta \beta} (\Gamma^{d b c})_{\alpha}^{\gamma}'), [Ex(r'C_{\alpha \beta}'), Ex(r'(\Gamma^{a})_{\alpha \beta}'), Ex(r'(\Gamma^{a b})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')], [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
print('\n', fierz10)