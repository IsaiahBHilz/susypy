import unittest
import sys
sys.path.append('../')

from cadabra2 import Ex
from susypy import fierz_expand, fierz_expand_2index, susy_env, load

class TestFierzExpand(unittest.TestCase):

	def test_fierz1(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex1  = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} - \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
		ex1  = fierz_expand(ex1, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex1p = load('test_exs/fierz1.p')

		self.assertEqual(ex1, ex1p)

	def test_fierz2(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex2  = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} + \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
		ex2  = fierz_expand(ex2, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex2p = load('test_exs/fierz2.p')

		self.assertEqual(ex2, ex2p)

	def test_fierz3(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex3  = Ex(r'(\Gamma^{a b})_{\alpha}^{\eta} (\Gamma_{a b})_{\beta}^{\gamma} + (\Gamma^{a b})_{\beta}^{\eta} (\Gamma_{a b})_{\alpha}^{\gamma}')
		ex3  = fierz_expand(ex3, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex3p = load('test_exs/fierz3.p')

		self.assertEqual(ex3, ex3p)

	def test_fierz4(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex4  = Ex(r'\delta_{\alpha}^{\eta} (\Gamma_{a})^{\rho \gamma} - \delta_{\alpha}^{\rho} (\Gamma_{a})^{\eta \gamma}')
		ex4  = fierz_expand(ex4, basis, [r'_{\alpha}', r'_{\gamma}'], [r'_{\eta}', r'_{\rho}'])
		ex4p = load('test_exs/fierz4.p')

		self.assertEqual(ex4, ex4p)

	def test_fierz5(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex5  = Ex(r'''(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}''')
		ex5  = fierz_expand(ex5, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\gamma}', r'_{\eta}'])
		ex5p = load('test_exs/fierz5.p')

		self.assertEqual(ex5, ex5p)

	def test_fierz2index1(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex6  = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} - \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
		ex6  = fierz_expand_2index(ex6, basis, [r'_{\alpha}', r'_{\beta}'])
		ex6p = load('test_exs/fierz2index1.p')

		self.assertEqual(ex6, ex6p)

	def test_fierz2index2(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex7  = Ex(r'\delta_{\alpha}^{\eta} \delta_{\beta}^{\gamma} + \delta_{\beta}^{\eta} \delta_{\alpha}^{\gamma}')
		ex7  = fierz_expand_2index(ex7, basis, [r'_{\alpha}', r'_{\beta}'])
		ex7p = load('test_exs/fierz2index2.p')

		self.assertEqual(ex7, ex7p)

	def test_fierz2index3(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex8  = Ex(r'(\Gamma^{a b})_{\alpha}^{\eta} (\Gamma_{a b})_{\beta}^{\gamma} + (\Gamma^{a b})_{\beta}^{\eta} (\Gamma_{a b})_{\alpha}^{\gamma}')
		ex8  = fierz_expand_2index(ex8, basis, [r'_{\alpha}', r'_{\beta}'])
		ex8p = load('test_exs/fierz2index3.p')

		self.assertEqual(ex8, ex8p)

	def test_fierz2index4(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex9  = Ex(r'\delta_{\alpha}^{\eta} (\Gamma_{a})^{\rho \gamma} - \delta_{\alpha}^{\rho} (\Gamma_{a})^{\eta \gamma}')
		ex9  = fierz_expand_2index(ex9, basis, [r'_{\alpha}', r'_{\gamma}'])
		ex9p = load('test_exs/fierz2index4.p')

		self.assertEqual(ex9, ex9p)

	def test_fierz2index5(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex10  = Ex(r'''(\Gamma^{d e})_{\alpha}^{\gamma} \delta_{\beta}^{\eta} + (\Gamma^{d e})_{\beta}^{\gamma} \delta_{\alpha}^{\eta}''')
		ex10  = fierz_expand_2index(ex10, basis, [r'_{\alpha}', r'_{\beta}'])
		ex10p = load('test_exs/fierz2index5.p')

		self.assertEqual(ex10, ex10p)

	def test_novelfierz1(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex11  = Ex(r'''(\Gamma^{d e})_{\eta \alpha} (\Gamma_{d e a b c})_{\beta}^{\gamma} + (\Gamma^{d e})_{\eta \beta} (\Gamma_{d e a b c})_{\alpha}^{\gamma}''')
		ex11  = fierz_expand(ex11, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex11p = load('test_exs/novelfierz1.p')

		self.assertEqual(ex11, ex11p)

	def test_novelfierz2(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex12  = Ex(r'''(\Gamma_{a b})_{\eta \alpha} (\Gamma^{a b c})_{\beta}^{\gamma} + (\Gamma_{a b})_{\eta \beta} (\Gamma^{a b c})_{\alpha}^{\gamma}''')
		ex12  = fierz_expand(ex12, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex12p = load('test_exs/novelfierz2.p')

		self.assertEqual(ex12, ex12p)

	def test_novelfierz3(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex13  = Ex(r'''(\Gamma_{a})_{\eta \alpha} (\Gamma^{a c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{a c})_{\alpha}^{\gamma}''')
		ex13  = fierz_expand(ex13, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex13p = load('test_exs/novelfierz3.p')

		self.assertEqual(ex13, ex13p)

	def test_novelfierz4(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex14  = Ex(r'''(\Gamma_{a})_{\eta \alpha} (\Gamma^{b c})_{\beta}^{\gamma} + (\Gamma_{a})_{\eta \beta} (\Gamma^{b c})_{\alpha}^{\gamma}''')
		ex14  = fierz_expand(ex14, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex14p = load('test_exs/novelfierz4.p')

		self.assertEqual(ex14, ex14p)

	def test_novelfierz5(self):
		__cdbkernel__ = susy_env()

		basis = [
			Ex(r'C_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d})_{\alpha \beta}'), 
			Ex(r'(\Gamma^{a b c d e})_{\alpha \beta}')
		]

		ex15  = Ex(r'''(\Gamma_{d a})_{\eta \alpha} (\Gamma^{d b c})_{\beta}^{\gamma} + (\Gamma_{d a})_{\eta \beta} (\Gamma^{d b c})_{\alpha}^{\gamma}''')
		ex15  = fierz_expand(ex15, basis, [r'_{\alpha}', r'_{\beta}'], [r'_{\eta}', r'_{\gamma}'])
		ex15p = load('test_exs/novelfierz5.p')

		self.assertEqual(ex15, ex15p)

if __name__ == '__main__':
    unittest.main(verbosity=2)