"""
SusyPy: a symbolic algebra system for supersymmetry calculations.
Copyright (C) 2022  Saul & Isaiah B. Hilsenrath <ihilsenr@umd.edu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sympy
import pickle
import cadabra2 as cdb2

from math import floor, ceil, factorial, cos, sin, pi
from typing import Any, Dict, List, Optional, Tuple, Union
from itertools import chain, permutations, product
from functools import partial, reduce
from collections import OrderedDict

from .node import nth_arg, n_indices
from .solvers import linsolve
from cadabra2 import Ex, ExNode, AntiCommuting, AntiSymmetric, Commuting, Coordinate, Depends, Derivative, EpsilonTensor, GammaMatrix, ImaginaryI, ImplicitIndex, IndexInherit, Indices, Integer, Kernel, KroneckerDelta, PartialDerivative, SelfNonCommuting, Symbol, Symmetric, TableauSymmetry, Trace
from cadabra2 import asym, canonicalise, collect_factors, collect_terms, decompose_product, distribute, eliminate_kronecker, expand_delta, expand_power, factor_in, factor_out, flatten_sum, integrate_by_parts, join_gamma, lower_free_indices, map_sympy, parent_rel_t, product_rule, rename_dummies, sort_product, sort_product, sort_sum, split_gamma, substitute, tree, untrace, unwrap
from cadabra2 import __cdbkernel__, sub, super

D: int = 11
desired_syms: List[int] = [0,0]
syms: List[int] = [0,0]
rep: int = -1
lorentz_indices: List[str] = []
spinor_indices: List[str] = []

def reorder_from_idx(idx: int, a: List[Any]) -> List[Any]:
	return a[idx:] + a[:idx]

def cyclic_perm(a: List[Any]) -> List[Any]:
	return [partial(reorder_from_idx, i) for i in range(len(a))]

def prod(array: List[Ex]) -> Ex:
	if not array:
		return Ex('1')
	else:
		return reduce((lambda x, y: x * y), array)

def ex_sum(array: List[Ex]) -> Ex:
	if not array:
		return Ex('0')
	else:
		return reduce((lambda x, y: x + y), array)

def reorder(x: List[Any], order: List[int]) -> List[Any]:
	return [x[i] for i in order]

def load(file_name: str) -> Any:
	with open(file_name, 'rb') as f:
		value = pickle.load(f)

	return value

def dump(data: Any, file_name: str) -> None:
	with open(file_name, 'wb') as f:
		pickle.dump(data, f)

def do_match(x: Optional[str], y: Optional[str], remove_parity: bool = False) -> bool:
	if x == None or y == None:
		return False
	elif remove_parity:
		return x[2:-1] == y[2:-1]
	else:
		return x == y

def mat_multiply(X: List[List[Ex]], Y: List[List[Ex]]) -> List[List[Ex]]:
	prod = []
	for row in X:
		prod_row = []
		for i in range(len(Y[0])):
			prod_entry = ex_sum([row[j]*Y[j][i] for j in range(len(row))])
			prod_row.append(prod_entry)

		prod.append(prod_row)

	return prod

def parity(perm0: List[Any], perm1: List[Any]) -> int:
	perm1 = perm1[:]
	transCount = 0
	for loc in range(len(perm0) - 1):
		p0 = perm0[loc]
		p1 = perm1[loc]
		if p0 != p1:
			sloc = perm1[loc:].index(p0)+loc
			perm1[loc], perm1[sloc] = p0, p1
			transCount += 1
	if (transCount % 2) == 0:
		return 1
	else:
		return -1

def gen_subs() -> List[Ex]:
	global D, rep, lorentz_indices

	m = D // 2
	sub_exs = []
	indices = lorentz_indices

	for r in range(m+1, 2*m+1):
		exstr = fr'''\Gamma^{{{' '.join(indices[0:r])}}}'''
		res_str_base = fr"({('I '*(m+1))[:-1]})/{factorial(D-r)} \epsilon^{{{' '.join(indices[:D])}}}" + (fr" \Gamma_{{{' '.join(indices[r:D][::-1])}}}" if r < D else '')
		if D % 2 == 0:
			sub_rule = fr"{exstr} -> ({(-1)**r}) {res_str_base} \Gamma'"
		else:
			sub_rule = f'{exstr} -> ({rep}) {res_str_base}'
		sub_exs.insert(0, Ex(sub_rule, False))

	if D % 2 == 0:
		exstr = fr"\Gamma^{{{' '.join(indices[0:m])}}} \Gamma'"
		res_str_base = fr"({(-1)**m}) ({('I '*(m+1))[:-1]})/{factorial(m)} \epsilon^{{{' '.join(indices[:D])}}} \Gamma_{{{' '.join(indices[m:D][::-1])}}}"
		sub_rule = f'{exstr} -> {res_str_base}'
		sub_exs.append(Ex(sub_rule, False))

	return sub_exs

def perform_subs(ex: Ex, subs: List[Ex]) -> Ex:
	for sub in subs:
		substitute(ex, sub)

	return ex

def substitute_to_spinor_metric(ex: Ex) -> Ex:
	global spinor_indices

	alpha = spinor_indices[0]
	beta = spinor_indices[1]

	for rel1, rel2 in product(['_', '^'], repeat=2):
		substitute(ex, Ex(fr'\indexbracket{{1}}{rel1}{{{alpha}}}{rel2}{{{beta}}} -> C{rel1}{{{alpha}}}{rel2}{{{beta}}}', False))

	return ex

def substitute_spinor_zeros(ex: Ex) -> Ex:
	global spinor_indices

	alpha = spinor_indices[0]
	beta = spinor_indices[1]

	for rel1, rel2 in product(['_', '^'], repeat=2):
		substitute(ex, Ex(fr'\indexbracket{{0}}{rel1}{{{alpha}}}{rel2}{{{beta}}} -> 0', False))

	return ex

def swap(ex: Ex, old_indices: List[str], new_indices: List[str]) -> Ex:
	if old_indices and new_indices:
		ex = ex.top().ex()
		for old, new in zip(old_indices, new_indices):
			name = old[2:-1]
			for index in ex[name]:
				if new[0] == '^':
					index.parent_rel = super
				elif new[0] == '_':
					index.parent_rel = sub

		swap_expr = ', '.join([f'{new[0]}{old[1:]} -> {new}' for old, new in zip(old_indices, new_indices)])
		substitute(ex, Ex(swap_expr, False))

	return ex

def swap_dummy(ex: Ex, old_indices: List[str], new_indices: List[str]) -> Ex:
	if old_indices and new_indices:
		swap_expr = ', '.join([f'_{old[1:]} -> _{new[1:]}, ^{old[1:]} -> ^{new[1:]}' for old, new in zip(old_indices, new_indices)])
		substitute(ex, Ex(swap_expr, False))

	return ex

def manipulate_indices(ex: Ex, indices: List[str], conjugate: bool = True) -> List[str]:
		new_indices = []

		for index in indices:
			index = index[2:-1]
			for ex_index in ex[index]:
				if (ex_index.parent_rel == super and conjugate) or (ex_index.parent_rel == sub and not conjugate):
					new_indices.append(f'_{{{ex_index.name}}}')
				elif (ex_index.parent_rel == sub and conjugate) or (ex_index.parent_rel == super and not conjugate):
					new_indices.append(f'^{{{ex_index.name}}}')
				break

		return new_indices

def index_with_rel(index: ExNode, conjugate: bool = False) -> str:
	is_raised = (not conjugate and index.parent_rel == super) or (conjugate and index.parent_rel == sub)
	return ('^' if is_raised else '_') + f'{{{index.name}}}'

def get_free_indices(node: Union[Ex, ExNode]) -> List[Ex]:
	# Find a monomial
	if isinstance(node, Ex):
		node = node.top()
	while node.name == r"\sum" or node.name == r"\equals":
		for child in node.children():
			if child.name != "1":
				node = child
				break
		else:
			return []

	indices = []
	inherits_indices = (IndexInherit.get(node) is not None)

	# Walk through children collecting indices recursively
	for child in node.children():
		if child.parent_rel == parent_rel_t.super or child.parent_rel == parent_rel_t.sub:
			# Discard integers, symbols and coordinates
			symb = Symbol.get(child)
			coord = Coordinate.get(child)
			integer = child.name == "1"
			if not (symb or coord or integer):
				index = child.ex()
				index.top().parent_rel = parent_rel_t.none
				indices.append(index)
		elif inherits_indices:
			indices.extend(get_free_indices(child))

	# Find all repeated indices
	dummies = []
	for pos, idx in enumerate(indices):
		if idx in dummies:
			continue
		try:
			other = indices.index(idx, pos + 1)
			dummies.append(idx)
		except ValueError:
			pass

	# Remove all instances of dummies from the set of indices and return
	for dummy in dummies:
		indices = [i for i in indices if i != dummy]
	return indices

def get_indices(node: Union[Ex, ExNode], rel: bool = False, free = False, dummy: bool = False, spinor: bool = False, lorentz: bool = False, remove_duplicates: bool = False, own_indices: bool = False, conjugate: bool = False) -> List[str]:
	if isinstance(node, Ex):
		node = node.top()
	elif isinstance(node, ExNode):
		pass
	else:
		raise ValueError('Input needs to be of type Ex or ExNode.')

	if free and dummy:
		ValueError('Cannot set both free and dummy to True.')

	if own_indices:
		indices = [(i.name, index_with_rel(i, conjugate=conjugate)) for i in node.own_indices()]

	else:
		if not spinor or lorentz:
			all_indices = [(i.name, index_with_rel(i, conjugate=conjugate)) for i in node.indices()]

		if spinor or lorentz:
			spinors = sum([[(i.name, index_with_rel(i, conjugate=conjugate)) for i in ib.own_indices()] for ib in chain(node[r'\indexbracket'], node['C'], node['D'])], [])

		if free or dummy:
			free_indices = [i.top().name for i in get_free_indices(node)]

		if spinor == lorentz:
			all_indices = [(i, iwr) for i, iwr in all_indices if i in free_indices] if free else all_indices
			all_indices = [(i, iwr) for i, iwr in all_indices if i not in free_indices] if dummy else all_indices
			indices = all_indices

		elif spinor:
			spinors = [(i, iwr) for i, iwr in spinors if i in free_indices] if free else spinors
			spinors = [(i, iwr) for i, iwr in spinors if i not in free_indices] if dummy else spinors
			indices = spinors

		elif lorentz:
			lorentzs = [i for i in all_indices if i not in spinors]
			lorentzs = [(i, iwr) for i, iwr in lorentzs if i in free_indices] if free else lorentzs
			lorentzs = [(i, iwr) for i, iwr in lorentzs if i not in free_indices] if dummy else lorentzs
			indices = lorentzs

	if indices:
		indices = list(list(zip(*indices))[int(rel)])
		indices = [*dict.fromkeys(indices)] if remove_duplicates else indices

	return indices

def index_list_to_str(indices: List[str]) -> str:
	unprimed = ','.join(indices)
	primed = ','.join([f"{i}'" for i in indices])
	double_primed = ','.join([f"{i}''" for i in indices])

	return fr'{{{unprimed},{primed},{double_primed}}}'

def apply_symmetry_to_element(element: str, spinors: List[str], symmetric = True) -> None:
	tableau_param = '{2}' if symmetric else '{1,1}'

	for rel1, rel2 in product(['_', '^'], repeat=2):
		ex = Ex(fr'''{element}{rel1}{{{spinors[0]}}}{rel2}{{{spinors[1]}}}''')
		num_of_indices = len(get_indices(ex))
		TableauSymmetry(ex, Ex(fr'''shape={tableau_param}, indices={{{num_of_indices-2},{num_of_indices-1}}}'''))

def gen_syms(D: int, desired_syms: List[int] = [0,0]) -> List[int]:
	m = D // 2
	syms = [-int(cos(m*pi/2)), -int(sin(m*pi/2))]

	if desired_syms not in [[0,0], [1,1], [1,-1], [-1,1], [-1,-1]]:
		raise ValueError(f'The given desired symmetries are impossible. They can only be [0,0] (null), [1,1], [1,-1], [-1,1], [-1,-1], but your desired_syms = {desired_syms}.')

	if D % 2 == 1 and syms[0] == 0:
		syms[0] = (-1)**(m % 2)*syms[1]
	elif D % 2 == 1 and syms[1] == 0:
		syms[1] = (-1)**(m % 2)*syms[0]
	elif D % 2 == 0 and syms[0] == 0 and syms[1] == desired_syms[1]:
		syms[0] = desired_syms[0]
	elif D % 2 == 0 and syms[1] == 0 and syms[0] == desired_syms[0]:
		syms[1] = desired_syms[1]
	elif D % 2 == 0 and syms[0] == 0 and desired_syms == [0,0]:
		syms[0] = (-1)**int(D % 8 != 0)*syms[1]
	elif D % 2 == 0 and syms[1] == 0 and desired_syms == [0,0]:
		syms[1] = (-1)**int(D % 8 != 0)*syms[0]
	#elif D % 2 == 0 and syms[0] == 0:
	#	raise ValueError(f'The given desired symmetries are impossible. desired_syms[1] must be syms[1] = {syms[1]}.')
	#elif D % 2 == 0 and syms[1] == 0:
	#	raise ValueError(f'The given desired symmetries are impossible. desired_syms[0] must be syms[0] = {syms[0]}.')

	if syms != desired_syms and desired_syms != [0,0]:
	#if D % 2 == 1 and syms != desired_syms and desired_syms != [0,0]:
		raise ValueError(f'The given desired symmetries are impossible. syms = {syms}, while desired_syms = {desired_syms}.')

	return syms

def antisymmetric_gammas(D: int, desired_syms: List[int] = [0,0]) -> List[int]:
	syms = gen_syms(D, desired_syms=desired_syms)
	antisym_gammas = []

	for i in range(1, 2*(D // 2) + 1):
		if (i % 4 == 0 and syms[0] == 1) or (i % 4 == 1 and syms[1] == 1) or (i % 4 == 2 and syms[0] == -1) or (i % 4 == 3 and syms[1] == -1):
			antisym_gammas.append(i)

	return antisym_gammas

def is_gamma_gamma_star_prod_sym(gamma_num_ind: int, D: int, desired_syms: List[int] = [0,0]) -> bool:
	syms = gen_syms(D, desired_syms=desired_syms)
	antisym_gammas = antisymmetric_gammas(D, desired_syms=desired_syms)

	const = 1
	const *= -1 if 2*(D // 2) in antisym_gammas else 1 # if \Gamma' is antisymmetric
	const *= -1 if syms[0] == 1 else 1 # if C is antisymmetric
	const *= -1 if gamma_num_ind in antisym_gammas else 1 # if the other \Gamma is antisymmetric
	const *= -1 if gamma_num_ind % 2 == 1 else 1 # if \Gamma' and the other \Gamma anticommute

	return const == 1

def susy_env(
	lorentz_indices: List[str] = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'],
	spinor_indices: List[str] = [r'\alpha', r'\beta', r'\gamma', r'\eta', r'\rho', r'\sigma'],
	D: int = 11,
	desired_syms: List[int] = [0,0],
	rep: int = -1 
) -> Kernel:
	if True in ["'" in i for i in lorentz_indices] or True in ["'" in i for i in spinor_indices]:
		raise ValueError("lorentz_indices and spinor_indices cannot contain indices with primes (i.e., the character ').")
	elif len(lorentz_indices) < D:
		raise ValueError(f'You need at least {D} Lorentz indices, but you have only {len(lorentz_indices)}.')
	elif len(spinor_indices) < 2:
		raise ValueError(f'You need at least 2 Spinor indices, but you have only {len(spinor_indices)}.')

	all_lorentz_indices = index_list_to_str(lorentz_indices)
	all_spinor_indices = index_list_to_str(spinor_indices)

	Indices(Ex(all_lorentz_indices), Ex(r'vector, position=free'))
	Integer(Ex(all_lorentz_indices), Ex(fr'0..{D}-1'))
	Indices(Ex(all_spinor_indices), Ex(r'''spinor, position=independent'''))

	ImaginaryI(Ex(r'''I'''))
	GammaMatrix(Ex(r'''\Gamma_{#}'''), Ex(r'''metric=\delta''') )
	KroneckerDelta(Ex(r'''\delta_{#}'''))
	EpsilonTensor(Ex(r'''\epsilon_{#}'''), Ex(r'''delta=\delta''') )

	Trace(Ex(r'''Tr{#}'''), Ex(r'''indices=spinor''') )
	Derivative(Ex(r'''D{#}'''))
	PartialDerivative(Ex(r'''\partial{#}'''))
	IndexInherit(Ex(r'''F1{#}'''))
	IndexInherit(Ex(r'''F2{#}'''))
	Symmetric(Ex(fr'''P^{{{' '.join(lorentz_indices[:2])}}}'''))

	syms = gen_syms(D, desired_syms=desired_syms)
	antisym_gammas = antisymmetric_gammas(D, desired_syms=desired_syms)

	apply_symmetry_to_element('C', spinor_indices[:2], syms[0] == -1)
	ImplicitIndex(Ex(r"\Gamma'"), Ex(fr'''\Gamma'_{{{' '.join(spinor_indices[:2])}}}'''))
	apply_symmetry_to_element(r'''\indexbracket{\Gamma'}''', spinor_indices[:2], 2*(D // 2) not in antisym_gammas)

	for i in range(1, 2*(D // 2) + 1):
		ith_gamma = fr'''\Gamma^{{{' '.join(lorentz_indices[:i])}}}'''

		if i % 2 == 0:
			Commuting(Ex(fr'''{{{ith_gamma}, \Gamma'}}'''))
		else:
			AntiCommuting(Ex(fr'''{{{ith_gamma}, \Gamma'}}'''))

		ImplicitIndex(Ex(ith_gamma), Ex(fr'''{ith_gamma}_{{{' '.join(spinor_indices[:2])}}}'''))
		apply_symmetry_to_element(fr'''\indexbracket{{{ith_gamma}}}''', spinor_indices[:2], i not in antisym_gammas)
		apply_symmetry_to_element(fr'''\indexbracket{{\Gamma' {ith_gamma}}}''', spinor_indices[:2], is_gamma_gamma_star_prod_sym(i, D, desired_syms=desired_syms))
		apply_symmetry_to_element(fr'''\indexbracket{{{ith_gamma} \Gamma'}}''', spinor_indices[:2], is_gamma_gamma_star_prod_sym(i, D, desired_syms=desired_syms))

	globals()['D'] = D
	globals()['desired_syms'] = desired_syms
	globals()['syms'] = syms
	globals()['rep'] = rep
	globals()['lorentz_indices'] = lorentz_indices
	globals()['spinor_indices'] = spinor_indices

	return __cdbkernel__

def is_delta_prod(ex: Ex) -> bool:
	return False not in [x.name == r'\delta' for x in ex.top().factors()]

def is_prod_of_gamma_and_gamma_star(ex: Ex) -> bool:
	return ex.top().name == r'\prod' and len(list(ex.top().factors())) == 2 and list(ex[r'\Gamma']) != [] and list(ex[r"\Gamma'"]) != []

def has_sym(ex: Ex) -> bool:
	global D
	gamma_indices = range(2*(D // 2))
	return (ex.top().name == r'\Gamma' and n_indices(ex.top()) in gamma_indices) or is_delta_prod(ex) or ex.top().name == 'C' or ex.top().name == r"\Gamma'" or is_prod_of_gamma_and_gamma_star(ex)

def is_antisym(ex: Ex) -> bool:
	global D, desired_syms
	antisym_gammas = antisymmetric_gammas(D, desired_syms=desired_syms)
	return (ex.top().name == r'\Gamma' and n_indices(ex.top()) in antisym_gammas) or is_delta_prod(ex) or ex.top().name == 'C' or (ex.top().name == r"\Gamma'" and 2*(D // 2) in antisym_gammas) or (is_prod_of_gamma_and_gamma_star(ex) and not is_gamma_gamma_star_prod_sym(n_indices(list(ex[r'\Gamma'])[0]), D, desired_syms=desired_syms))

def spinor_expand(ex: Ex) -> Ex:
	for indexbracket in ex[r'\indexbracket']:
		result = Ex('')
		inner_obj = nth_arg(indexbracket, 0)
		outer_indices = get_indices(indexbracket, rel=True, spinor=True)
		coeff = str(indexbracket.multiplier)
		inner_coeff = str(inner_obj.multiplier)
		inner_obj = Ex(f'1/({inner_coeff})') * inner_obj.ex()

		for term in inner_obj.top().terms():
			term_coeff = str(term.multiplier)
			normalized_term = Ex(f'1/({term_coeff})') * term.ex()
			factors = [x.ex() for x in normalized_term.top().factors() if KroneckerDelta.get(x) or EpsilonTensor.get(x) or ImaginaryI.get(x)]
			spinor_objs = [x.ex() for x in normalized_term.top().factors() if not (KroneckerDelta.get(x) or EpsilonTensor.get(x) or ImaginaryI.get(x))]
			factor_prod = Ex(term_coeff) * Ex(inner_coeff) * Ex(coeff) * prod(factors)
			spinor_obj_prod = prod(spinor_objs)
			expr = factor_prod * Ex(fr'''\indexbracket{{{spinor_obj_prod.input_form()}}}{''.join(outer_indices)}''')
			result = result + expr

		ex[indexbracket] = result

	substitute_to_spinor_metric(ex)
	flatten_prod(ex)

	return ex

def fourier_expand(ex: Ex) -> Ex:
	num = 1
	while len(list(ex[f'F{num}'])):
		for F in ex[f'F{num}']:
			inner_obj = nth_arg(F, 0).ex()
			coeff = str(inner_obj.top().multiplier)
			inner_obj *= Ex(f'1/({coeff})')
			ex[F] = Ex(f'({coeff}) F{num}({inner_obj.input_form()})')

		num += 1

	return ex

def flatten_prod(_ex: Union[Ex, ExNode]) -> Ex:
	if isinstance(_ex, ExNode):
		ex = _ex
	elif isinstance(_ex, Ex):
		ex = _ex.top()
	else:
		raise ValueError('Input needs to be of type Ex or ExNode.')

	if ex.name == r'\prod' or ex.name == '':
		factors = []

		for x in ex.children():
			if x.name == r'\prod' or x.name == '':
				if r'\prod' in [j.name for j in x.children()] or '' in [j.name for j in x.children()]:
					flatten_prod(x)

				coeff = Ex(str(x.multiplier))
				factors += [coeff] + [j.ex() for j in x.children()]
				x.erase()

			elif x.name in [r'\sum', r'\indexbracket', r'\partial', 'D', r'\int', 'Tr'] or (x.name[0] == 'F' and x.name[1:].isnumeric()):
				flatten_prod(x)

		for factor in factors:
			ex.append_child(factor)

	elif ex.name == r'\sum':
		for factor in ex.terms():
			flatten_prod(factor)

	elif ex.name in [r'\indexbracket', r'\partial', 'D', r'\int', 'Tr'] or (ex.name[0] == 'F' and ex.name[1:].isnumeric()):
		flatten_prod(nth_arg(ex, 0))

	return ex

def find_chain(a: List[Tuple[Optional[str], Optional[str], bool]]) -> List[int]:
	a = [list(x) + [i] for i, x in enumerate(a)]
	a.sort(key=lambda x: x[2])
	a.sort(key=lambda x: x[0] == None)
	a.sort(key=lambda x: x[1] == None)

	d = [[[a[0][3]], a[0][:2], a[0][2]]]

	for x in a[1:]:
		connected = -1
		for j, chain in enumerate(d):
			if chain[1][0] != None and ((x[2] and chain[1][0] in x[:2]) or (not x[2] and chain[1][0] == x[1])) and connected == -1:
				chain[0].insert(0, x[3])
				chain[1][0] = x[0] if chain[1][0] != x[0] else x[1]
				chain[2] = chain[2] and x[2]
				connected = j
			elif chain[1][1] != None and ((x[2] and chain[1][1] in x[:2]) or (not x[2] and chain[1][1] == x[0])) and connected == -1:
				chain[0].append(x[3])
				chain[1][1] = x[0] if chain[1][1] != x[0] else x[1]
				chain[2] = chain[2] and x[2]
				connected = j
			elif connected > -1 and do_match(d[connected][1][0], chain[1][0]) and chain[2]:
				d[connected] = [chain[0][::-1] + d[connected][0], [chain[1][1], d[connected][1][1]], chain[2] and d[connected][2]]
				del d[j]
			elif connected > -1 and do_match(d[connected][1][0], chain[1][0]) and d[connected][2]:
				d[connected] = [d[connected][0][::-1] + chain[0], [d[connected][1][1], chain[1][1]], chain[2] and d[connected][2]]
				del d[j]
			elif connected > -1 and do_match(d[connected][1][1], chain[1][1]) and chain[2]:
				d[connected] = [d[connected][0] + chain[0][::-1], [d[connected][1][0], chain[1][0]], chain[2] and d[connected][2]]
				del d[j]
			elif connected > -1 and do_match(d[connected][1][1], chain[1][1]) and d[connected][2]:
				d[connected] = [chain[0] + d[connected][0][::-1], [chain[1][0], d[connected][1][0]], chain[2] and d[connected][2]]
				del d[j]
			elif connected > -1 and do_match(d[connected][1][0], chain[1][1]):
				d[connected] = [chain[0] + d[connected][0], [chain[1][0], d[connected][1][1]], chain[2] and d[connected][2]]
				del d[j]
			elif connected > -1 and do_match(d[connected][1][1], chain[1][0]):
				d[connected] = [d[connected][0] + chain[0], [d[connected][1][0], chain[1][1]], chain[2] and d[connected][2]]
				del d[j]

		if connected == -1:
			d.append([[x[3]], x[:2], x[2]])

	for j, chain in enumerate(d):
		if chain[1][1] != None and chain[2] and (chain[1][0] == None or chain[1][1] < chain[1][0]):
			d[j] = [chain[0][::-1], chain[1][::-1], chain[2]]

	d.sort(key=lambda x: x[1][0] if x[1][0] else 'None')

	return [index for chain in d for index in chain[0]] 

def spinor_combine(ex: Ex) -> Ex:

	for indexbracket in ex[r'\indexbracket']:
		inner_obj = nth_arg(indexbracket, 0)
		coeff = str(inner_obj.multiplier)
		ex[inner_obj] = Ex(f'1/({coeff})')*inner_obj.ex()
		ex[indexbracket] = Ex(coeff)*indexbracket.ex()

	for term in ex.top().terms():
		coeff = Ex(str(term.multiplier))
		non_spinor_factors = [x.ex() for x in term.factors() if len(list(x[r'\indexbracket'])) != 1 and x.name != 'C']
		spinor_factors = [x.ex() for x in term.factors() if len(list(x[r'\indexbracket'])) == 1 or x.name == 'C']
		outer_indices = []

		for n, factor in enumerate(spinor_factors):
			factor_outer_indices = get_indices(factor, rel=True, spinor=True)
			if len(factor_outer_indices) == 1:
				factor_outer_indices.append(None)
				outer_indices.append(factor_outer_indices)
			elif len(factor_outer_indices) == 2:
				outer_indices.append(factor_outer_indices)
			else:
				non_spinor_factors.append(factor)
				del spinor_factors[n]

		num = len(spinor_factors)
		inners = [nth_arg(next(factor[r'\indexbracket']), 0).ex() if list(factor[r'\indexbracket']) else factor for factor in spinor_factors]

		if num >= 2:
			indices_with_syms = [(o[0][2:-1], None, True) if None in o else (o[0][2:-1], o[1][2:-1], has_sym(i)) for o, i in zip(outer_indices, inners)]
			order = find_chain(indices_with_syms)
			inners = reorder(inners, order)
			outer_indices = reorder(outer_indices, order)
			spinor_factors = reorder(spinor_factors, order)

			for i in range(1, num):

				const = 1

				if list(spinor_factors[i-1][r'\indexbracket']) and list(spinor_factors[i][r'\indexbracket']) and (None in outer_indices[i-1] + outer_indices[i] or spinor_factors[i-1].top().name != r'\indexbracket' or spinor_factors[i].top().name != r'\indexbracket'):
					continue

				elif do_match(outer_indices[i-1][1], outer_indices[i][0], remove_parity=True) and outer_indices[i-1][1][0] == '_' and outer_indices[i][0][0] == '^': # inner
					const *= -1 if syms[0] == 1 else 1
					outer_indices[i-1][1] = '^' + outer_indices[i-1][1][1:]
					outer_indices[i][0] = '_' + outer_indices[i][0][1:]

				elif do_match(outer_indices[i-1][0], outer_indices[i][0], remove_parity=True) and outer_indices[i-1][0][0] != outer_indices[i][0][0] and (has_sym(inners[i-1]) or None in outer_indices[i-1]): # front
					outer_indices[i-1] = outer_indices[i-1][::-1]
					const *= -1 if is_antisym(inners[i-1]) else 1
					const *= -1 if outer_indices[i-1][1][0] == '_' and outer_indices[i][0][0] == '^' and syms[0] == 1 else 1
					outer_indices[i-1][1] = '^' + outer_indices[i-1][1][1:]
					outer_indices[i][0] = '_' + outer_indices[i][0][1:]

				elif do_match(outer_indices[i-1][1], outer_indices[i][1], remove_parity=True) and outer_indices[i-1][1][0] != outer_indices[i][1][0] and (has_sym(inners[i]) or None in outer_indices[i]): # back
					outer_indices[i] = outer_indices[i][::-1]
					const *= -1 if is_antisym(inners[i]) else 1
					const *= -1 if outer_indices[i-1][1][0] == '_' and outer_indices[i][0][0] == '^' and syms[0] == 1 else 1
					outer_indices[i-1][1] = '^' + outer_indices[i-1][1][1:]
					outer_indices[i][0] = '_' + outer_indices[i][0][1:]

				elif do_match(outer_indices[i-1][0], outer_indices[i][1], remove_parity=True) and outer_indices[i-1][0][0] != outer_indices[i][1][0] and (has_sym(inners[i-1]) or None in outer_indices[i-1]) and (has_sym(inners[i]) or None in outer_indices[i]): # outer
					outer_indices[i-1] = outer_indices[i-1][::-1]
					outer_indices[i] = outer_indices[i][::-1]
					const *= -1 if is_antisym(inners[i-1]) else 1
					const *= -1 if is_antisym(inners[i]) else 1
					const *= -1 if outer_indices[i-1][1][0] == '_' and outer_indices[i][0][0] == '^' and syms[0] == 1 else 1
					outer_indices[i-1][1] = '^' + outer_indices[i-1][1][1:]
					outer_indices[i][0] = '_' + outer_indices[i][0][1:]

				if do_match(outer_indices[i-1][1], outer_indices[i][0], remove_parity=True) and outer_indices[i-1][1][0] == '^' and outer_indices[i][0][0] == '_':
					coeff1 = str(spinor_factors[i-1].top().multiplier)
					coeff2 = str(spinor_factors[i].top().multiplier)

					if spinor_factors[i-1].top().name == r'\indexbracket' and spinor_factors[i].top().name == r'\indexbracket':
						spinor_factors[i] = Ex(fr'({const}) ({coeff1}) ({coeff2}) \indexbracket{{{(inners[i-1]*inners[i]).input_form()}}}{outer_indices[i-1][0]}{outer_indices[i][1]}')
						inners[i] = inners[i-1]*inners[i]

					elif spinor_factors[i-1].top().name == 'C' and spinor_factors[i].top().name == 'C':
						spinor_factors[i] = Ex(fr'({const}) ({coeff1}) ({coeff2}) C{outer_indices[i-1][0]}{outer_indices[i][1]}')
						inners[i] = Ex(f'C{outer_indices[i-1][0]}{outer_indices[i][1]}')

					elif list(spinor_factors[i-1][r'\indexbracket']) and spinor_factors[i].top().name == 'C':
						indexbracket = list(spinor_factors[i-1][r'\indexbracket'])[0]
						coeff3 = str(indexbracket.multiplier)
						spinor_factors[i-1] = Ex(f'1/({coeff3})') * spinor_factors[i-1]
						spinor_factors[i-1][indexbracket] = Ex(fr'''({const}) ({coeff1}) ({coeff2}) ({coeff3}) \indexbracket{{{inners[i-1].input_form()}}}{outer_indices[i-1][0] or ''}{outer_indices[i][1] or ''}''')
						spinor_factors[i] = spinor_factors[i-1]
						inners[i] = inners[i-1]

					elif spinor_factors[i-1].top().name == 'C' and list(spinor_factors[i][r'\indexbracket']):
						indexbracket = list(spinor_factors[i][r'\indexbracket'])[0]
						coeff3 = str(indexbracket.multiplier)
						spinor_factors[i] = Ex(f'1/({coeff3})') * spinor_factors[i]
						spinor_factors[i][indexbracket] = Ex(fr'''({const}) ({coeff1}) ({coeff2}) ({coeff3}) \indexbracket{{{inners[i].input_form()}}}{outer_indices[i-1][0] or ''}{outer_indices[i][1] or ''}''')


					outer_indices[i][0] = outer_indices[i-1][0]
					inners[i-1] = None
					spinor_factors[i-1] = None
					outer_indices[i-1] = [None, None]

			spinor_factors = [x for x in spinor_factors if x is not None]
			factors = non_spinor_factors + spinor_factors
			ex[term] = coeff*prod(factors)

	return ex

def trace(ex: Ex) -> Ex:
	global D, lorentz_indices

	ex = Ex(f'Tr{{{ex.input_form()}}}')
	last_ex = Ex('1')

	while ex != last_ex:
		last_ex = ex.top().ex()
		join_gamma(ex)
		distribute(ex)
		sort_product(ex)
		perform_subs(ex, gen_subs())
		untrace(ex)
		substitute(ex, Ex(r"\Gamma' \Gamma' -> 1", False))

	canonicalise(ex)
	substitute(ex,  Ex(f'''Tr -> {2**(D // 2)}''', False) )
	substitute(ex, Ex(r'''Tr{\Gamma{#}} -> 0''', False))
	substitute(ex, Ex(r'''Tr{\Gamma'} -> 0''', False))
	substitute(ex, Ex(r'''Tr{\Gamma' \Gamma{#}} -> 0''', False))
	substitute(ex, Ex(r'''Tr{\Gamma{#} \Gamma'} -> 0''', False))

	dummies = [f'_{{{i}}}' for i in get_indices(ex, lorentz=True, dummy=True, remove_duplicates=True)]
	primed_dummies = [i[:-1] + "''}" for i in dummies]
	swap_dummy(ex, dummies, primed_dummies)

	return ex

def evaluate_traces(ex: Ex) -> Ex:
	for spinor_factor in chain(ex[r'\indexbracket'], ex['C']):
		outer_indices = get_indices(spinor_factor, rel=True, spinor=True)
		if (len(outer_indices) == 2) and (outer_indices[0][2:-1] == outer_indices[1][2:-1]):
			coeff = Ex(str(spinor_factor.multiplier))
			inner_obj = nth_arg(spinor_factor, 0).ex() if spinor_factor.name == r'\indexbracket' else Ex('1')
			if outer_indices[0][0] == '^' and outer_indices[1][0] == '_':
				ex[spinor_factor] = trace(Ex('(-1)')*coeff*inner_obj)
			elif outer_indices[0][0] == '_' and outer_indices[1][0] == '^':
				ex[spinor_factor] = trace(coeff*inner_obj)

	return ex

def epsilon_to_delta(ex: Ex, convention: int = -1) -> Ex:
	for p in ex[r'\prod']:
		epsilons = []
		for factor in p.factors():
			if EpsilonTensor.get(factor):
				epsilons.append(factor.ex())
				factor.erase()

		n = len(epsilons)
		epsilon_prod = prod(epsilons)

		for j in range(n // 2):
			old_free_inds = {i[2:-1]: i for i in get_indices(epsilon_prod, rel=True, free=True, remove_duplicates=True)}
			cdb2.epsilon_to_delta(epsilon_prod)
			new_free_inds = {i[2:-1]: i for i in get_indices(epsilon_prod, rel=True, free=True, remove_duplicates=True)}
			thing = [[obj, old_free_inds[key]] for key, obj in new_free_inds.items()]
			free_indices, indices = zip(*thing)
			free_indices, indices = list(free_indices), list(indices)
			epsilon_prod = swap(epsilon_prod, free_indices, indices)

		m = sum([1 for f in epsilon_prod.top().factors() if EpsilonTensor.get(f)])
		power = int((n-m)/2)
		epsilon_prod = Ex(f'({convention**power})')*epsilon_prod
		expand_delta(epsilon_prod)

		p.append_child(epsilon_prod)

def indexbracket_hex(ex: Ex, to_decompose_product=False) -> Ex:
	for indexbracket in ex[r'\indexbracket']:
		inner_obj = nth_arg(indexbracket, 0).ex()
		num = len(list(inner_obj.top().factors()))
		coeff1 = str(indexbracket.multiplier)
		coeff2 = str(inner_obj.top().multiplier)
		inner_obj *= Ex(f'1/({coeff2})')

		coeffs = []
		ibhds = []

		if inner_obj.top().name not in [r'\sum', r'\partial','D', r'Tr', r'\int']:
			for factor in inner_obj.top().factors():
				if factor.name not in [r'\prod', r'\sum', r'\partial','D', r'Tr', r'\int']:
					spinor_ind = get_indices(indexbracket, rel=True, spinor=True)
					lorentz_ind = get_indices(factor, rel=True, lorentz=True)

					obj_str = f'''{factor.name}&{str(len(lorentz_ind))}&{''.join(spinor_ind)}'''.encode('utf-8')
					hexd_obj = f'ibh{obj_str.hex()}'

					coeff = Ex(str(factor.multiplier))
					ibhd = Ex(hexd_obj+''.join(lorentz_ind))


					hexd_ex = Ex(f'{hexd_obj}{{#}}')
					if not SelfNonCommuting.get(hexd_ex):
						SelfNonCommuting(hexd_ex)
					if AntiSymmetric.get(factor) and not AntiSymmetric.get(hexd_ex):
						AntiSymmetric(hexd_ex)
					if Symmetric.get(factor) and not Symmetric.get(hexd_ex):
						Symmetric(hexd_ex)
					if KroneckerDelta.get(factor) and not KroneckerDelta.get(hexd_ex):
						KroneckerDelta(hexd_ex)
					if EpsilonTensor.get(factor) and not EpsilonTensor.get(hexd_ex):
						EpsilonTensor(hexd_ex)
					if ImaginaryI.get(factor) and not ImaginaryI.get(hexd_ex):
						ImaginaryI(hexd_ex)

					coeffs.append(coeff)
					ibhds.append(ibhd)

			if num == len(ibhds):
				ex[indexbracket] = Ex(coeff1) * Ex(coeff2) * prod(coeffs) * prod(ibhds)

	if to_decompose_product:
		decompose_product(ex)

	canonicalise(ex)

	if to_decompose_product:
		collect_terms(ex)

	for p in ex[r'\prod']:
		coeffs = []
		objects = dict()
		for factor in p.factors():
			if factor.name[:3] == 'ibh':
				obj_data = bytes.fromhex(factor.name[3:]).decode('utf-8')
				obj_data = obj_data.split('&')
				name = obj_data[0]
				spinor_ind = obj_data[2]

				lorentz_ind = get_indices(factor, rel=True, lorentz=True)
				lorentz_ind = ''.join(lorentz_ind)

				coeff = Ex(str(factor.multiplier))
				coeffs.append(coeff)

				obj = f'{name}{lorentz_ind}'

				if spinor_ind not in objects.keys():
					objects[spinor_ind] = [obj]
				else:
					objects[spinor_ind].append(obj)

				factor.erase()

		for spinor_ind, obj in objects.items():
			ib = Ex(fr"\indexbracket{{{' '.join(obj)}}}{spinor_ind}")
			p.append_child(ib)

		if coeffs:
			p.append_child(prod(coeffs))

	ex2 = ex.top().ex()

	for node in ex:
		if node.name[:3] == 'ibh':
			obj_data = bytes.fromhex(node.name[3:]).decode('utf-8')
			obj_data = obj_data.split('&')
			name = obj_data[0]
			spinor_ind = obj_data[2]

			lorentz_ind = get_indices(node, rel=True, lorentz=True)
			lorentz_ind = ''.join(lorentz_ind)

			coeff = str(node.multiplier)
			ex2[node] = Ex(coeff)*Ex(fr'\indexbracket{{{name}{lorentz_ind}}}{spinor_ind}')

	ex[ex.top()] = ex2
	sort_product(ex)
	return ex

def evaluate(ex: Ex, to_perform_subs: bool = True, to_join_gamma: bool = True) -> Ex:
	distribute(ex)

	prior_exs = []
	while ex not in prior_exs:
		prior_exs.append(ex.top().ex())

		product_rule(ex)
		spinor_combine(ex)
		evaluate_traces(ex)
		if to_join_gamma:
			join_gamma(ex)
		distribute(ex)
		unwrap(ex)
		sort_product(ex)
		sort_sum(ex)
		epsilon_to_delta(ex)
		#expand_delta(ex)
		eliminate_kronecker(ex)
		rename_dummies(ex)
		indexbracket_hex(ex)
		canonicalise(ex)
		substitute_spinor_zeros(ex)
		substitute(ex, Ex(r"\Gamma' \Gamma' -> 1", False))
		collect_terms(ex)
		spinor_expand(ex)
		fourier_expand(ex)

	if to_perform_subs:
		subs = gen_subs()
		perform_subs(ex, subs)
		spinor_expand(ex)
		canonicalise(ex)
		spinor_combine(ex)
		substitute(ex, Ex(r"\Gamma' \Gamma' -> 1", False))

	collect_factors(ex)

	return ex

def fierz_expand(ex: Ex, basis: List[Ex], commuted: List[str], uncommuted: List[str], to_perform_subs: bool = True) -> Ex:
	global D

	projections = []

	for element in basis:
		evaluate(element, to_perform_subs=False)

		outer_indices = get_indices(element, rel=True, spinor=True)
		conj_uncommuted = manipulate_indices(ex, uncommuted)
		unconj_uncommuted = manipulate_indices(ex, uncommuted, conjugate=False)
		unconj_commuted = manipulate_indices(ex, commuted, conjugate=False)
		lorentz_indices = [f"^{{{index}}}" for index in get_indices(element, lorentz=True, free=True, remove_duplicates=True)]
		primed_indices = [f"{index[:-1]}'}}" for index in lorentz_indices]

		element = swap(element, outer_indices, unconj_commuted)
		element_flipped = swap(element, unconj_commuted, conj_uncommuted[::-1])
		element_primed = swap(element_flipped, lorentz_indices, primed_indices)
		element_unconj = swap(lower_free_indices(element.top().ex()), unconj_commuted, unconj_uncommuted)
		element_unconj_primed = swap(lower_free_indices(element_primed.top().ex()), conj_uncommuted[::-1], unconj_uncommuted)
		element_squared = element_primed*element_unconj

		spinor_combine(element_squared)
		evaluate_traces(element_squared)
		rhs = evaluate(element_squared*element, to_perform_subs=False)
		inv_rhs_const = Ex(f'{2**(D // 2 - 1)}/({rhs.top().multiplier/element.top().multiplier})')

		lhs = evaluate(inv_rhs_const*element_primed*ex, to_perform_subs=False)
		result = evaluate(lhs*element_unconj_primed, to_perform_subs=to_perform_subs)

		projections.append(result)

	ex[ex.top()] = Ex(f'1/({2**(D // 2 - 1)})')*ex_sum(projections)

	return ex

def fierz_expand_2index(ex: Ex, basis: List[Ex], indices: List[str], to_perform_subs: bool = True) -> Ex:
	global D

	projections = []

	for element in basis:
		evaluate(element, to_perform_subs=False)

		outer_indices = get_indices(element, rel=True, spinor=True)
		lorentz_indices = [f"^{{{index}}}" for index in get_indices(element, lorentz=True, free=True, remove_duplicates=True)]
		primed_indices = [f"{index[:-1]}'}}" for index in lorentz_indices]
		conj_indices = manipulate_indices(ex, indices)
		unconj_indices = manipulate_indices(ex, indices, conjugate=False)
		primed_unconj_indices = [f"{index[:-1]}'}}" for index in unconj_indices]

		element = swap(element, outer_indices, unconj_indices)
		element_primed = swap(element, lorentz_indices, primed_indices)
		element_primed_conj = swap(lower_free_indices(element_primed.top().ex()), unconj_indices, conj_indices[::-1])
		element_squared = element*element_primed_conj
		#dummy_field = swap(lower_free_indices(element.top().ex()), unconj_indices, unconj_indices)
		dummy_field = lower_free_indices(element.top().ex())

		spinor_combine(element_squared)
		elsq_trace = evaluate_traces(element_squared)
		rhs = evaluate(elsq_trace*dummy_field, to_perform_subs=False)
		inv_rhs_const = Ex(f'{2**(D // 2 - 1)}/({rhs.top().multiplier/element.top().multiplier})')

		lhs = evaluate(inv_rhs_const*element_primed_conj*ex, to_perform_subs=False)
		lhs = swap_dummy(lhs, unconj_indices, primed_unconj_indices)
		result = evaluate(lhs*element_primed, to_perform_subs=to_perform_subs)

		projections.append(result)

	ex[ex.top()] = Ex(f'1/({2**(D // 2 - 1)})')*ex_sum(projections)

	return ex

def swap_D_and_partial(ex: Ex) -> Ex:
	for D in ex['D']:
		coeff1 = str(D.multiplier)
		D_indices = get_indices(D, rel=True, own_indices=True)
		D_inner_obj = nth_arg(D.ex(), 0).ex()
		if D_inner_obj.top().name == r'\partial':
			coeff2 = str(D_inner_obj.top().multiplier)
			partial_indices = get_indices(D_inner_obj, rel=True, own_indices=True)
			partial_inner_obj = nth_arg(D_inner_obj, 0).ex().input_form()
			ex[D] = Ex(fr'''({coeff1}) ({coeff2}) \partial{''.join(partial_indices)}(D{''.join(D_indices)}({partial_inner_obj}))''')

	return ex

def susy_expand(ex: Ex, susy: str) -> Ex:
	last_ex = Ex('1')
	while ex != last_ex:
		last_ex = ex.top().ex()
		substitute(ex, Ex(susy, False))
		distribute(ex)
		product_rule(ex)
		unwrap(ex)
		swap_D_and_partial(ex)

	return ex

def set_factor_to_value(ex: Ex, consts: List[str], value: str = '') -> str:
	factors = []
	coeff = Ex(str(ex.top().multiplier))

	for factor in ex.top().factors():
		if factor.name == r'\sum' or factor.name in consts or ImaginaryI.get(factor):
			factors.append(factor.ex())
		else:
			factor_coeff = str(factor.multiplier)
			coeff *= Ex(factor_coeff)

	factors.append(coeff)

	if len(factors) > 0:
		factors = prod(factors)
		return f'{factors.input_form()} = {value}' if value else factors.input_form()

def fourier(ex: Ex, fields: List[Ex]) -> Ex:
	field_names = [nth_arg(field, 0).name if field.top().name == r'\indexbracket' else field.top().name for field in fields]
	for term in ex.top().terms():
		num = 1

		for factor in term.children() if term.name == r'\prod' else [term]:

			if factor.name == r'\partial':
				for arg in nth_arg(factor, 0).factors():
					arg_name = nth_arg(arg, 0).name if arg.name == r'\indexbracket' else arg.name

					if arg_name in field_names:
						IndexInherit(Ex(f'F{num}{{#}}'))
						ex[arg] = Ex(f'F{num}({arg.ex().input_form()})')

				ks = []
				for index in factor.own_indices():
					rel = '^' if index.parent_rel == super else '_'
					k = Ex(f'I k{num}{rel}{{{index.name}}}')
					ks.append(k)

				coeff = Ex(str(factor.multiplier))
				new_prod = coeff*prod(ks)*prod([arg.ex() for arg in factor.args()])
				factor.insert(new_prod)

				factor.erase()
				num += 1

	flatten_prod(ex)
	return ex

def inverse_fourier(ex: Ex) -> Ex:
	for term in ex.top().terms():
		num = 1
		while len(list(term[f'k{num}'])):
			coeffs = []
			partial_indices = []
			partial_args = []

			for factor in term.children():

				if factor.name == f'k{num}':
					partial_indices += get_indices(factor.ex(), rel=True)

					coeffs.append(Ex(str(factor.multiplier))*Ex('(-I)'))
					factor.erase()

				elif factor.name == f'F{num}':
					partial_args.append(nth_arg(factor, 0).ex())
					coeffs.append(Ex(str(factor.multiplier)))
					factor.erase()

			term.append_child(prod(coeffs)*Ex(fr'\partial{"".join(partial_indices)}({prod(partial_args).input_form()})'))

			num += 1

	return ex

def find_gauge_param(gauge_trans: Ex) -> Tuple[Ex, Ex]:
	base = [t.ex() for t in gauge_trans.top().terms()][0]
	gauge_param = nth_arg(base, 0).ex()

	return base, gauge_param

def find_non_gauge_inv(ex: Ex, field: Ex, fields: List[Ex], gauge_trans: Ex) -> Ex:
	substitute(ex, Ex(f'{field.input_form()} -> {gauge_trans.input_form()}', False))
	evaluate(ex, to_perform_subs=False)

	if ex == Ex('0'):
		return ex

	_, gauge_param = find_gauge_param(gauge_trans)
	fourier(ex, fields + [gauge_param])

	kf_inds = []
	for term in ex.top().terms():
		t_l_inds = [get_indices(kf, lorentz=True, rel=True) for kf in chain(term[r'k1'], term['F1'])]
		t_l_inds = sum(t_l_inds, [])
		t_l_inds.sort()
		t_l_inds = tuple(t_l_inds)
		t_s_ind = get_indices(next(term['F1']), spinor=True, rel=True)
		t_s_ind = tuple(t_s_ind)
		kf_inds.append((t_l_inds, t_s_ind))

	kf_inds = list(OrderedDict.fromkeys(kf_inds))

	l_inds = get_indices(ex, lorentz=True, rel=True, remove_duplicates=True)
	s_inds = get_indices(ex, spinor=True, rel=True, remove_duplicates=True)

	template_ex = Ex(fr"(-I) \partial_{{a'}}({field.input_form()})")
	indices = get_indices(template_ex, lorentz=True, rel=True, remove_duplicates=True)
	s_index = get_indices(template_ex, spinor=True, rel=True, remove_duplicates=True)

	template_exg = Ex(fr"(-1) \partial_{{a'}}({gauge_trans.input_form()})")
	distribute(template_exg)

	fourier(template_ex, [field])
	fourier(template_exg, [gauge_param])

	sub_ex_rhs = Ex('0')
	zero_ex = Ex('(-1)')*ex
	coeffs = []
	coeffs_ex = []

	i = 0
	for t_l_inds, t_s_ind in kf_inds:
		non_kf_l_inds = [x for x in l_inds if x not in t_l_inds]
		non_kf_s_inds = [x for x in s_inds if x not in t_s_ind]

		for perm in permutations(t_l_inds):
			new_ex = swap(template_ex.top().ex(), indices + s_index, perm + t_s_ind)
			new_exg = swap(template_exg.top().ex(), indices + s_index, perm + t_s_ind)

			ci = Ex(fr"(c{i}{''.join(non_kf_l_inds)}){''.join(non_kf_s_inds)}")
			sub_ex_rhs += ci*new_ex
			zero_ex += ci*new_exg
			coeffs.append(f'c{i}')
			coeffs_ex.append(ci.input_form())
			i += 1

	distribute(zero_ex)
	canonicalise(zero_ex)

	kfs = [kf.ex().input_form() for kf in chain(zero_ex[r'k1'], zero_ex['F1'])]

	factor_out(zero_ex, Ex(', '.join(kfs), False), right=True)

	A = []
	B = []

	for term in zero_ex.top().terms():
		for factor in term.factors():
			factor_name = nth_arg(factor, 0).name if factor.name == r'\indexbracket' else factor.name
			if factor.name == r'\sum' or factor_name in coeffs:
				a = [0] * len(coeffs)
				for f_term in factor.terms():
					name = nth_arg(f_term, 0).name if f_term.name == r'\indexbracket' else f_term.name
					if name in coeffs:
						a[int(name[1:])] = f_term.multiplier
						f_term.erase()

				A.append(a)
				B.append([flatten_sum(Ex('(-1)')*factor.ex())])

	n, m = len(A), len(A[0])

	I = sympy.eye(n)
	A = sympy.Matrix(A)
	AI = sympy.Matrix([[A, I]])
	RP = sympy.Matrix(AI).rref()[0]
	R = RP.extract(range(n), range(m)).tolist()
	P = RP.extract(range(n), range(-n,0))
	P = [[Ex(str(p)) for p in row] for row in P.tolist()]
	PB = mat_multiply(P, B)

	subs0 = ', '.join([f"{coeffs_ex[row.index(1)]} -> {pb[0].input_form()}" for row, pb in zip(R, PB) if 1 in row])
	subs1 = ', '.join([f'{c} -> 0' for c in coeffs_ex])
	substitute(sub_ex_rhs, Ex(subs0, False))
	substitute(sub_ex_rhs, Ex(subs1, False))
	sub_ex_rhs = sub_ex_rhs*Ex('(-I)')
	distribute(sub_ex_rhs)
	inverse_fourier(sub_ex_rhs)
	evaluate(sub_ex_rhs, to_perform_subs=False)

	return sub_ex_rhs

def sunder_gamma(ex: Ex, desired_indices: List[str]) -> Tuple[Ex, Ex]:
	desired_indices.sort()
	lorentz_indices = get_indices(ex, rel=True, lorentz=True)
	desired_indices_pr = [('^' if next(ex[i]).parent_rel == super else '_') + f'{{{i}}}' for i in desired_indices]
	#new_lorentz_indices = sorted(list(set(lorentz_indices) - set(desired_indices_pr))) + desired_indices_pr
	new_lorentz_indices = sorted(list(set(lorentz_indices) - set(desired_indices_pr)), key=lambda x: x[2:-1]) + desired_indices_pr

	if set(lorentz_indices) == set(desired_indices_pr):
		return ex, Ex('0')

	sgn = parity(lorentz_indices, new_lorentz_indices)
	#ex = Ex(str(sgn))*swap(ex, lorentz_indices, new_lorentz_indices)
	ex = swap(ex, lorentz_indices, new_lorentz_indices)

	for i in range(len(desired_indices)):
		split_gamma(ex, on_back=True)
		distribute(ex)

	ex *= Ex(str(sgn))

	sundered_term = next(ex.top().terms()).ex()
	ex -= sundered_term
	collect_terms(ex)
	evaluate(ex, to_perform_subs=False)

	gamma_prod = []
	for gamma in sundered_term[r'\Gamma']:
		gamma_inds = get_indices(gamma.ex())
		if len(gamma_inds) == 1 and gamma_inds[0] in desired_indices:
			gamma_prod.append(gamma.ex())
			gamma.erase()

	gamma_prod = prod(gamma_prod)

	for i in range(len(desired_indices) - 1):
		join_gamma(gamma_prod)
		distribute(gamma_prod)

	sundered_term.top().append_child(gamma_prod)

	ex = sundered_term + ex
	distribute(ex)
	sundered_term = next(ex.top().terms()).ex()
	ex -= sundered_term
	collect_terms(ex)

	return sundered_term, ex

def filter_terms(ex: Ex, field: Ex, indices: List[str], identify_lorentz_proper: bool = True) -> Dict[str, List[Ex]]:
	desired_terms 		 = []
	gauge_terms 		 = []
	undesired_terms 	 = []
	lorentz_proper_terms = []

	indices = [i[2:-1] for i in indices]
	terms 	= [term.ex() for term in ex.top().terms()]

	n = 0
	while n < len(terms):
		term = terms[n]
		n += 1

		if term == Ex('0'):
			continue

		term_free_indices = get_indices(term, free=True)

		partial  = dict()
		field2 	 = dict()
		gamma_ab = dict()
		gamma_ey = dict()

		partial['obj']	 = next(term[r'\partial']).ex()
		partial['index'] = get_indices(partial['obj'], own_indices=True)[0]

		field2['obj']	  = nth_arg(partial['obj'], 0).ex()
		field2['name'] 	  = nth_arg(field2['obj'], 0).name if field2['obj'].top().name == r'\indexbracket' else field2['obj'].top().name
		field2['spinors'] = get_indices(field2['obj'], spinor=True)
		field2['lorentz'] = get_indices(field2['obj'], lorentz=True)
		
		for ib in term[r'\indexbracket']:
			gamma_dict = dict()
			inner_obj  = nth_arg(ib, 0)

			if inner_obj.name == r'\Gamma' or (inner_obj.name == r'\prod' and set([f.name for f in inner_obj.factors()]) == set([r'\Gamma', r"\Gamma'"])):
				gamma_dict['obj'] 	  = ib.ex()
				gamma_dict['name'] 	  = inner_obj.name
				gamma_dict['spinors'] = get_indices(ib, spinor=True)
				gamma_dict['lorentz'] = get_indices(ib, lorentz=True)

				if set(gamma_dict['spinors']) == set(indices): 
					gamma_ab = gamma_dict
				else:
					gamma_ey = gamma_dict

		if not gamma_ey and gamma_ab['name'] == r'\Gamma' and len(gamma_ab['lorentz']) == 1 and gamma_ab['lorentz'][0] == partial['index'] and canonicalise(field2['obj']) == canonicalise(field):
			desired_terms.append(term)

		elif identify_lorentz_proper and gamma_ey and not set(field2['lorentz'] + [partial['index']]) - set(gamma_ey['lorentz']):
			desired_indices = field2['lorentz'] + [partial['index']]
			remainder 		= term.top().ex()
			sundered_term 	= None
			remainder_terms = None

			for ib in term[r'\indexbracket']:
				if set([i.name for i in ib.own_indices()]) == set(gamma_ey['spinors']):
					inner_obj 	= nth_arg(ib, 0)
					gamma 		= next(inner_obj[r'\Gamma'])
					sundered_term, remainder_terms = sunder_gamma(gamma.ex(), desired_indices)
					term[gamma] = sundered_term

			for ib in remainder[r'\indexbracket']:
				if set([i.name for i in ib.own_indices()]) == set(gamma_ey['spinors']):
					inner_obj 		 = nth_arg(ib, 0)
					gamma 			 = next(inner_obj[r'\Gamma'])
					remainder[gamma] = remainder_terms

			evaluate(remainder, to_perform_subs=False)

			lorentz_proper_terms.append(term)
			terms += [t.ex() for t in remainder.top().terms()]

		elif identify_lorentz_proper and gamma_ey and set(field2['lorentz'] + [partial['index']]) & set(gamma_ey['lorentz']):
			index_intersect = set(field2['lorentz'] + [partial['index']]) & set(gamma_ey['lorentz'])
			desired_indices = [sorted(list(index_intersect))[-1]]
			remainder 		= term.top().ex()
			sundered_term 	= None
			remainder_terms = None

			for ib in term[r'\indexbracket']:
				if set([i.name for i in ib.own_indices()]) == set(gamma_ey['spinors']):
					inner_obj 	= nth_arg(ib, 0)
					gamma 		= next(inner_obj[r'\Gamma'])
					sundered_term, remainder_terms = sunder_gamma(gamma.ex(), desired_indices)
					term[gamma] = sundered_term

			for ib in remainder[r'\indexbracket']:
				if set([i.name for i in ib.own_indices()]) == set(gamma_ey['spinors']):
					inner_obj 		 = nth_arg(ib, 0)
					gamma 			 = next(inner_obj[r'\Gamma'])
					remainder[gamma] = remainder_terms

			evaluate(remainder, to_perform_subs=False)

			lorentz_proper_terms.append(term)
			terms += [t.ex() for t in remainder.top().terms()]

		elif partial['index'] in term_free_indices:
			gauge_terms.append(term)

		else:
			undesired_terms.append(term)

	return {
		'desired_terms'		   : desired_terms,
		'gauge_terms'		   : gauge_terms,
		'undesired_terms'	   : undesired_terms,
		'lorentz_proper_terms' : lorentz_proper_terms
	}

def distill_constrs(system: Ex, consts: List[str]) -> Ex:
	system_str = system.input_form()
	prods = [(i, j) for i, j in product(consts, repeat=2) if f'{i} {j}' in system_str]

	forward, inverse, cross_consts = zip(*[(f'{i}*{j} -> {i}{j}', f'{i}{j} -> {i}*{j}', f'{i}{j}') for i, j in prods])

	substitute(system, Ex(', '.join(forward), False))
	sol = linsolve(system, Ex(', '.join(cross_consts), False))

	for arrow in sol[r'\arrow']:
		children = [c.ex() for c in arrow.children()]
		left, right = children[0], children[1]
		if left == right:
			arrow.erase()

	substitute(sol, Ex(', '.join(inverse), False))

	return sol

def susy_solve(bosons: List[Ex], fermions: List[Ex], gauge_transs: List[Ex], susy: str, basis: List[Ex], consts: List[str], indices: List[str], comm_coef: float = 1) -> Tuple[Ex, Dict[Ex, Dict[str, Ex]]]:

	susy_dict = dict()
	system = []
	consts_ex = Ex(', '.join(consts), False)

	for i, field in enumerate(fermions+bosons):
		gauge_trans = gauge_transs[i] if i < len(gauge_transs) else None
		field_name = nth_arg(field, 0).name if field.top().name == r'\indexbracket' else field.top().name

		ex = Ex(rf'D{indices[0]}(D{indices[1]}({field.input_form()})) + D{indices[1]}(D{indices[0]}({field.input_form()}))')
		susy_expand(ex, susy)
		evaluate(ex, to_perform_subs=False)
		fierz_expand_2index(ex, basis, indices, to_perform_subs=False)
		#ex = load(f'../tests/test_exs/save4_{field.input_form()}.p')

		filtered_terms = filter_terms(ex, field, indices, identify_lorentz_proper=list(field[r'\indexbracket']) != [])

		desired_terms 		 = filtered_terms['desired_terms']
		gauge_terms 		 = filtered_terms['gauge_terms']
		undesired_terms 	 = filtered_terms['undesired_terms']
		lorentz_proper_terms = filtered_terms['lorentz_proper_terms']

		pre_desired_terms 	= evaluate(ex_sum(desired_terms))
		pre_gauge_terms		= evaluate(ex_sum(gauge_terms))
		pre_undesired_terms = evaluate(ex_sum(undesired_terms))
		lorentz_proper_exp	= evaluate(ex_sum(lorentz_proper_terms), to_perform_subs=False, to_join_gamma=False)
		
		factor_in(pre_desired_terms, consts_ex)
		factor_in(pre_gauge_terms, consts_ex)
		factor_in(pre_undesired_terms, consts_ex)

		if list(field[r'\indexbracket']) != []:
			added_terms			 = find_non_gauge_inv(lorentz_proper_exp.top().ex(), field, bosons + fermions, gauge_trans) if gauge_trans != Ex('0') else Ex('0')
			filtered_added_terms = filter_terms(added_terms, field, indices, identify_lorentz_proper=False)

			desired_terms 	+= filtered_added_terms['desired_terms']
			gauge_terms 	+= filtered_added_terms['gauge_terms']
			undesired_terms += filtered_added_terms['undesired_terms']

			desired_terms 	= evaluate(ex_sum(desired_terms))
			gauge_terms 	= evaluate(ex_sum(gauge_terms))
			undesired_terms = evaluate(ex_sum(undesired_terms))

			factor_in(desired_terms, consts_ex)
			factor_in(gauge_terms, consts_ex)
			factor_in(undesired_terms, consts_ex)
		else:
			desired_terms, gauge_terms, undesired_terms = pre_desired_terms, pre_gauge_terms, pre_undesired_terms

		lorentz_proper_exp = evaluate(lorentz_proper_exp, to_join_gamma=False)
		factor_in(lorentz_proper_exp, consts_ex)

		system.append(set_factor_to_value(desired_terms, consts, f'({comm_coef}) I'))

		if undesired_terms != Ex('0'):
			for term in undesired_terms.top().terms():
				system.append(set_factor_to_value(term.ex(), consts, '0'))

		susy_dict[field_name] = {
			'desired_terms'		   : pre_desired_terms,
			'gauge_terms'		   : pre_gauge_terms,
			'undesired_terms'	   : pre_undesired_terms,
			'lorentz_proper_terms' : lorentz_proper_exp
		}

	system = Ex(', '.join(system))
	#print('system', system)
	sol = distill_constrs(system, consts)

	return sol, susy_dict

def holoraumy(fields: List[Ex], susy: str, basis: List[Ex], subs: Ex, indices: List[str]) -> Dict[Ex, Dict[str, Ex]]:

	holoraumy_dict = dict()

	for field in fields:
		field_name = nth_arg(field, 0).name if field.top().name == r'\indexbracket' else field.top().name

		ex = Ex(rf'D{indices[0]}(D{indices[1]}({field.input_form()})) - D{indices[1]}(D{indices[0]}({field.input_form()}))')
		susy_expand(ex, susy)
		substitute(ex, subs)
		distribute(ex)
		evaluate(ex, to_perform_subs=False)
		fierz_expand_2index(ex, basis, indices, to_perform_subs=False)

		filtered_terms = filter_terms(ex, field, indices, identify_lorentz_proper=list(field[r'\indexbracket']) != [])

		regular_terms 		 = filtered_terms['desired_terms'] + filtered_terms['undesired_terms']
		gauge_terms 		 = filtered_terms['gauge_terms']
		lorentz_proper_terms = filtered_terms['lorentz_proper_terms']

		regular_terms 		 = evaluate(ex_sum(regular_terms))
		gauge_terms 		 = evaluate(ex_sum(gauge_terms))
		lorentz_proper_terms = evaluate(ex_sum(lorentz_proper_terms), to_join_gamma=False)

		holoraumy_dict[field_name] = {
			'regular_terms'		   : regular_terms,
			'gauge_terms'		   : gauge_terms,
			'lorentz_proper_terms' : lorentz_proper_terms
		}

	return holoraumy_dict

def make_action_susy_inv(lagrangian: Ex, susy: str, consts: List[str], fermion_names: List[str]) -> Ex: 
	susy_expand(lagrangian, susy)
	evaluate(lagrangian, to_perform_subs=False)

	action = Ex(fr'\int{{{lagrangian.input_form()}}}{{obh}}')
	for name in fermion_names:
		integrate_by_parts(action, Ex(fr'\indexbracket{{{name}{{#}}}}{{#}}', False))

	integrand = nth_arg(action, 0).ex()
	evaluate(integrand)
	factor_in(integrand, Ex(', '.join(consts), False))

	system = []
	for term in integrand.top().terms():
		system.append(set_factor_to_value(term.ex(), consts, '0'))

	system = Ex(', '.join(system))
	sol = distill_constrs(system, consts)

	return sol

def substitute_KleinGordon(ex: Ex) -> Ex:
	global lorentz_indices

	index = lorentz_indices[0] + "''"
	substitute(ex, Ex(f'(k1^{{{index}}} k1_{{{index}}})**(-1) -> KleinGordon', False))
	substitute(ex, Ex(f'KleinGordon k1^{{{index}}} k1_{{{index}}} -> 1', False))

	return ex

def proj_coef(r: int, s: float) -> str:
	if r == 0:
		return 1

	quotient = 1
	for i in range(1, r+1):
		quotient *= (2*s - 2*i + 1)

	return f'({(-1)**r * factorial(s)})/({2**r * (factorial(r)*factorial(s - 2*r)*quotient)})'

def sym_proj_decomp(s: float, to_perform_subs: bool = False) -> Ex:
	if int(2*s) % 2 == 0:
		summand = Ex('0')
		for r in range(0, s // 2 + 1):
			summand += Ex(f'({proj_coef(r, s)})' + ' '.join([f'P^{{a{2*i-1} a{2*i}}} P^{{b{2*i-1} b{2*i}}}' for i in range(1, r+1)]) + ' '.join([f'P^{{a{i} b{i}}}' for i in range(2*r+1, s+1)]))
		#map_sympy(summand, "simplify")
		if s > 1:
			asym(summand, Ex(', '.join([f'^{{a{i}}}' for i in range(1, s+1)]), False), antisymmetric=False)
			asym(summand, Ex(', '.join([f'^{{b{i}}}' for i in range(1, s+1)]), False), antisymmetric=False)

		evaluate(summand, to_perform_subs=False)

		return Ex(f"P^{{{' '.join([f'a{i}' for i in range(1, s+1)] + [f'b{i}' for i in range(1, s+1)])}}} -> {summand.input_form()}", False)

	elif int(s*2) % 2 == 1:
		Q = Ex(fr"({ceil(s)})/({int(2*s+2)}) (\Gamma_{{a}} \Gamma_{{b}})_{{\alpha \beta}} P^{{a {' '.join([f'a{i}' for i in range(1, ceil(s))])} b {' '.join([f'b{i}' for i in range(1, ceil(s))])}}}")
		substitute(Q, sym_proj_decomp(ceil(s)))
		evaluate(Q, to_perform_subs=False)
		substitute(Q, Ex(r'P^{a}_{a} -> 3', False))
		substitute(Q, Ex(r'P^{a b} P_{b c} -> P^{a}_{c}', False))
		substitute(Q, Ex(r'P^{a b} P_{c b} -> P^{a}_{c}', False))
		substitute(Q, Ex(r'P^{b a} P_{b c} -> P^{a}_{c}', False))
		substitute(Q, Ex(r'P^{b a} P_{c b} -> P^{a}_{c}', False))
		evaluate(Q, to_perform_subs=to_perform_subs)

		return Ex(fr"(Q^{{{' '.join([f'a{i}' for i in range(1, ceil(s))] + [f'b{i}' for i in range(1, ceil(s))])}}})_{{\alpha \beta}} -> {Q.input_form()}", False)

	else:
		raise ValueError(f'The spin s must be an integer or half-integer, while you inputted s = {s}.')

def sym_prop(s: float, to_perform_subs: bool = False) -> Ex:
	if int(2*s) % 2 == 0:
		D = Ex(f"KleinGordon P^{{{' '.join([f'a{i}' for i in range(1, s+1)] + [f'b{i}' for i in range(1, s+1)])}}}")

	elif int(s*2) % 2 == 1:
		D = Ex(fr"I KleinGordon (\Gamma^{{b}})_{{\alpha}}^{{\gamma}} k1_{{b}} (Q^{{{' '.join([f'a{i}' for i in range(1, ceil(s))] + [f'b{i}' for i in range(1, ceil(s))])}}})_{{\gamma \beta}}")

	else:
		raise ValueError(f'The spin s must be an integer or half-integer, while you inputted s = {s}.')

	#print('\n D', D)
	#print('\n subs', sym_proj_decomp(s))
	substitute(D, sym_proj_decomp(s))
	#print('\n D', D)
	substitute(D, Ex(r'P^{a1 b1} -> \delta^{a1 b1} - (k1^{a1} k1^{b1}) KleinGordon'))
	#print('\n D', D)
	evaluate(D, to_perform_subs=to_perform_subs)
	#print('\n D', D)
	expand_power(D)
	substitute_KleinGordon(D)
	collect_terms(D)

	return D

def rarita_schwinger_prop() -> Ex:
	global D, lorentz_indices

	if D > 2:
		a, b, c = tuple(lorentz_indices[:3])

		return Ex(fr'I KleinGordon \delta_{{{a} {b}}} (\Gamma^{{{c}}})_{{\alpha \beta}} k1_{{{c}}} + (1/{D-2}) I KleinGordon (\Gamma_{{{a}}} \Gamma^{{{c}}} \Gamma_{{{b}}})_{{\alpha \beta}} k1_{{{c}}}')
	else:
		raise ValueError(f'This function only returns Rarita-Schwinger Feynman propagators for dimension D > 2. You have D = {D}.')	

def coupling_current_name(field: Ex) -> str:
	field_name = nth_arg(field, 0).name if field.top().name == r'\indexbracket' else field.top().name
	current_name =  f'''J{field_name.encode('utf-8').hex()}'''
	current = Ex(current_name + '{#}')
	if not Symmetric.get(current) and Symmetric.get(field):
		Symmetric(current)
	elif not AntiSymmetric.get(current) and AntiSymmetric.get(field):
		AntiSymmetric(current)

	return current_name

def gen_prop_sub(field: Ex, propagator: Ex, is_boson: bool = False) -> Ex:
	substitute_KleinGordon(propagator)
	evaluate(propagator, to_perform_subs=False)

	field_spinors = get_indices(field, rel=True, free=True, spinor=True, remove_duplicates=True)
	propagator_spinors = get_indices(propagator, rel=True, free=True, spinor=True, remove_duplicates=True)

	field_lorentzs = get_indices(field, rel=True, free=True, lorentz=True, remove_duplicates=True)
	propagator_lorentzs = get_indices(propagator, rel=True, free=True, lorentz=True, remove_duplicates=True)

	primed_last_propagator_lorentzs = [i[:-1] + "'}" for i in propagator_lorentzs[len(field_lorentzs):]] # new

	for term in propagator.top().terms():
		propagator_lorentzs_wo_rel = get_indices(term, free=True, lorentz=True, remove_duplicates=True)
		k_inds = sum([get_indices(k) for k in term['k1']], [])

		if is_boson and set(k_inds) and not set(k_inds) - set(propagator_lorentzs_wo_rel):
			term.erase()
		elif not is_boson and set(k_inds) & set(propagator_lorentzs_wo_rel):
			term.erase()

	propagator_dummies = [f'_{{{i}}}' for i in get_indices(propagator, dummy=True, remove_duplicates=True)]
	primed_propagator_dummies = [i[:-1] + "''}" for i in propagator_dummies]
	swap_dummy(propagator, propagator_dummies, primed_propagator_dummies)

	if len(field_spinors) == 1:
		#old_prop_inds = propagator_spinors + propagator_lorentzs[:len(field_lorentzs)]
		#new_prop_inds = [field_spinors[0], '^' + propagator_spinors[1][1:]] + field_lorentzs
		old_prop_inds = propagator_spinors + propagator_lorentzs
		new_prop_inds = [field_spinors[0], '^' + propagator_spinors[1][1:]] + field_lorentzs + primed_last_propagator_lorentzs
	else:
		#old_prop_inds = propagator_lorentzs[:len(field_lorentzs)]
		#new_prop_inds = field_lorentzs
		old_prop_inds = propagator_lorentzs
		new_prop_inds = field_lorentzs + primed_last_propagator_lorentzs

	propagator = swap(propagator, old_prop_inds, new_prop_inds)

	coupling_current = coupling_current_name(field)
	#coupling_current_indices = manipulate_indices(propagator, propagator_lorentzs[len(field_lorentzs):])
	coupling_current_indices = manipulate_indices(propagator, new_prop_inds[2+len(field_lorentzs):])

	if len(field_spinors) == 1:
		sub_rule = Ex(rf'''F1({field.input_form()}) -> (-1) ({propagator.input_form()}) F1(({coupling_current}{''.join(coupling_current_indices)})_{{{propagator_spinors[1][1:]}}})''', False)
	else:
		sub_rule = Ex(rf'''F1({field.input_form()}) -> ({propagator.input_form()}) F1(({coupling_current}{''.join(coupling_current_indices)})''', False)

	return sub_rule

def remove_coeffs(ex: Ex) -> List[Ex]:
	new_ex = ex.top().ex()
	term_list = []

	for term in new_ex.top().terms():
		coeff = str(term.multiplier)

		for factor in term.factors():
			if ImaginaryI.get(factor):
				factor.erase()

		term_list.append(Ex(f'1/({coeff})')*term.ex())


	return term_list

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def filter_terms_momentum_space(ex: Ex, field: Ex, propagator_sub: Ex, indices: List[str]) -> Dict[str, List[Ex]]:
	global D

	desired_terms 			   = []
	ignored_terms 			   = []
	gauge_terms 			   = []
	central_charge_terms 	   = []
	undesired_terms 		   = []
	normalized_undesired_terms = []
	propagator_coeff		   = Ex('1')

	indices = [i[2:-1] for i in indices]

	Ixfield = Ex(f'F1({field.input_form()})')
	#print(Ixfield, propagator_sub)
	substitute(Ixfield, propagator_sub)
	distribute(Ixfield)
	eliminate_kronecker(Ixfield)
	#print('\nIxfield', Ixfield)
	Ixfield_dummies = [f'_{{{i}}}' for i in get_indices(Ixfield, spinor=True, dummy=True, remove_duplicates=True)]
	primed_Ixfield_dummies = [i[:-1] + "''}" for i in Ixfield_dummies]
	swap_dummy(Ixfield, Ixfield_dummies, primed_Ixfield_dummies)

	Ixfield_remainder = []

	#print(Ixfield)

	for term in Ixfield.top().terms():
		term_free_indices = get_indices(term, free=True)
		#print(term.ex(), coupling_current_name(field))

		if len(list(term[r'\Gamma'])) == 1 and len(list(term[r'k1'])) == 1 and get_indices(next(term[r'\Gamma']))[0] == get_indices(next(term[r'k1']))[0] and not set(get_indices(next(term[coupling_current_name(field)]), free=True)) - set(term_free_indices):
			propagator_coeff = Ex(set_factor_to_value(term.ex(), []))
			propagator_coeff *= Ex(f'''1/({propagator_coeff.top().multiplier})''')*Ex(f'''1/({propagator_coeff.top().multiplier})''')
			substitute(propagator_coeff, Ex('I -> (-1) I', False))
			#print('\nprop', propagator_coeff)
		else:
			Ixfield_remainder.append(term.ex())

	#print('\nIxfield:', Ixfield)
	#print('propagator_coeff:', propagator_coeff)
	#print('Ixfield_remainder:', Ixfield_remainder, '\n')

	terms = [term.ex() for term in ex.top().terms()]

	n = 0
	while n < len(terms):
		term = terms[n]
		n += 1

		if term == Ex('0'):
			continue

		term_free_indices = get_indices(term, free=True)

		ks 		 = dict()
		current  = dict()
		gamma_ab = dict()
		gamma_ey = dict()

		ks['obj'] 		  = [k.ex() for k in term['k1']]
		ks['indices'] 	  = sum([get_indices(k) for k in ks['obj']], [])
		ks['indiceswrel'] = sum([get_indices(k, rel=True) for k in ks['obj']], [])

		current['obj'] 		   = nth_arg(next(term[r'F1']), 0).ex()
		current['name'] 	   = nth_arg(current['obj'], 0).name if current['obj'].top().name == r'\indexbracket' else current['obj'].top().name
		current['spinors'] 	   = get_indices(current['obj'], spinor=True)
		current['lorentz'] 	   = get_indices(current['obj'], lorentz=True)
		current['lorentzwrel'] = get_indices(current['obj'], rel=True, lorentz=True)

		for ib in term[r'\indexbracket']:
			gamma_dict = dict()
			inner_obj = nth_arg(ib, 0)

			if inner_obj.name == r'\Gamma' or (inner_obj.name == r'\prod' and set([f.name for f in inner_obj.factors()]) == set([r'\Gamma', r"\Gamma'"])):
				gamma_dict['obj'] 		  = ib.ex()
				gamma_dict['name'] 		  = inner_obj.name
				gamma_dict['spinors'] 	  = get_indices(ib, spinor=True)
				gamma_dict['lorentz'] 	  = get_indices(ib, lorentz=True)
				gamma_dict['lorentzwrel'] = get_indices(ib, rel=True, lorentz=True)

				if set(gamma_dict['spinors']) == set(indices): 
					gamma_ab = gamma_dict
				else:
					gamma_ey = gamma_dict

		
		if gamma_ab and gamma_ey and len(gamma_ab['lorentz']) == 1 and len(gamma_ey['lorentz']) == 1 and set(gamma_ab['lorentz'] + gamma_ey['lorentz']) == set(ks['indices']) and not set(current['lorentz']) - set(term_free_indices):
			#print('desired term:', term)

			desired_term = term.top().ex()
			substitute(desired_term, Ex(fr"k1_{{{gamma_ey['lorentz'][0]}}} KleinGordon {gamma_ey['obj'].input_form()} F1({current['obj'].input_form()}) -> 1", False))
			desired_term = propagator_coeff*desired_term*Ex(f'F1({field.input_form()})')
			desired_terms.append(desired_term)
			#print('edited desired term:', desired_term)

			coefficient = term.top().ex()
			substitute(coefficient, Ex(fr"k1_{{{gamma_ab['lorentz'][0]}}} k1_{{{gamma_ey['lorentz'][0]}}} KleinGordon {gamma_ab['obj'].input_form()} {gamma_ey['obj'].input_form()} F1({current['obj'].input_form()}) -> 1", False))
			
			terms += [coefficient*propagator_coeff*Ex(fr"(-1) k1_{{{gamma_ab['lorentz'][0]}}} {gamma_ab['obj'].input_form()}")*r for r in Ixfield_remainder]

		elif set(ks['indices']) & set(current['lorentz']):
			ignored_terms.append(term)
			#print('ignored term:', term)

		elif set(ks['indices']) & set(term_free_indices):
			gauge_terms.append(term)
			#print('gauge term:', term)

		elif not list(term['KleinGordon']):
			central_charge_terms.append(term)
			#print('central charge:', term)

		elif gamma_ey and len(gamma_ey['lorentz']) == D and set(ks['indices']) - set(gamma_ey['lorentz']):
			normalized_undesired_terms += remove_coeffs(term)

			term_coeff 		= str(term.top().multiplier)
			term_inv_coeff 	= Ex(f'1/({term_coeff})')
			term_normalized = term_inv_coeff*term

			indices = gamma_ey['lorentzwrel'] + [[i for i in ks['indiceswrel'] if i[2:-1] not in gamma_ey['lorentz']][0]]
			antisymmetrised_term = asym(term_normalized.top().ex(), Ex(', '.join(indices), False))
			indexbracket_hex(antisymmetrised_term)
			collect_terms(antisymmetrised_term)
			substitute_KleinGordon(antisymmetrised_term)

			original_term = next(antisymmetrised_term.top().terms()).ex()
			remainder = antisymmetrised_term - original_term
			collect_terms(remainder)
			substitute(original_term, Ex(f'{term_normalized.input_form()} -> 1', False))

			if intersection(remove_coeffs(remainder), normalized_undesired_terms):
				undesired_terms.append(term)
				#print('undesired term 2:', term)
			else:
				remainder *= Ex(f'(-1) ({term_coeff})/({original_term.input_form()})')
				terms += [t.ex() for t in remainder.top().terms()]
				#print('Schouten undesired term:', term)
				#print('remainder:', remainder)
				#print('antisymmetrised_term:', antisymmetrised_term)
		else:
			undesired_terms.append(term)
			#print('undesired term:', term)

	return {
		'desired_terms': desired_terms,
		'ignored_terms': ignored_terms,
		'gauge_terms': gauge_terms,
		'central_charge_terms': central_charge_terms,
		'undesired_terms': undesired_terms
	}

def susy_solve_propagator(bosons: List[Ex], fermions: List[Ex], fermion_propagators: List[Ex], susy: str, basis: List[Ex], consts: List[str], indices: List[str], comm_coef: float = 1) -> Ex:

	bosons_sol, _ = susy_solve(bosons, [], [], susy, basis, consts, indices, comm_coef=comm_coef)
	system = bosons_sol[0].input_form()[1:-1].replace('->', '=').split(',')
	consts_ex = Ex(', '.join(consts), False)

	for field in fermions:

		ex = Ex(rf'D{indices[0]}(D{indices[1]}({field.input_form()})) + D{indices[1]}(D{indices[0]}({field.input_form()}))')
		susy_expand(ex, susy)
		evaluate(ex, to_perform_subs=False)
		fierz_expand_2index(ex, basis, indices, to_perform_subs=False)
		#ex = load(f'../tests/test_exs/save4_{field.input_form()}.p')

		fourier(ex, fermions)

		propagator_for_field = None

		for f, p in zip(fermions, fermion_propagators):
			sub_rule = gen_prop_sub(f, p)

			if f == field:
				propagator_for_field = sub_rule

			substitute(ex, sub_rule)

		evaluate(ex, to_perform_subs=False)
		substitute_KleinGordon(ex)

		filtered_terms = filter_terms_momentum_space(ex, field, propagator_for_field, indices)

		desired_terms 		 = filtered_terms['desired_terms']
		undesired_terms 	 = filtered_terms['undesired_terms']

		desired_terms 	= evaluate(ex_sum(desired_terms))
		undesired_terms = evaluate(ex_sum(undesired_terms))

		factor_in(desired_terms, consts_ex)
		factor_in(undesired_terms, consts_ex)

		system.append(set_factor_to_value(desired_terms, consts, f"(-1) ({comm_coef})"))

		if undesired_terms != Ex('0'):
			for term in undesired_terms.top().terms():
				system.append(set_factor_to_value(term.ex(), consts, '0'))

	system = Ex(', '.join(system))
	sol = distill_constrs(system, consts)

	return sol