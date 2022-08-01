from cadabra2 import *

def n_indices(ex_node):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	return sum(1 for i in ex_node.indices())

def nth_index(ex_node, n):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	for i, index in enumerate(ex_node.indices()):
		if i == n:
			return index
	raise IndexError

def n_args(ex_node):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	return sum(1 for i in ex_node.args())

def nth_arg(ex_node, n):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	for i, arg in enumerate(ex_node.args()):
		if i == n:
			return arg
	raise IndexError

def arg_tuple(ex_node, n):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	return (nth_arg(ex_node, i)._latex_() for i in range(n))

def n_children(ex_node):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()
	return sum(1 for i in ex_node.children())

def nth_child(ex_node, n):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()

	for i, child in enumerate(ex_node.children()):
		if i == n:
			return child
	raise IndexError

def child_tuple(ex_node, n):
	if type(ex_node) is Ex:
		ex_node = ex_node.top()

	return (nth_child(ex_node, i)._latex_() for i in range(n))

def empty_sum():
	ex  = Ex(r''' a + b''')#; _=ex 
	while n_children(ex.top()) != 0:
		nth_child(ex.top(), 0).erase()
	return ex

def empty_product():
	ex  = Ex(r''' a * b''')#; _=ex 
	while n_children(ex.top()) != 0:
		nth_child(ex.top(), 0).erase()
	return ex