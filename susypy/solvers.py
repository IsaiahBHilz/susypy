from cadabra2 import *
import sympy

def _iu_to_I(ex):
   # Cadabra2 >= 2.5 maps sympy's imaginary unit to its built-in \iu node;
   # SusyPy uses I (declared ImaginaryI in susy_env), so convert back.
   if len(list(ex[r'\iu'])) > 0:
      substitute(ex, Ex(r'\iu -> I', False))
   return ex

def dsolve(ex, obj=None):
   # FIXME: If ex is a list, convert to an Ex when it is pulled
   # into an expression like below!
   ret = Ex(r'''@(ex)''', False)
   sb = SympyBridge(ret)
   if obj is None:
      sb.from_sympy( sympy.dsolve( sb.to_sympy() ) )
   else:
      sbo = SympyBridge(obj)
      sb.from_sympy( sympy.dsolve( sb.to_sympy(), sbo.to_sympy() ) )
   _iu_to_I(ret)
   return ret

def linsolve(exs, objs):
   ret = Ex(r'''@(exs)''', False)
   obl = Ex(r'''@(objs)''', False)
   if ret.head()!=r"\comma":
      ret = Ex(r'''\comma{ @(ret) }''', False)
   if obl.head()!=r"\comma":
      obl = Ex(r'''\comma{ @(obl) }''', False)
   sb = SympyBridge(ret)
   sb.from_sympy( sympy.linsolve( sb.to_sympy(), obl._sympy_()) )
   # Cadabra2 >= 2.5 converts the sympy result directly (no longer via a
   # string), giving \comma{ \comma{ sol_1, sol_2, ... } }. Turn this into
   # a list of rules (same layout as upstream cdb.sympy.solvers.linsolve).
   nret = Ex(r'''\comma{}''', False)
   for i in range(len(ret[0])):
      tmpv=obl[i]
      tmps=ret[0][i]
      nret.top().append_child(Ex(r''' @(tmpv) -> @(tmps) ''', False))
   ret=Ex(r''' \comma{ @(nret) }''', False)
   _iu_to_I(ret)
   return ret

from sympy.solvers.solveset import nonlinsolve as sympy_nonlinsolve

def nonlinsolve(exs, objs):
   ret = Ex(r'''@(exs)''', False)
   obl = Ex(r'''@(objs)''', False)
   if ret.head()!=r"\comma":
      ret = Ex(r'''\comma{ @(ret) }''', False)
   if obl.head()!=r"\comma":
      obl = Ex(r'''\comma{ @(obl) }''', False)
   sb = SympyBridge(ret)
   sb.from_sympy( sympy_nonlinsolve( sb.to_sympy(), objs._sympy_() ) )
   # The result is a list of solutions. Turn this into a list of
   # rules (same layout as upstream cdb.sympy.solvers.nonlinsolve).
   nret = Ex(r'''\comma{}''', False)
   for i in range(len(ret)):
      tsol = Ex(r'''\comma{}''', False)
      for v in range(len(obl)):
         tmpv=obl[v]
         tmps=ret[i][v]
         tsol.top().append_child(Ex(r''' @(tmpv) -> @(tmps) ''', False))
      nret.top().append_child(tsol)
   ret=nret
   _iu_to_I(ret)
   return ret
