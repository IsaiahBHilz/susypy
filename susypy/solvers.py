from cadabra2 import *
import sympy

def dsolve(ex, obj=None):
   # FIXME: If ex is a list, convert to an Ex when it is pulled
   # into an expression like below!
   ret = Ex(r'''@(ex)''', False)
   sb = SympyBridge(ret)
   if obj is None:
      sb.from_sympy( str(sympy.dsolve( sb.to_sympy() )) )
   else:
      sbo = SympyBridge(obj)
      sb.from_sympy( str(sympy.dsolve( sb.to_sympy(), sbo.to_sympy() )) )
   return ret

def linsolve(exs, objs):
   ret = Ex(r'''@(exs)''', False)
   obl = Ex(r'''@(objs)''', False)
   if ret.head()!=r"\comma":
      ret = Ex(r'''\comma{ @(ret) }''', False)
   if obl.head()!=r"\comma":
      obl = Ex(r'''\comma{ @(obl) }''', False)
   sb = SympyBridge(ret)
   sb.from_sympy( str(sympy.linsolve( sb.to_sympy(), obl._sympy_())) )
   # The result is a list of solutions. Turn this into a list of
   # rules. Newer sympy versions return a FiniteSet, convert this too.
   if ret.head()==r"FiniteSet":
      if len(ret)==1:
         ret = ret[0]
   if ret.head()!=r"\comma":
      ret = Ex(r''' \comma{ \comma{ @(objs) -> @(ret) } } ''', False)
   else:
      nret = Ex(r'''\comma{}''', False)
      for i in range(len(ret)):
         tmpv=objs[i]
         tmps=ret[i]
         nret.top().append_child(Ex(r''' @(tmpv) -> @(tmps) ''', False))
      ret=Ex(r''' \comma{ @(nret) }''', False)
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
   sb.from_sympy( str(sympy_nonlinsolve( sb.to_sympy(), objs._sympy_())) )
   # The result is a list of solutions. Turn this into a list of
   # rules.
   if len(obl)==1:
      nret = Ex(r'''\comma{}''', False)
      for i in range(len(ret)):
         tmps=ret[i]
         nret.top().append_child(Ex(r''' \comma{ @(objs) -> @(tmps) } ''', False))
      ret=nret
   else:
      nret = Ex(r'''\comma{}''', False)
      for i in range(len(ret)):
         tsol=Ex(r'''\comma{}''', False)
         for v in range(len(objs)):
            tmpv=objs[v]
            tmps=ret[i][v]
            tsol.top().append_child(Ex(r''' @(tmpv) -> @(tmps) ''', False))
         nret.top().append_child(tsol)
      ret=nret
   return ret