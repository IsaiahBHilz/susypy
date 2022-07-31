#!/usr/local/bin/cadabra2
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
#__cdbtmp__ = Coordinate(Ex(r'''x'''), Ex(r'')); display(__cdbtmp__)
#__cdbtmp__ = Depends(Ex(r'''{f,g}'''), Ex(r'''x''') ); display(__cdbtmp__)
#__cdbtmp__ = PartialDerivative(Ex(r'''\partial{#}'''), Ex(r'')); display(__cdbtmp__)
#ex = Ex(r''' \partial_{x x}{f} + \partial_{x}{f} = \sin(x)'''); _=ex; display(ex)
#_ = sol=dsolve(ex, Ex(r'''f''', False)); display(_)
#_ = substitute(ex, sol); display(_)
#_ = map_sympy(ex, "simplify"); display(_)
#st = Ex(r''' { \partial_{x}{f} = f g \sin(x), \partial_{x}{g} = g**2 \sin(x) }'''); _=st; display(st)
#_ = dsolve(st); display(_)
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
#ex = Ex(r'''  x + a = 0'''); _=ex; display(ex)
#_ = sol=linsolve(ex, Ex(r'''x''', False)); display(_)
#_ = substitute(ex, sol[0]); display(_)
#ex = Ex(r''' x + 1 = 0, y + 4/3 x + 2 a = 0'''); _=ex; display(ex)
#_ = sol=linsolve(ex, Ex(r'''x,y''', False)); display(_)
#_ = substitute(ex, sol[0]); display(_)
#ex = Ex(r''' x A_{m} C^{m} + b D_{n p} G^{n p} = 0'''); _=ex; display(ex)
#_ = linsolve(ex, Ex(r'''x''', False)); display(_)
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
#ex = Ex(r''' x**2-1'''); _=ex; display(ex)
#_ = sol=nonlinsolve(ex, Ex(r'''x''', False)); display(_)
#substitute(ex, sol[0])
#_ = simplify(ex); display(_)
#ex = Ex(r''' x**2 + A_{m} B^{m} = 0'''); _=ex; display(ex)
#_ = sol=nonlinsolve(ex, Ex(r'''x''', False)); display(_)
#ex = Ex(r''' { x**2 - 1 = 0, y**2 - 2 =0 }'''); _=ex; display(ex) 
#_ = nonlinsolve(ex, Ex(r'''x,y''', False)); display(_)
