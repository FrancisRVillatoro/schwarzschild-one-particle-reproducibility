"""Symbolic checks for paper I, Appendix A and Section 5 (exact rational arithmetic, sympy).
 1. Adiabatic iterates of the friction-free Kasner mode equation at y = 1, with mass mu.
 2. Angular identity (5.6) for the second-order subtraction and the residue C_IR = 1/(240 pi^2).
 3. The function J_mu(lambda) of Eq. (A.5), its massless limit -1/(63 lambda^3) and the residue Delta_IR = -523/(498960 pi^2).
 4. Vanishing of the massless scalar fourth-order average and the sixth-order tails -1/2835, 4/2835.
 5. The polynomials B_2(u) and C_4(u) used in Appendix C.
 6. Coupling shift (delta Box - nabla nabla)(y^-2) on Kasner, the same on Schwarzschild for r^-3, and Kretschmann scalars."""
import sympy as sp
y,u,lam,mu,t=sp.symbols('y u lambda mu t',positive=True)
g=sp.Function('g')(y); om=sp.sqrt(lam**2*g+mu**2); S=1/(4*y**2)
def Q(W): return sp.diff(W,y,2)/W-sp.Rational(3,2)*(sp.diff(W,y)/W)**2
W2=sp.simplify((S-Q(om)/2)/(2*om))
Q1=sp.diff(Q(om+t*W2),t).subs(t,0); W4=(-Q1/2-W2**2)/(2*om)
Q2=sp.diff(Q(om+t*W2+t**2*W4),t,2).subs(t,0)/2; W6=(-Q2/2-2*W2*W4)/(2*om)
gy=u**2*y**sp.Rational(2,3)+(1-u**2)*y**sp.Rational(-4,3)
subs={sp.Derivative(g,(y,k)):sp.diff(gy,y,k).subs(y,1) for k in range(6,0,-1)}
def at1(e): return sp.simplify(e.subs(subs).subs(g,1).subs(y,1))
w2,w4=at1(W2),at1(W4); v=sp.sqrt(lam**2+mu**2)
A0=1/(2*v); A2=-w2/(2*v**2); A4=-w4/(2*v**2)+w2**2/(2*v**3)
I=lambda e,wt=1: sp.simplify(sp.integrate(sp.expand(wt*e),(u,0,1)))
print('(5.6): int S2 du - [1/(2v) - mu^4/(16 v^7)] =',sp.simplify(I(A0+A2)-(1/(2*v)-mu**4/(16*v**7))))
x=sp.symbols('x',positive=True)
print('C_IR = 1/(32 pi^2) int x^2/(1+x^2)^(7/2) =',sp.integrate(x**2/(1+x**2)**sp.Rational(7,2),(x,0,sp.oo))/(32*sp.pi**2),' (1/(240 pi^2) =',1/(240*sp.pi**2),')')
Jmu=I(A0+A2+A4,3*u**2-1)
Jref=-lam**2*(720*lam**8+4680*lam**6*mu**2+9720*lam**4*mu**4+8280*lam**2*mu**6+2520*mu**8)/(45360*v**13)+lam**2*mu**2*(224*lam**4+1498*lam**2*mu**2+12089*mu**4)/(45360*v**13)
print('J_mu - Eq.(A.5) =',sp.simplify(Jmu-Jref),'   J_0 =',sp.simplify(Jmu.subs(mu,0)))
DIR=-sp.integrate(x**6*(224*x**4+1498*x**2+12089)/(1+x**2)**sp.Rational(13,2),(x,0,sp.oo))/(4*sp.pi**2*45360)
print('Delta_IR =',DIR,' (-523/(498960 pi^2) =',-sp.Rational(523,498960)/sp.pi**2,')')
print('massless scalar 4th-order angular average =',sp.simplify(I(A4).subs(mu,0)))
w6=at1(W6.subs(mu,0)); w2m,w4m=w2.subs(mu,0),w4.subs(mu,0)
A6=-w6/(2*lam**2)+w2m*w4m/lam**3-w2m**3/(2*lam**4)
print('6th-order tails: scalar',sp.simplify(I(A6)),'  weighted',sp.simplify(I(A6,3*u**2-1)))
print('B_2(u) =',sp.expand(sp.simplify(A2.subs(mu,0)*lam**3)),'   C_4(u) =',sp.expand(sp.simplify(A4.subs(mu,0)*lam**5)))
# ---- coupling shifts and Kretschmann scalars
def shift(metric,coords,w,yvar):
    gi=metric.inv(); n=len(coords)
    G=[[[sum(gi[a_,d]*(sp.diff(metric[d,b_],coords[c_])+sp.diff(metric[d,c_],coords[b_])-sp.diff(metric[b_,c_],coords[d])) for d in range(n))/2 for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    H=sp.Matrix(n,n,lambda b_,c_: sp.diff(w,coords[b_],coords[c_])-sum(G[a_][b_][c_]*sp.diff(w,coords[a_]) for a_ in range(n)))
    Hm=sp.simplify(gi*H); box=sp.simplify(Hm.trace())
    return [sp.simplify((box-Hm[i,i])) for i in range(n)],G,gi
X,Z1,Z2,th=sp.symbols('X Z1 Z2 theta',positive=True)
gK=sp.diag(-1,y**sp.Rational(-2,3),y**sp.Rational(4,3),y**sp.Rational(4,3)); cK=[y,X,Z1,Z2]
sh,G,gi=shift(gK,cK,y**-2,y); print('Kasner (delta Box - nabla nabla) y^-2  * y^4 =',[sp.simplify(e*y**4) for e in sh])
r,tt,ph=sp.symbols('r t phi',positive=True); M=1
gS=sp.diag(-1/(2*M/r-1),2*M/r-1,r**2,r**2*sp.sin(th)**2); cS=[r,tt,th,ph]
shS,_,_=shift(gS,cS,r**-3,r); print('Schwarzschild (delta Box - nabla nabla) r^-3 * r^6, r->0 =',[sp.limit(sp.simplify(e*r**6),r,0) for e in shS])
def kretschmann(metric,coords):
    gi=metric.inv(); n=len(coords)
    G=[[[sum(gi[a_,d]*(sp.diff(metric[d,b_],coords[c_])+sp.diff(metric[d,c_],coords[b_])-sp.diff(metric[b_,c_],coords[d])) for d in range(n))/2 for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    R=lambda a_,b_,c_,d: sp.diff(G[a_][b_][d],coords[c_])-sp.diff(G[a_][b_][c_],coords[d])+sum(G[a_][c_][e]*G[e][b_][d]-G[a_][d][e]*G[e][b_][c_] for e in range(n))
    Rl=[[[[sum(metric[a_,e]*R(e,b_,c_,d) for e in range(n)) for d in range(n)] for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    return sp.simplify(sum(Rl[a_][b_][c_][d]**2*gi[a_,a_]*gi[b_,b_]*gi[c_,c_]*gi[d,d] for a_ in range(n) for b_ in range(n) for c_ in range(n) for d in range(n)))
print('Kretschmann Kasner * y^4 =',sp.simplify(kretschmann(gK,cK)*y**4),'  anomaly coefficient K/(2880 pi^2) =',sp.simplify(kretschmann(gK,cK)*y**4/(2880*sp.pi**2)),'= 1/(1215 pi^2)')
print('Kretschmann Schwarzschild (M=1) * r^6 =',sp.simplify(kretschmann(gS,cS)*r**6))
