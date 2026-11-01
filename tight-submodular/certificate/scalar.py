"""Independent Arb scalar formulas; no historical proof package is imported."""
from fractions import Fraction as Q
from flint import arb, fmpq, ctx
ctx.prec=256


def ar(x):
    x=Q(x)
    return arb(fmpq(x.numerator,x.denominator))


def endpoint(x, upper=False):
    a=x.upper() if upper else x.lower()
    n,e=a.man_exp()
    return Q(int(n))*Q(2)**int(e)


def decimal(q, digits=32, upper=False):
    q=Q(q); scale=10**digits
    n=q.numerator*scale//q.denominator
    if upper and Q(n,scale)<q: n+=1
    sign='-' if n<0 else ''
    n=abs(n)
    return f'{sign}{n//scale}.{n%scale:0{digits}d}'


class Home:
    def __init__(self):
        self.z=arb(-1).exp(); self.beta=1-self.z
        self.k=1/(4*self.beta); self.ell=arb(2).log()
        self.rs=(1/self.ell).log()/self.ell
        self.eps=1+self.rs-1/self.ell
        self.c=self.z-self.eps
        self.a=ar('.001'); self.ba=(1-(-10*self.a).exp())/10
        self.ea=self.a.exp()
        self.L=ar('1.00032357').log()
        self.halfK=7*self.L/12
        self.bc=43*self.L/75
        self.geometries={}

    def logM(self,D):
        D=ar(D); x=ar(Q(9,16))*D
        return (x*x/(1+x)).log1p().min((ar(Q(5,16))*D).log1p())

    def H(self,s):
        r=(s-self.beta).max(arb(0))
        e=(self.ell*r).exp(); a=1+r-self.eps
        if r<self.rs: return e
        if r>self.rs: return a
        return e.union(a)

    def geometry(self,t,D):
        key=(t,D)
        if key in self.geometries: return self.geometries[key]
        if not D: out=[(arb(0),arb(0))]
        else:
            tb,Db=ar(t),ar(D); qm=self.z*tb
            def clip(v): return v.max(arb(0)).min(qm)
            inv=ar(1/D)
            onset=Db*(self.beta+inv.max(arb(1)).min(1/self.ell).log()/self.ell
                      +(inv-1/self.ell).max(arb(0)))
            level_linear=self.ea-self.c*Db
            level_exp=Db*(self.beta+(self.a-Db.log())/self.ell)
            if Db<self.ell*self.ea: level=level_linear
            elif Db>self.ell*self.ea: level=level_exp
            else: level=level_linear.union(level_exp)
            qs=[arb(0),qm,clip(tb-onset),clip(tb-Db*(self.beta+self.rs)),clip(tb-level)]
            md=self.logM(D);out=[]
            for q in qs:
                s=(Db*self.H((tb-q)/Db)).max(arb(1)).log()
                local=s-((1-(-10*s).exp())/10).min(self.ba)
                out.append((q,md.min(ar(Q(3,5))*s).min(local)))
        self.geometries[key]=out
        return out

    @staticmethod
    def credit(geometry,mu):
        vals=[ar(mu)*q+c for q,c in geometry]
        out=vals[0]
        for v in vals[1:]: out=out.min(v)
        return out
