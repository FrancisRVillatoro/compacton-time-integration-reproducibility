#!/usr/bin/env python3
"""
Petrov-Galerkin semidiscretization for |K|(p,p) plus three L-stable integrators:
M2 (Padé R12 fractional-step), TR-BDF2, and two-stage Radau IIA.

Equation:
    U_t + A^{-1}[(mu*h^2 D-c0 B)U + (B+C)g(U)] = 0,
    g(U)=|U|^(p-1)U.
"""
import numpy as np
from scipy.sparse import coo_matrix, bmat
from scipy.sparse.linalg import splu, spsolve

def periodic_matrix(N, offsets, coeffs):
    i=np.arange(N)
    rr=[]; cc=[]; dd=[]
    for off,c in zip(offsets,coeffs):
        rr.append(i); cc.append((i+off)%N); dd.append(np.full(N,c))
    return coo_matrix((np.concatenate(dd),
                       (np.concatenate(rr),np.concatenate(cc))),
                      shape=(N,N)).tocsc()

class AbsKPetrovGalerkin:
    def __init__(self,N,L=100.0,p=2.0,c0=1.0,mu=0.0):
        self.N=N; self.L=L; self.h=L/N
        self.p=p; self.c0=c0; self.mu=mu
        h=self.h
        self.A=periodic_matrix(N,[-2,-1,0,1,2],
                              [1/120,26/120,66/120,26/120,1/120])
        self.B=periodic_matrix(N,[-2,-1,1,2],
                              [-1/(24*h),-10/(24*h),10/(24*h),1/(24*h)])
        self.C=periodic_matrix(N,[-2,-1,1,2],
                              [-1/(2*h**3),2/(2*h**3),-2/(2*h**3),1/(2*h**3)])
        self.D=periodic_matrix(N,[-2,-1,0,1,2],
                              [1/h**4,-4/h**4,6/h**4,-4/h**4,1/h**4])
        self.BC=(self.B+self.C).tocsc()
        self.Lop=(mu*h*h*self.D-c0*self.B).tocsc()
        self.Alu=splu(self.A)

    def g(self,U):
        return np.abs(U)**(self.p-1)*U

    def gp(self,U):
        return self.p*np.abs(U)**(self.p-1)

    def G(self,U):
        return self.Lop@U+self.BC@self.g(U)

    def F(self,U):
        return self.Alu.solve(self.G(U))

    def JG(self,U):
        return (self.Lop+self.BC.multiply(self.gp(U))).tocsc()

def damped_newton(residual_jacobian,x0,tol=1e-10,maxit=35):
    x=x0.copy()
    for it in range(maxit+1):
        R,J=residual_jacobian(x)
        nrm=float(np.linalg.norm(R,np.inf))
        if nrm<tol:
            return x,it,True,nrm
        if it==maxit:
            return x,it,False,nrm
        dx=spsolve(J,-R)
        if not np.all(np.isfinite(dx)):
            return x,it,False,nrm
        best_nrm=np.inf; best_x=None
        for s in (1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125):
            xt=x+s*dx
            if not np.all(np.isfinite(xt)):
                continue
            rt,_=residual_jacobian(xt)
            nt=float(np.linalg.norm(rt,np.inf))
            if nt<best_nrm:
                best_nrm=nt; best_x=xt
            if nt<nrm:
                break
        if best_x is None:
            return x,it,False,nrm
        x=best_x
    return x,maxit,False,np.inf

def step_M2(model,U,dt,Uprev=None):
    FU=model.F(U)
    W=U-dt*FU/3.0
    N=model.N
    X=U-dt*FU if Uprev is None else U+(U-Uprev)
    V=X+dt*model.F(X)/4.0
    Y0=np.r_[X,V]
    def RJ(Y):
        X=Y[:N]; V=Y[N:]
        R1=model.A@((V-X)/dt)-0.25*model.G(X)
        R2=model.A@((X-W)/dt)+(2.0/3.0)*model.G(V)
        J=bmat([
            [-model.A/dt-0.25*model.JG(X), model.A/dt],
            [ model.A/dt,                  (2.0/3.0)*model.JG(V)]
        ],format="csc")
        return np.r_[R1,R2],J
    Y,it,ok,nrm=damped_newton(RJ,Y0)
    return Y[:N],it,ok,nrm

def step_TRBDF2(model,U,dt,Uprev=None):
    gamma=2.0-np.sqrt(2.0)
    # Stage 1: trapezoidal rule over gamma*dt.
    guess1=U-gamma*dt*model.F(U)
    def RJ1(Y):
        R=model.A@((Y-U)/(gamma*dt))+0.5*(model.G(U)+model.G(Y))
        J=model.A/(gamma*dt)+0.5*model.JG(Y)
        return R,J.tocsc()
    Y,it1,ok1,nrm1=damped_newton(RJ1,guess1)
    if not ok1:
        return Y,it1,False,nrm1
    # Stage 2: variable-step BDF2 continuation.
    a=1.0/(gamma*(2.0-gamma))
    b=-(1.0-gamma)**2/(gamma*(2.0-gamma))
    d=(1.0-gamma)/(2.0-gamma)
    guess2=Y-(1.0-gamma)*dt*model.F(Y)
    def RJ2(X):
        R=model.A@((X-a*Y-b*U)/dt)+d*model.G(X)
        J=model.A/dt+d*model.JG(X)
        return R,J.tocsc()
    X,it2,ok2,nrm2=damped_newton(RJ2,guess2)
    return X,it1+it2,(ok1 and ok2),max(nrm1,nrm2)

def step_RadauIIA2(model,U,dt,Uprev=None):
    # c=(1/3,1), A=[[5/12,-1/12],[3/4,1/4]], stiffly accurate.
    N=model.N
    FU=model.F(U)
    Y1=U-(dt/3.0)*FU
    Y2=U-dt*FU if Uprev is None else U+(U-Uprev)
    Z0=np.r_[Y1,Y2]
    def RJ(Z):
        Y1=Z[:N]; Y2=Z[N:]
        G1=model.G(Y1); G2=model.G(Y2)
        R1=model.A@((Y1-U)/dt)+(5.0/12.0)*G1-(1.0/12.0)*G2
        R2=model.A@((Y2-U)/dt)+(3.0/4.0)*G1+(1.0/4.0)*G2
        J=bmat([
            [model.A/dt+(5.0/12.0)*model.JG(Y1), -(1.0/12.0)*model.JG(Y2)],
            [(3.0/4.0)*model.JG(Y1), model.A/dt+(1.0/4.0)*model.JG(Y2)]
        ],format="csc")
        return np.r_[R1,R2],J
    Z,it,ok,nrm=damped_newton(RJ,Z0)
    return Z[N:],it,ok,nrm

def scalar_validation():
    # Integrate y'=-y from 0 to 1 with a scalar stand-in solved analytically
    # through the known stability functions; this validates the expected orders.
    def Rtr(z):
        g=2-np.sqrt(2)
        # Apply the two scalar stages exactly.
        y=1.0
        Y=y*(1+g*z/2)/(1-g*z/2)
        a=1/(g*(2-g)); b=-(1-g)**2/(g*(2-g)); d=(1-g)/(2-g)
        return (a*Y+b*y)/(1-d*z)
    def Rrad(z):
        return (1+z/3)/(1-2*z/3+z*z/6)
    def Rm2(z):
        return (1+z/3)/(1-2*z/3+z*z/6)
    exact=np.exp(-1.0)
    rows=[]
    for name,R in [("TR-BDF2",Rtr),("RadauIIA2",Rrad),("M2",Rm2)]:
        prev=None
        for dt in [0.2,0.1,0.05,0.025]:
            n=int(round(1/dt)); val=R(-dt)**n
            e=abs(val-exact)
            p=np.nan if prev is None else np.log(prev/e)/np.log(2)
            rows.append((name,dt,e,p))
            prev=e
    return rows
