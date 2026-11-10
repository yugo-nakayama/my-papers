#!/usr/bin/env python3
"""Reproduce the seven numerical experiments in the manuscript.

The hard-margin dual is solved directly in its N-dimensional form.  We optimize
``gamma = d * alpha`` against the normalized signed Gram matrix, which avoids
poor scaling when d is large.  All random generators are spawned from one
reported seed.  Use ``--quick`` for a smoke test; omit it for manuscript settings.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
RES = ROOT / "results"


@dataclass
class Theory:
    F: float
    a: float
    b: float
    s1: float = 0.0
    s2: float = 0.0
    c: float = 0.0
    P: float = 0.0
    regime: str = ""


def one_theory(Delta, tau1, tau2, kappa):
    p = Delta / 4.0
    Ds = Delta + tau1 + tau2
    E = p * (tau1 + tau2) + tau1 * tau2
    if kappa >= p:
        den = kappa * Ds + E
        a = 2 * (Delta / 4 + tau2 + kappa) / den
        b = 2 * (kappa - p) / den
        s = 2 * (Delta / 2 + tau2) / den
        regime = "N"
    else:
        a = s = 8 / (Delta + 4 * tau1 + 4 * kappa)
        b = 0.0
        regime = "B2"
    c = 1 - kappa * s
    P = p * (a + b)
    return Theory(a, a, b, 0.0, s, c, P, regime)


def beta(i, x, Delta, tau1, tau2):
    p = Delta / 4
    other_tau = tau2 if i == 1 else tau1
    return p * x / (x + p + other_tau)


def two_theory(Delta, tau1, tau2, k1, k2):
    p = Delta / 4
    Ds = Delta + tau1 + tau2
    E = p * (tau1 + tau2) + tau1 * tau2
    if k2 < beta(2, k1, Delta, tau1, tau2):
        rho1 = 1 / (1 / (p + tau1) + 1 / k1)
        c = (rho1 - k2) / (rho1 + k2)
        a = (1 + c) / (p + tau1)
        b = 0.0
        s1, s2 = (1 + c) / k1, (1 - c) / k2
        F = 2 / (rho1 + k2)
        return Theory(F, a, b, s1, s2, c, p * a, "B2")
    if k1 < beta(1, k2, Delta, tau1, tau2):
        rho2 = 1 / (1 / (p + tau2) + 1 / k2)
        c = -(rho2 - k1) / (rho2 + k1)
        a = 0.0
        b = (1 - c) / (p + tau2)
        s1, s2 = (1 + c) / k1, (1 - c) / k2
        F = 2 / (rho2 + k1)
        return Theory(F, a, b, s1, s2, c, p * b, "B1")
    H = k1 * k2 * Ds + (k1 + k2) * E
    c = (k1 * k2 * (tau1 - tau2) + (k1 - k2) * E) / H
    P = (Delta / 2) * (2 * k1 * k2 + k1 * tau2 + k2 * tau1) / H
    a = (1 + c - P) / tau1
    b = (1 - c - P) / tau2
    s1, s2 = (1 + c) / k1, (1 - c) / k2
    F = 2 * ((k1 + p + tau1) * (k2 + p + tau2) - p**2) / H
    return Theory(F, a, b, s1, s2, c, P, "N")


def _noise(rng, shape, kind="gaussian", scale=1.0):
    if kind == "gaussian":
        return rng.normal(size=shape) * scale
    if kind == "t5":
        return rng.standard_t(5, size=shape) / np.sqrt(5 / 3) * scale
    raise ValueError(kind)


def make_data(rng, d, Delta, n1=10, n2=10, m1=0, m2=0,
              k1=None, k2=None, sigma0=0.1, noise="gaussian"):
    """Generate orthogonal mean contrast and outlier offsets."""
    mu = np.zeros(d)
    mu[: d // 2] = np.sqrt(Delta / (d // 2))
    mu1, mu2 = mu / 2, -mu / 2
    X1 = mu1 + _noise(rng, (n1, d), noise)
    X2 = mu2 + _noise(rng, (n2, d), noise)
    rows, labels, groups = [X1, X2], [-np.ones(n1), np.ones(n2)], ["I1"]*n1 + ["I2"]*n2

    total = m1 + m2
    if total:
        last = np.arange(d // 2, d)
        blocks = np.array_split(last, total)
        out_no = 0
        for cls, m, kap, lab in [(1, m1, k1, -1.0), (2, m2, k2, 1.0)]:
            if m == 0:
                continue
            D = m * kap
            offset2 = D - d * sigma0**2
            if offset2 <= 0:
                raise ValueError(f"target D={D:g} must exceed noise trace {d*sigma0**2:g}")
            out = []
            for _ in range(m):
                nu = np.zeros(d)
                block = blocks[out_no]
                nu[block] = np.sqrt(offset2 / len(block))
                out_no += 1
                out.append(nu + _noise(rng, d, noise, sigma0))
            rows.append(np.asarray(out))
            labels.append(np.full(m, lab))
            groups.extend([f"O{cls}"] * m)
    return np.vstack(rows), np.concatenate(labels), np.asarray(groups), mu1, mu2


def solve_dual(X, t, groups):
    d = X.shape[1]
    Xc = X - X.mean(axis=0, keepdims=True)
    Z = t[:, None] * Xc
    H = (Z @ Z.T) / d
    nneg = np.sum(t < 0)
    npos = np.sum(t > 0)
    x0 = np.where(t < 0, 1 / nneg, 1 / npos)

    def fun(g): return -(g.sum() - 0.5 * g @ H @ g)
    def jac(g): return -(np.ones_like(g) - H @ g)
    cons = {"type": "eq", "fun": lambda g: t @ g, "jac": lambda g: t}
    opt = minimize(fun, x0, jac=jac, method="SLSQP", bounds=[(0, None)]*len(t),
                   constraints=cons, options={"ftol": 1e-11, "maxiter": 1200, "disp": False})
    if not opt.success:
        raise RuntimeError(opt.message)
    g = np.maximum(opt.x, 0)
    alpha = g / d
    w = (alpha * t) @ X
    support = g > max(1e-7, g.max() * 1e-6)
    b = float(np.median(t[support] - X[support] @ w))
    sums = {name: float(g[groups == name].sum()) for name in ["I1", "I2", "O1", "O2"]}
    return {"dL": -opt.fun, "da": sums["I1"], "db": sums["I2"],
            "ds1": sums["O1"], "ds2": sums["O2"], "w": w, "b": b}


def mean_runs(records):
    keys = [k for k in records[0] if k not in {"w"}]
    return {k: float(np.mean([r[k] for r in records])) for k in keys}


def scalar_test_rates(rng, fit, mu1, mu2, ntest=100):
    w, b = fit["w"], fit["b"]
    sd = np.linalg.norm(w)
    y1 = rng.normal(w @ mu1 + b, sd, ntest)
    y2 = rng.normal(w @ mu2 + b, sd, ntest)
    return y1.mean(), y2.mean(), np.mean(y1 >= 0), np.mean(y2 < 0)


def save(df, name):
    df.to_csv(RES / f"{name}.csv", index=False)


def experiment1(rng, quick=False):
    d, Delta, sigma0 = (4000 if quick else 20000), None, 0.2
    Delta = d / 2
    tau = d / 10
    Ds = Delta + 2*tau
    ratios = [0.1, .3, 1, 3, 10, 30]
    reps = 2 if quick else 10
    rows=[]
    for m in [1,3]:
        for r in ratios:
            th=one_theory(Delta,tau,tau,r*Ds)
            rr=[]
            for _ in range(reps):
                X,t,g,_,_=make_data(rng,d,Delta,m2=m,k2=r*Ds,sigma0=sigma0)
                rr.append(solve_dual(X,t,g))
            z=mean_runs(rr)
            rows.append(dict(m=m,kappa_over_Ds=r,dL_sim=z['dL'],dF_th=d*th.F,
                             db_sim=z['db'],db_th=d*th.b,ds_sim=z['ds2'],ds_th=d*th.s2))
    df=pd.DataFrame(rows); save(df,"experiment1")
    x=np.logspace(-1.3,1.5,200); y=[d*one_theory(Delta,tau,tau,q*Ds).F for q in x]
    fig,ax=plt.subplots(figsize=(6.4,4.2),layout="constrained")
    ax.plot(x,y,'k-',label='Theory'); ax.axhline(d*2/Ds,ls='--',c='gray',label='Clean SVM')
    ax.axvline((Delta/4)/Ds,ls=':',c='red',label='Breakdown')
    for m,mk in [(1,'o'),(3,'s')]:
        q=df[df.m==m]; ax.plot(q.kappa_over_Ds,q.dL_sim,mk,mfc='none',label=f'Simulation, m={m}')
    ax.set_xscale('log'); ax.set_xlabel(r'$\kappa/\Delta_*$'); ax.set_ylabel(r'$d\,\max L$'); ax.legend(fontsize=8)
    fig.savefig(FIG/'experiment1.pdf'); plt.close(fig)


def experiment2(rng, quick=False):
    dims=[1000,2000,5000,10000,20000,40000]
    if quick: dims=[500,1000,2000]
    reps=3 if quick else 20
    rows=[]
    for kind in ['gaussian','t5']:
        for d in dims:
            Delta=d/2; tau=d/10; Ds=Delta+2*tau; th=one_theory(Delta,tau,tau,Ds)
            eF=[]; eb=[]
            for _ in range(reps):
                X,t,g,_,_=make_data(rng,d,Delta,m2=1,k2=Ds,sigma0=.2,noise=kind)
                z=solve_dual(X,t,g); eF.append(abs(z['dL']/(d*th.F)-1)); eb.append(abs(z['db']/(d*th.b)-1))
            rows.append(dict(noise=kind,d=d,rel_F=np.mean(eF),rel_b=np.mean(eb),se_F=np.std(eF,ddof=1)/np.sqrt(reps),se_b=np.std(eb,ddof=1)/np.sqrt(reps)))
    df=pd.DataFrame(rows); save(df,"experiment2")
    fig,ax=plt.subplots(figsize=(6.4,4.2),layout="constrained")
    for kind,ls in [('gaussian','-'),('t5','--')]:
        q=df[df.noise==kind]
        ax.loglog(q.d,q.rel_F,ls+'o',label=f'{kind}: dual'); ax.loglog(q.d,q.rel_b,ls+'s',label=f'{kind}: inlier mass')
    ref=np.array(dims,float); ax.loglog(ref,.6/np.sqrt(ref),'k:',label=r'$d^{-1/2}$')
    ax.set_xlabel('$d$'); ax.set_ylabel('mean absolute relative error'); ax.legend(fontsize=8)
    fig.savefig(FIG/'experiment2.pdf'); plt.close(fig)


def experiment3(rng, quick=False):
    d=4000 if quick else 20000; Delta=d/2; tau=d/10; D=2*Delta
    ms=[1,2,4,6,8,10,12,16]; reps=2 if quick else 6; rows=[]
    for m in ms:
        k=D/m; th=one_theory(Delta,tau,tau,k); rr=[]
        for _ in range(reps):
            X,t,g,_,_=make_data(rng,d,Delta,m2=m,k2=k,sigma0=.1); rr.append(solve_dual(X,t,g))
        z=mean_runs(rr); rows.append(dict(m=m,kappa=k,dL_sim=z['dL'],dF_th=d*th.F,db_sim=z['db'],db_th=d*th.b,ds_sim=z['ds2'],ds_th=d*th.s2))
    save(pd.DataFrame(rows),"experiment3")


def experiment4(rng, quick=False, expno=4, two_sided=False):
    d=4000 if quick else 20000; Delta=.15*d; tau=d/10; Ds=Delta+2*tau; sigma0=.05
    ratios=[.03,.06,.1,.15,.19,.25,.4,.7,1.5,5]; reps=2 if quick else (8 if two_sided else 10)
    rows=[]
    for r in ratios:
        k=r*Ds; rr=[]
        for _ in range(reps):
            if two_sided:
                X,t,g,mu1,mu2=make_data(rng,d,Delta,m1=1,m2=1,k1=k,k2=k,sigma0=sigma0)
                th=two_theory(Delta,tau,tau,k,k)
            else:
                X,t,g,mu1,mu2=make_data(rng,d,Delta,m2=1,k2=k,sigma0=sigma0)
                th=one_theory(Delta,tau,tau,k)
            fit=solve_dual(X,t,g); y1,y2,e1,e2=scalar_test_rates(rng,fit,mu1,mu2)
            rr.append(dict(y1=y1,y2=y2,e1=e1,e2=e2))
        z={k0:np.mean([q[k0] for q in rr]) for k0 in rr[0]}
        rows.append(dict(kappa_over_Ds=r,y1_sim=z['y1'],y1_th=-th.P+th.c,y2_sim=z['y2'],y2_th=th.P+th.c,e1=z['e1'],e2=z['e2']))
    df=pd.DataFrame(rows); save(df,f"experiment{expno}")
    return df,(d,Delta,tau,Ds)


def plot_exp4(df, pars):
    d,Delta,tau,Ds=pars; x=np.logspace(-1.6,.8,250); th=[one_theory(Delta,tau,tau,r*Ds) for r in x]
    y1=np.array([-q.P+q.c for q in th]); y2=np.array([q.P+q.c for q in th]); kc=(tau*tau/Delta)/Ds
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    axs[0].semilogx(x,y1,'b-',label=r'Theory $y_1^*$'); axs[0].semilogx(x,y2,'r-',label=r'Theory $y_2^*$')
    axs[0].plot(df.kappa_over_Ds,df.y1_sim,'bo',mfc='none'); axs[0].plot(df.kappa_over_Ds,df.y2_sim,'rs',mfc='none')
    axs[0].axhline(0,c='k',lw=.7); axs[0].axvline(kc,c='green',ls=':'); axs[0].set_xlabel(r'$\kappa/\Delta_*$'); axs[0].set_ylabel('discriminant'); axs[0].legend(fontsize=8)
    axs[1].semilogx(df.kappa_over_Ds,df.e1,'bo-',label='$e(1)$'); axs[1].semilogx(df.kappa_over_Ds,df.e2,'rs-',label='$e(2)$'); axs[1].axvline(kc,c='green',ls=':'); axs[1].set_xlabel(r'$\kappa/\Delta_*$'); axs[1].set_ylabel('misclassification rate'); axs[1].legend(fontsize=8)
    fig.savefig(FIG/'experiment4.pdf'); plt.close(fig)


def experiment5(rng, quick=False):
    one,pars=experiment4(rng,quick,expno=5,two_sided=False); two,_=experiment4(rng,quick,expno="5_symmetric",two_sided=True)
    df=one.merge(two,on='kappa_over_Ds',suffixes=('_one','_two')); save(df,'experiment5_combined')
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    x=df.kappa_over_Ds
    axs[0].semilogx(x,df.y1_sim_one,'o-',label='class 1, one-sided'); axs[0].semilogx(x,df.y1_sim_two,'s-',label='class 1, symmetric'); axs[0].semilogx(x,df.y2_sim_one,'o--',label='class 2, one-sided'); axs[0].semilogx(x,df.y2_sim_two,'s--',label='class 2, symmetric'); axs[0].axhline(0,c='k',lw=.7); axs[0].set_xlabel(r'$\kappa_2/\Delta_*$'); axs[0].set_ylabel('discriminant'); axs[0].legend(fontsize=7)
    axs[1].semilogx(x,df.e1_one,'o-',label='$e(1)$ one-sided'); axs[1].semilogx(x,df.e1_two,'s-',label='$e(1)$ symmetric'); axs[1].semilogx(x,df.e2_one,'o--',label='$e(2)$ one-sided'); axs[1].semilogx(x,df.e2_two,'s--',label='$e(2)$ symmetric'); axs[1].set_xlabel(r'$\kappa_2/\Delta_*$'); axs[1].set_ylabel('misclassification rate'); axs[1].legend(fontsize=7)
    fig.savefig(FIG/'experiment5.pdf'); plt.close(fig)


def experiment6(rng, quick=False):
    d=4000 if quick else 20000; Delta=.15*d; tau=d/10; Ds=Delta+2*tau; k2=.06*Ds; sigma0=.05
    ratios=[.01,.015,.02,.03,.06,.1,.15,.2,.25,.4,1,5]; reps=2 if quick else 8; rows=[]
    for r in ratios:
        k1=r*Ds; th=two_theory(Delta,tau,tau,k1,k2); rr=[]
        for _ in range(reps):
            X,t,g,mu1,mu2=make_data(rng,d,Delta,m1=1,m2=1,k1=k1,k2=k2,sigma0=sigma0)
            fit=solve_dual(X,t,g); y1,y2,e1,e2=scalar_test_rates(rng,fit,mu1,mu2); rr.append((y1,y2,e1,e2))
        z=np.mean(rr,axis=0); rows.append(dict(k1_over_Ds=r,y1_sim=z[0],y1_th=-th.P+th.c,y2_sim=z[1],y2_th=th.P+th.c,c_th=th.c,e1=z[2],e2=z[3],regime=th.regime))
    df=pd.DataFrame(rows); save(df,'experiment6')
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    x=df.k1_over_Ds; axs[0].semilogx(x,df.y1_th,'b-',label='$y_1^*$'); axs[0].semilogx(x,df.y2_th,'r-',label='$y_2^*$'); axs[0].semilogx(x,df.c_th,'k:',label='$c^*$'); axs[0].plot(x,df.y1_sim,'bo',ms=3); axs[0].plot(x,df.y2_sim,'rs',ms=3); axs[0].axhline(0,c='gray',lw=.7); axs[0].set_xlabel(r'$\kappa_1/\Delta_*$'); axs[0].set_ylabel('value'); axs[0].legend(fontsize=8)
    axs[1].semilogx(x,df.e1,'bo-',label='$e(1)$'); axs[1].semilogx(x,df.e2,'rs-',label='$e(2)$'); axs[1].set_xlabel(r'$\kappa_1/\Delta_*$'); axs[1].set_ylabel('misclassification rate'); axs[1].legend(fontsize=8)
    fig.savefig(FIG/'experiment6.pdf'); plt.close(fig)
    return df,(Delta,tau,Ds,k2)


def phase_diagram(exp6,pars):
    Delta,tau,Ds,k2=pars; x=np.logspace(-3,1,400)*Ds
    b2=np.array([beta(2,q,Delta,tau,tau) for q in x]); kc=[]
    for q in x:
        if tau>=Delta/2: kc.append(tau*tau*q/(tau*tau+Delta*q+Delta*tau))
        else: kc.append(max((tau-Delta/4)*q/(q+Delta/4+tau),0))
    kc=np.asarray(kc)
    fig,ax=plt.subplots(figsize=(6.2,4.5),layout="constrained")
    ax.loglog(x/Ds,b2/Ds,'b--',label='class 2 collapse'); ax.loglog(x/Ds,kc/Ds,'b-',label='class 1 lost')
    ax.loglog(b2/Ds,x/Ds,'r--',label='class 1 collapse'); ax.loglog(kc/Ds,x/Ds,'r-',label='class 2 lost')
    colors=np.where(exp6.e1>.5,'blue',np.where(exp6.e2>.5,'red','green'))
    ax.scatter(exp6.k1_over_Ds,np.full(len(exp6),k2/Ds),c=colors,s=18,zorder=5)
    ax.plot([1e-3,10],[1e-3,10],'k:',lw=.8,label='symmetric'); ax.set_xlim(1e-3,10); ax.set_ylim(1e-3,10); ax.set_xlabel(r'$\kappa_1/\Delta_*$'); ax.set_ylabel(r'$\kappa_2/\Delta_*$'); ax.legend(fontsize=7)
    fig.savefig(FIG/'phase_diagram.pdf'); plt.close(fig)


def experiment7(rng, quick=False):
    d=4000 if quick else 20000; Delta=d/2; tau=d/10; Ds=Delta+2*tau; sigma0=.1
    pairs=[(1,1),(.3,.3),(.1,.1),(3,.1),(1,.1),(.3,.1),(.1,.3),(.1,1),(.1,3),(.05,.05)]
    reps=2 if quick else 6; rows=[]
    for r1,r2 in pairs:
        k1,k2=r1*Ds,r2*Ds; th=two_theory(Delta,tau,tau,k1,k2); rr=[]
        for _ in range(reps):
            X,t,g,_,_=make_data(rng,d,Delta,m1=1,m2=1,k1=k1,k2=k2,sigma0=sigma0); rr.append(solve_dual(X,t,g))
        z=mean_runs(rr); rows.append(dict(k1_over_Ds=r1,k2_over_Ds=r2,dL_sim=z['dL'],dF_th=d*th.F,da_sim=z['da'],da_th=d*th.a,db_sim=z['db'],db_th=d*th.b,ds1_sim=z['ds1'],ds1_th=d*th.s1,ds2_sim=z['ds2'],ds2_th=d*th.s2,regime=th.regime))
    save(pd.DataFrame(rows),'experiment7')


def make_latex_tables():
    specs={
      'experiment1':['m','kappa_over_Ds','dL_sim','dF_th','db_sim','db_th','ds_sim','ds_th'],
      'experiment3':['m','kappa','dL_sim','dF_th','db_sim','db_th','ds_sim','ds_th'],
      'experiment4':['kappa_over_Ds','y1_sim','y1_th','y2_sim','y2_th','e1','e2'],
      'experiment5_combined':['kappa_over_Ds','y1_sim_one','y1_th_one','e1_one','e2_one','y1_sim_two','y1_th_two','e1_two','e2_two'],
      'experiment6':['k1_over_Ds','y1_sim','y1_th','y2_sim','y2_th','c_th','e1','e2'],
      'experiment7':['k1_over_Ds','k2_over_Ds','dL_sim','dF_th','da_sim','da_th','db_sim','db_th','ds1_sim','ds1_th','ds2_sim','ds2_th','regime']}
    for name,cols in specs.items():
        df=pd.read_csv(RES/f'{name}.csv')[cols]
        pretty={
          'm':r'$m$', 'kappa_over_Ds':r'$\kappa/\Delta_\ast$', 'kappa':r'$\kappa$',
          'k1_over_Ds':r'$\kappa_1/\Delta_\ast$', 'k2_over_Ds':r'$\kappa_2/\Delta_\ast$',
          'dL_sim':r'$dL$ (sim.)', 'dF_th':r'$dF$ (th.)',
          'da_sim':r'$da$ (sim.)', 'da_th':r'$da^\ast$ (th.)',
          'db_sim':r'$db$ (sim.)', 'db_th':r'$db^\ast$ (th.)',
          'ds_sim':r'$ds$ (sim.)', 'ds_th':r'$ds^\ast$ (th.)',
          'ds1_sim':r'$ds_1$ (sim.)', 'ds1_th':r'$ds_1^\ast$ (th.)',
          'ds2_sim':r'$ds_2$ (sim.)', 'ds2_th':r'$ds_2^\ast$ (th.)',
          'y1_sim':r'$\bar y_1$ (sim.)', 'y1_th':r'$y_1^\ast$ (th.)',
          'y2_sim':r'$\bar y_2$ (sim.)', 'y2_th':r'$y_2^\ast$ (th.)',
          'c_th':r'$c^\ast$', 'e1':r'$e(1)$', 'e2':r'$e(2)$', 'regime':'regime',
          'y1_sim_one':r'$\bar y_1$ one', 'y1_th_one':r'$y_1^\ast$ one',
          'e1_one':r'$e_1$ one', 'e2_one':r'$e_2$ one',
          'y1_sim_two':r'$\bar y_1$ sym.', 'y1_th_two':r'$y_1^\ast$ sym.',
          'e1_two':r'$e_1$ sym.', 'e2_two':r'$e_2$ sym.'}
        df=df.rename(columns=pretty)
        (RES/f'{name}_table.tex').write_text(df.to_latex(index=False,float_format=lambda x:f'{x:.3f}',escape=False),encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--seed',type=int,default=20261029); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    FIG.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)
    master=np.random.SeedSequence(args.seed); rngs=[np.random.default_rng(s) for s in master.spawn(8)]
    experiment1(rngs[0],args.quick); experiment2(rngs[1],args.quick); experiment3(rngs[2],args.quick)
    e4,p4=experiment4(rngs[3],args.quick); plot_exp4(e4,p4)
    experiment5(rngs[4],args.quick); e6,p6=experiment6(rngs[5],args.quick); phase_diagram(e6,p6)
    experiment7(rngs[6],args.quick); make_latex_tables()
    meta={'seed':args.seed,'quick':args.quick,'numpy':np.__version__,'experiments':'1--7'}
    (RES/'run_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta))

if __name__=='__main__':
    main()
