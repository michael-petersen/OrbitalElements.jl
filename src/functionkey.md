## Where are all the functions?
And more importantly, what is the full list of arguments?

### Frequencies.jl:
- αβFromFrequencies(Ω1,Ω2,Ω₀)
- FrequenciesFromαβ(α,β,Ω₀)
- FrequenciesDerivsFromαβDerivs(α,β,dα,dβ,Ω₀)
- αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
- ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
- ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,TOLA)
- ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,VERBOSE,NINT,EDGE,Ω₀)
- ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
- ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)
- ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e,da,de,TOLA,TOLECC,NINT)
- ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,e,eps,maxiter,da,de,TOLA,TOLECC)
- JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,NINT,Ω₀)
- JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT,EDGE,TOLA)

### Circular/CircularFrequencies.jl
- Ω1circular(dψ,d2ψ,a)
- Ω1circular(dψ,d2ψ,d3ψ,d4ψ,a,e)
- dΩ1circular(dψ,d2ψ,d3ψ,a)
- αcircular(dψ,d2ψ,d3ψ,d4ψ,a,e,Ω₀)
- Ω2circular(dψ,r)
- Ω2circular(dψ,d2ψ,r)
- dΩ2circular(dψ,d2ψ,d3ψ,a)
- βcircular2ndorderExpansionCoefs(ψ,dψ,d2ψ,d3ψ,d4ψ,a)
- βcircular(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
- βcirc(αcirc,dψ,d2ψ,Ω₀,rmin,rmax)
- RcircFromΩ1circ(Ω₁,dψ,d2ψ,rmin,rmax,tolx,tolf)
- RcircFromΩ2circ(Ω₂,dψ,d2ψ,rmin,rmax)

### Henon/Frequencies.jl:
If one wanted to write a new frequency anomaly, this is the template to use.
- HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLECC)
- αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
- DαβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
- DFrequenciesHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

### Henon/Ufunc.jl:
- henon_f(u)
- henon_df(u)
- henon_d2f(u)
- henon_d3f(u)
- henon_d4f(u)
- ru(u,a,e)
- drdu(u,a,e)
- ψeff(ψ,r,L)
- dψeffdr(dψ,r,L)
- d2ψeffdr2(d2ψ,r,L)
- ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)
- ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC)
- dΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,da,de,TOLA,TOLECC,EDGE)

### Resonance/ABtoUV.jl:
- αβFromUV(u,v,n1,n2,ωmin,ωmax)
- UVFromαβ(α,β,n1,n2,ωmin,ωmax)
- JacαβToUV(n1,n2,ωmin,ωmax,v)

### Resonance/UVbounds.jl:
- Getϖ(ω₀,n₁,n₂,dψ,d2ψ,Ω₀,rmin,rmax)
- Getϖ(ω,ωmin,ωmax)
- Findωminωmax(n₁,n₂,dψ,d2ψ,αmin,αmax,Ω₀,rmin,rmax)
- FindVminVmax(u,n₁,n₂,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀,rmin,rmax)
- HUFunc(u,ωmin,ωmax)
- FindVbound(n₁,n₂,dψ,d2ψ,Ω₀,rmin,rmax)

### Utils/ComputeEL.jl:
- EccentricityTolerance(a,TOLA,TOLECC)
- EFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
- LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
- ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
- Erad(ψ,a)
- Ecirc(ψ,dψ,a,e)
- EcircExpansion(ψ,dψ,d2ψ,d3ψ,a,e)
- Lcirc(dψ,a)
- Lcirc2ndorderExpansionCoefs(dψ,d3ψ,a)
- LcircExpansion(dψ,d3ψ,a,e)
- dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)
- dELcircExpansion(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
- JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)
- EFromRpRa(ψ,dψ,d2ψ,d3ψ,rp,ra,TOLA,TOLECC)
- LFromRpRa(ψ,dψ,d2ψ,d3ψ,rp,ra,TOLA,TOLECC)
- ELFromRpRa(ψ,dψ,d2ψ,d3ψ,rp,ra,TOLA,TOLECC)
- RcircFromL(L,dψ,rmin,rmax,tolx,tolf)

### Utils/NumericalInversion.jl:
- inverse2Dlinear(a,b,c,d,y1,y2)
- nextguess(acur,ecur,adir,edir)
- AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε,Ω₀)
- AEFromJLBrute(J,L,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)

### Utils/OrbitDefinitions.jl:
- AEFromRpRa(rp,ra)
- RpRaFromAE(a,e)
- Ecirc(ψ,dψ,r)
- Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC)