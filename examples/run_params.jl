

using BenchmarkTools

using OrbitalElements

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64    = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64  = OrbitalElements.d4ψIsochrone(r,bc,M,G)

Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)


a,e = 1.,1.0



da,de,TOLA,TOLECC,NINT,EDGE = 0.0001,0.0001,0.001,0.001,32,0.01

#@benchmark OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)



const O = OrbitalElements.OrbitsParametersCreate(dψ,d2ψ,Ω₀)
O2 = OrbitalElements.OrbitsParametersCreate(dψ,d2ψ,Ω₀)

#@benchmark OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,O2.TOLA,O2.TOLECC,O2.NINT,O2.EDGE)

#@benchmark OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,O.TOLA,O.TOLECC,O.NINT,O.EDGE)

da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε,eps = 0.0001,0.0001,0.001,32,0.01,0.001,6,1.e-6,1.e6,1.e-8,1.e-12
Ω₁,Ω₂ = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)

Ω₁,Ω₂ = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,O)


#@benchmark OrbitalElements.AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)

gooda,goode,iter,tolval = OrbitalElements.AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)

# test out some edge cases
Ω₁,Ω₂ = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,0.0,TOLA,TOLECC,NINT,EDGE,Ω₀)
#@benchmark OrbitalElements.AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)
# Memory estimate: 77.41 KiB, allocs estimate: 1316. (338.383 μs)

Ω₁,Ω₂ = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,1.0,TOLA,TOLECC,NINT,EDGE,Ω₀)
#@benchmark OrbitalElements.AEFromΩ1Ω2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)
# Memory estimate: 8.50 KiB, allocs estimate: 146. (40.811 μs)

n1,n2=2,2
u = 0.0
#OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,αmin,αmax,Ω₀,rmin,rmax)
ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,O)
#OrbitalElements.FindVminVmax(u,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀,rmin,rmax)
βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin,rmax)
OrbitalElements.FindVminVmax(u,n1,n2,dψ,d2ψ,ωmin,ωmax,βc,O)
