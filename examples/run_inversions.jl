"""
test some basic inversions to make sure all definitions are equivalent

julia --compile=min run_inversions.jl

"""

using BenchmarkTools
import OrbitalElements

# define easy potentials to pass to frequency calculators
#bc, M, G = 1.,1. ,1.
const bc, M, G = 1.,1. ,1.  # these can be constant for optimization flags

ψ(r::Float64)::Float64   = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψIsochrone(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)

# select an (a,e) value for the orbit
a,e = 100.0, 0.9
#a,e = 2.6636008542719694,0.9824936963436619
#a,e = 100.,0.001
println("Input     a=$a e=$e")

da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε,eps = 0.001,0.001,0.001,32,0.01,0.001,20,1.e-6,1.e6,1.e-10,1.e-12

# compute rperi and rapo
#println("Compute rp,ra...")
#rp,ra = OrbitalElements.RpRaFromAE(a,e)
#println("rp=$rp ra=$ra")

# compute exact frequencies
Ω₁e,Ω₂e = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
println("Exact       Ω₁=$Ω₁e, Ω₂=$Ω₂e")

# compute approximate frequencies
@time Ω₁a,Ω₂a = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
println("Approximate Ω₁=$Ω₁a, Ω₂=$Ω₂a")

@time α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = OrbitalElements.ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
println("$α, $β, $∂α∂a, $∂β∂a, $∂α∂e, $∂β∂e")

@time f1,f2,df1da,df2da,df1de,df2de = OrbitalElements.ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
println("$f1,$f2,$df1da,$df2da,$df1de,$df2de")

# try an inversion back to (a,e) using third derivative
@time aa,ea,ITERMAX,tol = OrbitalElements.AEFromΩ1Ω2Brute(Ω₁a,Ω₂a,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε,Ω₀)
println("Recovered a=$aa e=$ea")
println("with $ITERMAX iterations and $tol tolerance")


Je,Le = OrbitalElements.IsochroneActionsFromAE(a,e,bc,M,G)
println("Exact       Jr=$Ω₁e, L=$Ω₂e")

println("Lcirc1=",OrbitalElements.Lcirc(dψ,a))
println("Lcircmin=",OrbitalElements.Lcirc(dψ,1.e-6))
println("Lcircmax=",OrbitalElements.Lcirc(dψ,1.e6))

Lt = a * (1-(e)^(2)) * sqrt( (ψ(a*(1+e)) - ψ(a*(1-e))) / (2e) )
println("Ltest=$Lt")

Lv = OrbitalElements.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
println("Lv=$Lv")

tolx=1000.0*eps
tolf=1000.0*eps
bs = OrbitalElements.bisection(r -> Lv - OrbitalElements.Lcirc(dψ,r),rmin,rmax,tolx=tolx,tolf=tolf)

println("rcirc",OrbitalElements.RcircFromL(Lv,dψ,rmin,rmax))


@time Ja,La = OrbitalElements.ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)
println("Approximate Jr=$Ja, L=$La")


@time aa,ea,ITERMAX,tol = OrbitalElements.AEFromJLBrute(Ja,La,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)
println("Recovered a=$aa e=$ea")
println("with $ITERMAX iterations and $tol tolerance")
