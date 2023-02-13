"""
check a bunch of quantities against the isochrone case, so we can be confident we are doing the numerical work correctly!
"""


import OrbitalElements
using Printf

# define easy potentials to pass to frequency calculators
const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64    = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64   = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
d3ψ(r::Float64)::Float64  = OrbitalElements.d3ψIsochrone(r,bc,M,G)
d4ψ(r::Float64)::Float64  = OrbitalElements.d4ψIsochrone(r,bc,M,G)

Ω₀ = OrbitalElements.Ω₀Isochrone(bc,M,G)

# params to specify
da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε,eps = 0.001,0.001,0.001,32,0.01,0.001,100,1.e-6,1.e6,1.e-10,1.e-12


x = -1.0
n1,n2 = 1,2
vc = n1*OrbitalElements.Ω1circular(dψ,d2ψ,x) + n2*OrbitalElements.Ω2circular(dψ,x)

println(vc)

# select an (a,e) value for the orbit
a,e = 0.05, 1.0

# compute rperi and rapo
rp,ra = OrbitalElements.RpRaFromAE(a,e)
println("rp=$rp ra=$ra",)

# test frequency computation
Ω₁e,Ω₂e = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
Jrr = OrbitalElements.IsochroneJrRpRa(rp,ra,bc,M,G)
#println("Ω₁r=$Ω₁r,Ω₂r=$Ω₂r")

#rcirc1 = OrbitalElements.Omega1circ_to_radius(Ω₁e,dψ,d2ψ)#,Ziter=32,verbose=false)
#println("Ω₁ Bisect r=$rcirc1")

#rcirc0 = OrbitalElements.Omega2circ_to_radius(Ω₂e,dψ)
#println("Ω₂ Bisect r=$rcirc1")

println("truth Ω₁=$Ω₁e,Ω₂=$Ω₂e")

# make a HIGH RES version of the frequencies
Ω₁r,Ω₂r = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
println("oldae O1=$Ω₁r O2=$Ω₂r")

@time Ec,Lc,dEda,dEde,dLda,dLde = OrbitalElements.dELFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)
Em,Lm = OrbitalElements.IsochroneELFromAE(a,e,bc,M,G)
#Em /= OrbitalElements.isochrone_E0()
println("estimated E=$Ec,L=$Lc")
println("true      E=$Em,L=$Lm")


# plot Omega1,Omega2 vs e for the analytic and expansion cases
# is there a problem in the expansion?
# any discontinuities from switching to beta

#=
alpha,beta = Ω₁c/Ω₀,Ω₂c/Ω₁c
println("alpha=$alpha,beta=$beta")

J_EL_ab = OrbitalElements.IsochroneJacELtoAlphaBeta(alpha,beta,bc,M,G)
println("Jacobian(EL,ab):$J_EL_ab")

J_EL_abT = OrbitalElements.JacELToAlphaBetaAE(a,e,ψ,dψ,d2ψ)
println("TJacobian(EL,ab):$J_EL_abT")
=#


# get the numerical frequency derivatives at this point
f1c,f2c,df1da,df2da,df1de,df2de = OrbitalElements.ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)


f1h,f2h,df1dah,df1deh,df2dah,df2deh = OrbitalElements.DFrequenciesHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)


# check isochrone numerical diff for frequencies
da = 1.e-6
de = 1.e-6

if e+de > 1.0
    de *= -1.0
end
f1m,f2m = OrbitalElements.IsochroneOmega12FromAE(a,e,bc,M,G)
f1a,f2a = OrbitalElements.IsochroneOmega12FromAE(a+da,e,bc,M,G)
f1e,f2e = OrbitalElements.IsochroneOmega12FromAE(a,e+de,bc,M,G)
df1da2,df1de2,df2da2,df2de2 = (f1a-f1m)/da,(f1e-f1m)/de,(f2a-f2m)/da,(f2e-f2m)/de

println("Compare derivatives:")
println("NDiff : df1da=$df1da,df2da=$df2da,df1de=$df1de,df2de=$df2de")
println("DΘ    : df1da=$df1dah,df2da=$df2dah,df1de=$df1deh,df2de=$df2deh")
println("Truth : df1da=$df1da2,df2da=$df2da2,df1de=$df1de2,df2de=$df2de2")
