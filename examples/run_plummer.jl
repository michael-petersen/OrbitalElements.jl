"""
In contrast to the isochrone case, we do not have ground truth for Plummer. So all comparisons will be made against high-resolution quantities.
"""


import OrbitalElements

# define easy potentials to pass to frequency calculators
bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
d3ψ(r::Float64)::Float64 = OrbitalElements.d3ψPlummer(r,bc,M,G)
d4ψ(r::Float64)::Float64 = OrbitalElements.d4ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

println("Central pot=$(ψ(0.0))")
println("Central dpot=$(dψ(0.0))")
println("Central d2pot=$(d2ψ(0.0))")

rmin,rmax = 0.,1000000.0
βc(αc::Float64)::Float64 = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin,rmax)

αmin,αmax = OrbitalElements.αminmax(dψ,d2ψ,rmin,rmax,Ω₀)
println("(αmin,αmax)=($αmin,$αmax)")

n1,n2 = 1,-1
ωmin,ωmax = OrbitalElements.Findωminωmax(n1,n2,dψ,d2ψ,αmin,αmax,Ω₀,rmin,rmax)
println("(ωmin,ωmax)=($ωmin,$ωmax)")

uval = 0.5
vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀,rmin,rmax)
println("(vmin,vmax)=($vmin,$vmax)")

vval = 0.5*(vmax-vmin) + vmin


α,β = OrbitalElements.αβFromUV(uval,vval,n1,n2,ωmin,ωmax)
println("(u,v)=($uval,$vval)")
println("(α,β)=($α,$β)")

Ω1,Ω2 = α*Ω₀,α*β*Ω₀
# (Ω1,Ω2) -> (a,e)
da,de,TOLECC,NINT,EDGE,TOLA,ITERMAX,rmin,rmax,invε,eps = 0.001,0.001,0.001,32,0.01,0.001,100,1.e-6,1.e6,1.e-10,1.e-12
a,e = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε,Ω₀)
println("(a,e)=($a,$e)")


# select an (a,e) value for the orbit
a,e = 0.1, 0.5
rp,ra = OrbitalElements.RpRaFromAE(a,e)
println("rp=$rp,ra=$ra")

sp,sa = OrbitalElements.SpSaFromRpRa(rp,ra,bc)
println("sp=$sp,sa=$sa")

Eval,Lval = OrbitalElements.PlummerELFromSpSa(sp, sa, bc=bc ,M=M,G=G)
println("anomaly E=$Eval,L=$Lval")

spo,sao = OrbitalElements.SpSaFromEL(Eval,Lval, bc=bc ,M=M,G=G)
println("sp=$spo,sa=$sao")


Eval,Lval = OrbitalElements.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e)
println("brute   E=$Eval,L=$Lval")

u=-1.0
testrp = OrbitalElements.RFromURpRa(u,rp,ra,bc)
println("rp=$rp,rp=$testrp")

O1,O2,Jr = OrbitalElements.ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
println("O1=$O1,O2=$O2,Jr=$Jr")
O1,O2,Jr = OrbitalElements.PlummerOmega12FromRpRa(rp,ra,bc,M,G,action=true)
println("O1=$O1,O2=$O2,Jr=$Jr")

α,β = OrbitalElements.PlummerAlphaBetaFromRpRa(rp,ra,bc,M,G)
println("α=$α,β=$β")


dJrdE, dJrdL, d2JrdE2, d2JrdEL, d2JrdL2 = OrbitalElements.GradJrELWrap(Eval,Lval,bc=bc,M=M,G=G)

Etest,Ltest = OrbitalElements.ELFromAlphaBeta(α,β,bc=bc,M=M,G=G)
println("output  E=$Etest,L=$Ltest")
