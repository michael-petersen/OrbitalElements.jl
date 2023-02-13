"""

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

functions:
αβFromFrequencies(Ω1,Ω2,Ω₀)
FrequenciesFromαβ(α,β,Ω₀)
FrequenciesDerivsFromαβDerivs(α,β,dα,dβ,Ω₀)
αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,VERBOSE,NINT,EDGE,Ω₀)
ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)
ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e,da,de,TOLA,TOLECC,NINT)
ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,e,eps,maxiter,da,de,TOLA,TOLECC)
JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,NINT,Ω₀)
JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT,EDGE,TOLA)

anything that takes theta needs EDGE

"""


# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")

# bring in the frequency mapping
include("Henon/Frequencies.jl")

# bring in the frequency inversion
include("Utils/NumericalInversion.jl")

"""
    αβFromFrequencies(Ω1,Ω2,Ω₀)
"""
function αβFromFrequencies(Ω1::Float64,Ω2::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω1/Ω₀, Ω2/Ω1
end

"""
    FrequenciesFromαβ(α,β,Ω₀)
"""
function FrequenciesFromαβ(α::Float64,β::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*α, Ω₀*α*β
end

"""
    FrequenciesDerivsFromαβDerivs(α,β,dα,dβ,Ω₀)
"""
function FrequenciesDerivsFromαβDerivs(α::Float64,β::Float64,
                                    dα::Float64,dβ::Float64,
                                    Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*dα, Ω₀*(dα*β + α*dβ)
end


"""
    αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
"""
function αβFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                  a::Float64,e::Float64,
                  TOLA::Float64,TOLECC::Float64,
                  NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
end


"""
    ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE)
"""
function ComputeFrequenciesAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    α,β = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
    return FrequenciesFromαβ(α,β,Ω₀)
end

"""
    ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,TOLA)
"""
function ComputeFrequenciesJAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    Ω1, Ω2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLECC)

    return Ω1, Ω2, J
end

"""
    ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
"""
function ComputeαβWithDerivAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              da::Float64,de::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,da,de,TOLA,TOLECC)

    # Derivation outside the integral
    α, β = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)

    # For a derivatives
    αap, βap = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,ap,e,TOLA,TOLECC,NINT,EDGE,Ω₀)

    # For e derivatives
    αep, βep = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,ep,TOLA,TOLECC,NINT,EDGE,Ω₀)

    ∂α∂a = (αap-α)/da
    ∂β∂a = (βap-β)/da

    ∂α∂e = (αep-α)/de
    ∂β∂e = (βep-β)/de

    return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
end

"""
    ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                       a::Float64,e::Float64,
                                       da::Float64,de::Float64,
                                       TOLA::Float64,TOLECC::Float64,
                                       NINT::Int64,EDGE::Float64,Ω₀::Float64=1.0)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end


"""
    ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)
"""
function ComputeActionsAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                          a::Float64,e::Float64,
                          TOLA::Float64,TOLECC::Float64,
                          NINT::Int64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLECC)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
    return J, L
end


"""
    ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e,da,de,TOLA,TOLECC,NINT)
"""
function ComputeActionsAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                                   a::Float64,e::Float64,
                                   da::Float64,de::Float64,
                                   TOLA::Float64,TOLECC::Float64,
                                   NINT::Int64)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}


        # first, check for values that need to be expanded

    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,da,de,TOLA,TOLECC)

    # Derivation outside the integral
    J, L = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)

    # For a derivatives
    Jap, Lap = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,ap,e,TOLA,TOLECC,NINT)

    # For e derivatives
    Jep, Lep = ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,ep,TOLA,TOLECC,NINT)

    ∂J∂a = (Jap-J)/da
    ∂L∂a = (Lap-L)/da

    ∂J∂e = (Jep-J)/de
    ∂L∂e = (Lep-L)/de

    return J, L, ∂J∂a, ∂L∂a, ∂J∂e, ∂L∂e
end



"""
    ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,e,eps,maxiter,da,de,TOLA,TOLECC)
"""
function ComputeAEFromFrequencies(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                                  Ω1::Float64,Ω2::Float64,
                                  eps::Float64=1*10^(-12),
                                  maxiter::Int64=1000,
                                  da::Float64=0.0001,de::Float64=0.0001,
                                  TOLA::Float64=0.0001,TOLECC::Float64=ELTOLECC)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

        # @IMPROVE
        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,_,_ = AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ;
                                                    eps=eps,ITERMAX=maxiter,
                                                    TOLECC=TOLECC,TOLA=TOLA,da=da,de=de)

        return a,e
end



"""
    JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,NINT,Ω₀)

compute the jacobian J = |d(E,L)/d(α,β)| = |d(E,L)/d(a,e)|/|d(α,β)/d(a,e)|
"""
function JacELToαβAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                     a::Float64,e::Float64,
                     TOLA::Float64=ELTOLECC,TOLECC::Float64=ELTOLECC,
                     NINT::Int64=64,EDGE::Float64=0.02,Ω₀::Float64=1.0)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    JacELae = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)

    # the (α,β) -> (a,e) Jacobian (below)
    Jacαβae = JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT,EDGE,TOLA)

    # compute the Jacobian
    Jac = JacELae/Jacαβae

    # @IMPROVE: better adaptive checks here
    # do some cursory checks for quality
    if Jac < 0.0
        return 0.0
    end

    # does this throw an allocation?
    if isnan(Jac)
        return 0.0
    end

    return Jac

end

"""
    JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT,EDGE,TOLA)

"""
function JacαβToAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64;
                   NINT::Int64=64,EDGE::Float64=0.02,Ω₀::Float64=1.0)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # calculate the frequency derivatives
    α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e = OrbitalElements.DHenonΘFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=NINT,EDGE=EDGE,Ω₀=Ω₀)

    # return the Jacobian
    Jacαβae = abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)

end
