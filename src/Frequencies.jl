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
    αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
"""
function αβFromAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                  a::Float64,e::Float64,
                  TOLA::Float64,TOLECC::Float64,
                  NINT::Int64,EDGE::Float64,Ω₀::Float64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return αβHenonΘAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
end

function αβFromAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                             a::Float64,e::Float64,
                             params::OrbitsParameters)::Tuple{Float64,Float64}

    return αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end


"""
    ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
    ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params)
"""
function ComputeFrequenciesAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    α,β = αβFromAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
    return FrequenciesFromαβ(α,β,Ω₀)
end

function ComputeFrequenciesAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64}

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end



"""
    ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,TOLA,Ω₀)
"""
function ComputeFrequenciesJAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64)::Tuple{Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    Ω1, Ω2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC,NINT,EDGE,Ω₀)
    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLA,TOLECC)

    return Ω1, Ω2, J
end


function ComputeFrequenciesJAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                               a::Float64,e::Float64,
                               params::OrbitsParameters)::Tuple{Float64,Float64,Float64}

    return ComputeFrequenciesJAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.TOLA)
end

"""
    ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
"""
function ComputeαβWithDerivAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              da::Float64,de::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              NINT::Int64,EDGE::Float64,Ω₀::Float64)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

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

function ComputeαβWithDerivAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              a::Float64,e::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.da,params.de,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end


"""
    ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                       a::Float64,e::Float64,
                                       da::Float64,de::Float64,
                                       TOLA::Float64,TOLECC::Float64,
                                       NINT::Int64,EDGE::Float64,Ω₀::Float64)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end

function ComputeFrequenciesAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                       a::Float64,e::Float64,
                                       params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.da,params.de,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end


"""
    ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC,NINT)
"""
function ComputeActionsAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                          a::Float64,e::Float64,
                          TOLA::Float64,TOLECC::Float64,
                          NINT::Int64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    J = HenonJFromAE(ψ,dψ,d2ψ,d3ψ,a,e,NINT,TOLA,TOLECC)
    L = LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)
    return J, L
end


function ComputeActionsAE(ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                          a::Float64,e::Float64,
                          params::OrbitsParameters)::Tuple{Float64,Float64}

    return ComputeActionsAE(ψ,dψ,d2ψ,d3ψ,a,e,params.TOLA,params.TOLECC,params.NINT)
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

function ComputeActionsAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                                   a::Float64,e::Float64,
                                   params::OrbitsParameters)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    return ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e,params.da,params.de,params.TOLA,params.TOLECC,params.NINT)
end



"""
    ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε,Ω₀)
"""
function ComputeAEFromFrequencies(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                  Ω1::Float64,Ω2::Float64,
                                  ITERMAX::Int64,
                                  da::Float64,de::Float64,
                                  TOLA::Float64,TOLECC::Float64,
                                  EDGE::Float64,NINT::Int64,
                                  rmin::Float64,rmax::Float64,
                                  invε::Float64,
                                  Ω₀::Float64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

        # @IMPROVE
        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)
        a,e,_,_ = AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε,Ω₀)

        return a,e
end

function ComputeAEFromFrequencies(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                                  Ω1::Float64,Ω2::Float64,
                                  params::OrbitsParameters)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

        return ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1,Ω2,params.ITERMAX,params.da,params.de,params.TOLA,params.TOLECC,params.EDGE,params.NINT,params.rmin,params.rmax,params.invε,params.Ω₀)
end

"""
    ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,d4ψ,J,L,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)
"""
function ComputeAEFromActions(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              J::Float64,L::Float64,
                              ITERMAX::Int64,
                              da::Float64,de::Float64,
                              TOLA::Float64,TOLECC::Float64,
                              EDGE::Float64,NINT::Int64,
                              rmin::Float64,rmax::Float64,
                              invε::Float64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

        # @IMPROVE
        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)
        a,e,_,_ = AEFromJLBrute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,d4ψ,ITERMAX,da,de,TOLA,TOLECC,EDGE,NINT,rmin,rmax,invε)

        return a,e
end

function ComputeAEFromActions(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                              J::Float64,L::Float64,
                              params::OrbitsParameters)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

        return ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,d4ψ,J,L,params.ITERMAX,params.da,params.de,params.TOLA,params.TOLECC,params.EDGE,params.NINT,params.rmin,params.rmax,params.invε)
end


"""
    JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

compute the jacobian J = |d(E,L)/d(α,β)| = |d(E,L)/d(a,e)|/|d(α,β)/d(a,e)|
"""
function JacELToαβAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                     a::Float64,e::Float64,
                     da::Float64,de::Float64,
                     TOLA::Float64,TOLECC::Float64,
                     NINT::Int64,EDGE::Float64,Ω₀::Float64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    JacELae = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLA,TOLECC)

    # the (α,β) -> (a,e) Jacobian (below)
    Jacαβae = JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

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

function JacELToαβAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                     a::Float64,e::Float64,
                     params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.da,params.de,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end



"""
    JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

"""
function JacαβToAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64,
                   da::Float64,de::Float64,
                   TOLA::Float64,TOLECC::Float64,
                   NINT::Int64,EDGE::Float64,Ω₀::Float64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # calculate the frequency derivatives
    α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e = OrbitalElements.ComputeαβWithDerivAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLA,TOLECC,NINT,EDGE,Ω₀)

    # return the Jacobian
    Jacαβae = abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)

end

function JacαβToAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,
                   a::Float64,e::Float64,
                   params::OrbitsParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # return the Jacobian
    return JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,params.da,params.de,params.TOLA,params.TOLECC,params.NINT,params.EDGE,params.Ω₀)
end
