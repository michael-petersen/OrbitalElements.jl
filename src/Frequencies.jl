"""

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

"""


# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")

# bring in the frequency mapping
include("Henon/Frequencies.jl")

# bring in the frequency inversion
include("Utils/NumericalInversion.jl")

########################################################################
#
# (α,β) -> (Ω1,Ω2) mapping
#
########################################################################
"""
    αβFromFrequencies(Ω1,Ω2,Ω₀)

converts frequencies to frequencies ratios
"""
function αβFromFrequencies(Ω1::Float64,Ω2::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω1/Ω₀, Ω2/Ω1
end

"""
    FrequenciesFromαβ(α,β,Ω₀)

converts frequencies ratios to frequencies
"""
function FrequenciesFromαβ(α::Float64,β::Float64,
                           Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*α, Ω₀*α*β
end

"""
    FrequenciesDerivsFromαβDerivs(α,β,dα,dβ,Ω₀)

converts frequencies ratios derivatives to frequencies derivatives.
"""
function FrequenciesDerivsFromαβDerivs(α::Float64,β::Float64,
                                       dα::Float64,dβ::Float64,
                                       Ω₀::Float64)::Tuple{Float64,Float64}

    return Ω₀*dα, Ω₀*(dα*β + α*dβ) 
end

########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : Wrappers
#
########################################################################
"""
    αβFromAE(ψ,dψ,d2ψ,a,e,params)
"""
function αβFromAE(ψ::F0,dψ::F1,d2ψ::F2,
                  a::Float64,e::Float64,
                  params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    return αβHenonΘAE(ψ,dψ,d2ψ,a,e,params)
end

"""
    ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeFrequenciesAE(ψ::F0,dψ::F1,d2ψ::F2,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    α, β = αβFromAE(ψ,dψ,d2ψ,a,e,params)
    return FrequenciesFromαβ(α,β,params.Ω₀)
end




"""
    ComputeFrequenciesJAE(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeFrequenciesJAE(ψ::F0,dψ::F1,d2ψ::F2,
                               a::Float64,e::Float64,
                               params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    Ω1, Ω2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,params)
    J = HenonJFromAE(ψ,dψ,d2ψ,a,e,params)

    return Ω1, Ω2, J
end

########################################################################
#
# (a,e) -> (Ω1,Ω2), (α,β) mapping : derivatives wrappers
#
########################################################################

"""
    ComputeαβWithDerivAE(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeαβWithDerivAE(ψ::F0,dψ::F1,d2ψ::F2,
                              a::Float64,e::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    # Numerical derivative points
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    ap, da, ep, de = NumericalDerivativePoints(a,e,params.da,params.de,tola,tole)

    # Derivation outside the integral
    α, β = αβFromAE(ψ,dψ,d2ψ,a,e,params)

    # For a derivatives
    αap, βap = αβFromAE(ψ,dψ,d2ψ,ap,e,params)
    
    # For e derivatives
    αep, βep = αβFromAE(ψ,dψ,d2ψ,a,ep,params)

    ∂α∂a = (αap-α)/da
    ∂β∂a = (βap-β)/da

    ∂α∂e = (αep-α)/de
    ∂β∂e = (βep-β)/de

    return α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e
end

"""
    ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,
                                       a::Float64,e::Float64,
                                       params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}
    
    Ω₀ = params.Ω₀
    α, β, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,a,e,params)

    Ω1, Ω2          = FrequenciesFromαβ(α,β,Ω₀)
    ∂Ω1∂a, ∂Ω2∂a    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂a,∂β∂a,Ω₀)
    ∂Ω1∂e, ∂Ω2∂e    = FrequenciesDerivsFromαβDerivs(α,β,∂α∂e,∂β∂e,Ω₀)

    return Ω1, Ω2, ∂Ω1∂a, ∂Ω2∂a, ∂Ω1∂e, ∂Ω2∂e
end


########################################################################
#
# (a,e) -> (J,L) mapping : Wrappers
#
########################################################################

"""
    ComputeActionsAE(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeActionsAE(ψ::F0,dψ::F1,d2ψ::F2,
                          a::Float64,e::Float64,
                          params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    J = HenonJFromAE(ψ,dψ,d2ψ,a,e,params)
    L = LFromAE(ψ,dψ,a,e,params)
    return J, L
end

########################################################################
#
# (a,e) -> (J,L) mapping : derivatives wrappers
#
########################################################################

"""
    ComputeActionsAEWithDeriv(ψ,dψ,d2ψ,a,e,params)
"""
function ComputeActionsAEWithDeriv(ψ::F0,dψ::F1,d2ψ::F2,
                                   a::Float64,e::Float64,
                                   params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    # Numerical derivative points
    tola, tole = params.TOLA, EccentricityTolerance(a,params.rc,params.TOLECC)
    ap, da, ep, de = NumericalDerivativePoints(a,e,params.da,params.de,tola,tole)

    # Derivation outside the integral
    J, L = ComputeActionsAE(ψ,dψ,d2ψ,a,e,params)

    # For a derivatives
    Jap, Lap = ComputeActionsAE(ψ,dψ,d2ψ,ap,e,params)
    
    # For e derivatives
    Jep, Lep = ComputeActionsAE(ψ,dψ,d2ψ,a,ep,params)

    ∂J∂a = (Jap-J)/da
    ∂L∂a = (Lap-L)/da

    ∂J∂e = (Jep-J)/de
    ∂L∂e = (Lep-L)/de

    return J, L, ∂J∂a, ∂L∂a, ∂J∂e, ∂L∂e
end


########################################################################
#
# (Ω1,Ω2) -> (a,e) mapping : Wrappers
#
########################################################################

"""
    ComputeAEFromFrequencies(ψ,dψ,d2ψ,Ω1,Ω2,params)
"""
function ComputeAEFromFrequencies(ψ::F0,dψ::F1,d2ψ::F2,
                                  Ω1::Float64,Ω2::Float64,
                                  params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

        a, e, _, _ = AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,params)

        return a, e
end


"""
    ComputeAEFromActions(ψ,dψ,d2ψ,J,L,params)
"""
function ComputeAEFromActions(ψ::F0,dψ::F1,d2ψ::F2,
                              J::Float64,L::Float64,
                              params::OrbitalParameters=OrbitalParameters())::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

        a, e, _, _ = AEFromJLBrute(J,L,ψ,dψ,d2ψ,params)

        return a, e
end



########################################################################
#
# (E,L) -> (α,β) mapping : Jacobian
#
########################################################################

"""
    JacELToαβAE(ψ,dψ,d2ψ,a,e,da,de,params)

Jacobian of the (α,β) ↦ (E,L) mapping, i.e. |∂(E,L)/∂(α,β)| = |∂(E,L)/∂(a,e)| / |∂(α,β)/∂(a,e)|
"""
function JacαβToELAE(ψ::F0,dψ::F1,d2ψ::F2,
                     a::Float64,e::Float64,
                     params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}


    # the (a,e) -> (E,L) Jacobian (in Utils/ComputeEL.jl)
    Jac_AE_To_EL = JacAEToEL(ψ,dψ,d2ψ,a,e,params)

    # the (a,e) -> (α,β) Jacobian (below)
    Jac_AE_To_αβ = JacAEToαβ(ψ,dψ,d2ψ,a,e,params)

    # compute the Jacobian
    Jac_αβ_To_EL = Jac_AE_To_EL/Jac_AE_To_αβ

    # do some cursory checks for quality
    if (Jac_αβ_To_EL < 0.0) || isinf(Jac_αβ_To_EL) || isnan(Jac_αβ_To_EL)
        Jac_αβ_To_EL = 0.0
    end

    return Jac_αβ_To_EL
end

"""
    JacαβToAE(ψ,dψ,d2ψ,a,e,params)
Jacobian of the (a,e) ↦ (α,β) mapping, i.e. |∂(α,β)/∂(a,e)|
"""
function JacAEToαβ(ψ::F0,dψ::F1,d2ψ::F2,
                   a::Float64,e::Float64,
                   params::OrbitalParameters=OrbitalParameters())::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}

    # calculate the frequency derivatives
    _, _, ∂α∂a, ∂β∂a, ∂α∂e, ∂β∂e = ComputeαβWithDerivAE(ψ,dψ,d2ψ,a,e,params)

    # return the Jacobian
    return abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)
end