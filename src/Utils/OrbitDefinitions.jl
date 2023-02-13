"""
basic orbit transformations

"""

"""
    AEFromRpRa(rp,ra)

function to translate pericentre and apocentre to semi-major axis and eccentricity

"""
function AEFromRpRa(rp::Float64,ra::Float64)::Tuple{Float64,Float64}

    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

"""
    RpRafromAE(a,e)

function to translate semi-major axis and eccentricity to pericentre and apocentre

"""
function RpRaFromAE(a::Float64,e::Float64)::Tuple{Float64,Float64}

    return a*(1-e),a*(1+e)
end




"""
    Ecirc(ψ,dψ,r)

compute the energy of a circular orbit at some radius
must define the potential and potential derivative a priori

"""
function Ecirc(ψ::Function,dψ::Function,r::Float64)::Float64

    return  ψ(r) + 0.5*r*dψ(r)
end


"""
    Vrad(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC)

vr, radial velocity for computing action
as a function of (a,e)

used in action computation

@IMPROVE this still has some allocations associated with it if no @inline.
"""
function Vrad(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
              u::Float64,a::Float64,e::Float64,
              TOLECC::Float64=ELTOLECC)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC)

    r = ru(u,a,e)

    vrSQ = 2*(E - ψeff(ψ,r,L))

    if (vrSQ < 0.0) || isnan(vrSQ) || isinf(vrSQ)
        return 0.0
    else
        return sqrt(vrSQ)
    end
end
