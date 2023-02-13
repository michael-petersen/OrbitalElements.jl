"""
Ufunc.jl

collection of functions that are used to compute quantities when integrating over an orbit angle.

Specific to Henon mapping.

@IMPROVE are any of these helped by function specialisation?

functions:
henon_f(u)
henon_df(u)
henon_d2f(u)
henon_d3f(u)
henon_d4f(u)
ru(u,a,e)
drdu(u,a,e)
ψeff(ψ,r,L)
dψeffdr(dψ,r,L)
d2ψeffdr2(d2ψ,r,L)
ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)
ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC)
ΘAEdade(ψ,dψ,d2ψ,d3ψ,u,a,e,da,de,TOLA,TOLECC,EDGE)

"""

########################################################################
#
# Hénon anomaly function and derivatives
#
########################################################################

"""
    henon_f(u)

the henon anomaly increment
"""
function henon_f(u::Float64)
    u*(3/2 - (u^2)/2)
end

"""
    henon_df(u)

the derivative of the henon anomaly increment
"""
function henon_df(u::Float64)
    1.5*(1.0 - u^(2))
end

"""
    henon_d2f(u)

the second derivative of the henon anomaly increment
"""
function henon_d2f(u::Float64)
    return -3u
end

"""
    henon_d3f(u)

the third derivative of the henon anomaly increment
"""
function henon_d3f(u::Float64)
    return -3.
end

"""
    henon_d4f(u)

the fourth derivative of the henon anomaly increment
"""
function henon_d4f(u::Float64)
    return 0.
end


"""
    ru(u,a,e)

mapping from u->r in Henon variables
"""
function ru(u::Float64,a::Float64,e::Float64)::Float64
    fu = henon_f(u)
    return a*(1+e*fu)
end


"""
    drdu(u,a,e)

mapping from u->r in Henon variables
"""
function drdu(u::Float64,a::Float64,e::Float64)::Float64
    dfu = henon_df(u)
    return a*e*dfu
end


########################################################################
#
# Effective potential and derivatives w.r.t. the radius
#
########################################################################

"""
    ψeff(ψ,r,L)

the effective potential: note the relationship to Q
"""
function ψeff(ψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return ψ(r)
    else
        return ψ(r) + (1/2) * (L/r)^(2)
    end
end

"""
    dψeffdr(dψ,r,L)

the derivative of the effective potential
"""
function dψeffdr(dψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return dψ(r)
    else
        return dψ(r) - (L)^(2) / (r^3)
    end
end

"""
    d2ψeffdr2(d2ψ,r,L)

the second derivative of the effective potential
"""
function d2ψeffdr2(d2ψ::Function,r::Float64,L::Float64)::Float64
    if L == 0.
        return d2ψ(r)
    else
        return d2ψ(r) + 3 * (L)^(2) / (r^4)
    end
end


"""
    ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC,EDGE)

Θ, the anomaly for computing orbit averages as a function of (a,e)
equivalent to Θ = (dr/du)(1/Vrad)
"""
function ΘAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
             u::Float64,a::Float64,e::Float64,
             TOLA::Float64,TOLECC::Float64,EDGE::Float64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # use the expanded approximation
    if ((1-abs(u))<EDGE)
        return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC)

    # if not close to the boundary, can calculate as normal
    else
        E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)

        r = ru(u,a,e)

        # this can somehow be negative: do we need an extra check?
        denomsq = 2*(E - ψeff(ψ,r,L))

        if denomsq < 0.0
            # go back to the expansion -- or should we return 0.0?
            #return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLECC,f=f,d2f=d2f,d3f=d3f,d4f=d4f)
            return ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC)
        end

        # do the standard return
        return a * e * henon_df(u) / sqrt(denomsq)
    end
end



"""
    ΘExpansionAE(ψ,dψ,d2ψ,d3ψ,u,a,e,TOLA,TOLECC)

ΘExpansion, the anomaly for computing orbit averages as a function of (a,e)

Used when u is sufficiently close to +1,-1

BIG ALLOCATIONS here from ELTOLECC not being specified.
Downside is that this guarantees allocations if TOLECC not specified.

"""
function ΘExpansionAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                      u::Float64,a::Float64,e::Float64,
                      TOLA::Float64,TOLECC::Float64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # which boundary are we close to?
    ul = (u > 0) ? 1.0 : -1.0

    # compute the corresponding radius value
    rl = ru(ul,a,e)

    # compute energy and angular momentum from the potential (allow for expansions)
    E, L = ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLA,TOLECC)

    # compute the derivatives of the effective potential
    dψeffl, d2ψeffl = dψeffdr(dψ,rl,L), d2ψeffdr2(d2ψ,rl,L)

    # compute the derivatives of the Henon f function
    d2fl, d3fl, d4fl = henon_d2f(ul), henon_d3f(ul), henon_d4f(ul)

    # define the prefactor
    pref = - ul * a * e

    # this denominator can be negative?
    combination = - a * e * dψeffl * d2fl

    # switch to safety: don't contribute anything at this point
    # In particular for radial orbits with ψ(r) = - Inf in r = 0
    # combination = - Inf
    if combination <= 0.
        return 0.0
    end

    # if >0, sqrt is safe, proceed
    denom = sqrt(combination)

    zeroorder   = d2fl
    firstorder  = d3fl / 3.0
    secondorder = (3.0 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * (d2fl)^(3)) -  dψeffl * (d3fl)^(2)) / (24.0 * dψeffl * d2fl)

    return pref / denom * ( zeroorder + firstorder * (u - ul) + secondorder * (u - ul)^(2) )
end



"""
    dΘAE(ψ,dψ,d2ψ,d3ψ,u,a,e,da,de,TOLA,TOLECC,EDGE)

numerical differentiation of Θ w.r.t. semimajor axis and eccentricity


"""
function dΘAE(ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
              u::Float64,a::Float64,e::Float64,
              da::Float64=1.0e-8,de::Float64=1.0e-8,
              TOLA::Float64=0.001,TOLECC::Float64=ELTOLECC,
              EDGE::Float64=0.01)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # Numerical derivative points
    ap, da, ep, de = NumDerivPoints(a,e,da,de,TOLA,TOLECC)

    # current point
    Θloc = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a   ,e,TOLA,TOLECC,EDGE)


    Θap  = ΘAE(ψ,dψ,d2ψ,d3ψ,u,ap,e,TOLA,TOLECC,EDGE)


    Θep  = ΘAE(ψ,dψ,d2ψ,d3ψ,u,a,ep,TOLA,TOLECC,EDGE)

    ∂Θ∂a = (Θap-Θloc)/da
    ∂Θ∂e = (Θep-Θloc)/de

    return ∂Θ∂a, ∂Θ∂e
end
