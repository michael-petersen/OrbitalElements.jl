
########################################################################
#
# Effective potential and derivatives w.r.t. the radius
#
########################################################################
"""
    _ψeff(r, L, model)

the effective potential.

@IMPROVE: not safe close to L=0 and r=0.
"""
function _ψeff(r::Float64, L::Float64, model::Potential)::Float64
    if L == 0.
        return ψ(r, model)
    end
        
    return ψ(r, model) + 0.5 * (L/r)^2
end

"""
    _dψeffdr(r, L, model)

the derivative of the effective potential.

@IMPROVE: not safe close to L=0 and r=0.
"""
function _dψeffdr(r::Float64, L::Float64, model::Potential)::Float64
    if L == 0. # Not good, just prevents errors (r=0 implies L=0)
        return dψ(r, model)
    end
        
    return dψ(r, model) - L^2 / (r^3)
end

"""
    _d2ψeffdr2(r, L, model)

the second derivative of the effective potential.

@IMPROVE: not safe close to L=0 and r=0.
"""
function _d2ψeffdr2(r::Float64, L::Float64, model::Potential)::Float64
    if L == 0. # Not good, just prevents errors (r=0 implies L=0)
        return d2ψ(r, model)
    end
        
    return d2ψ(r, model) + 3 * L^2 / (r^4)
end

########################################################################
#
# Radial velocity
#
########################################################################
"""
    _radial_velocity(u, a, e, model[, params])

radial velocity as a function of the orbital constants (a,e) and the anomaly u.

@IMPROVE: Quite naive negative/infinite velocity handling.
"""
function radial_velocity(
    u::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64

    E, L = EL_from_ae(a, e, model, params)

    r = radius_from_anomaly(u, a, e)

    vrSQ = 2(E - _ψeff(r, L, model))

    if vrSQ < 0.0 || isnan(vrSQ) || isinf(vrSQ)
        return 0.0
    end
        
    return sqrt(vrSQ)
end


########################################################################
#
# Θ(u) (frequencies integrand anomaly)
#
########################################################################
"""
    Θ(u, a, e, model[, params])

cured inverse radial velocity, Θ(u) = (dr/du)/v_rad, at anomaly `u` on orbit `(a,e)`.

@IMPROVE: right now, Hénon anomaly is hard-coded.
@IMPROVE: the expansion is trying both Taylor expansion and if fails, extrapolation with 
particular care need at boundary switch. Fix it to use only one.
@IMPROVE: find a better name
"""
function Θ(
    u::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64

    # use the expanded approximation
    # CAUTION: 1-(1-EDGE) < EDGE is true ...
    # To prevent this → EDGE - eps(Float64)
    if 1.0 - abs(u) < params.EDGE - eps(Float64)
        return _Θedge(u, a, e, model, params)
    end

    dr = radius_from_anomaly_derivative(u, a, e)
    vr = radial_velocity(u, a, e, model, params)

    if vr == 0.0
        # go back to the expansion -- or should we return 0.0?
        return 0.0
    end 
        
    return dr / vr
end

"""
    _Θedge(u, a, e, model, params)

same as `Θ(...)` for u close to +1,-1 (curing the 0/0 limit at peri/apocentre)

@IMPROVE: right now, Hénon anomaly is hard-coded.
@IMPROVE: the expansion is trying both Taylor expansion and, if fails, extrapolation with 
particular care need at boundary switch. Fix it to use only one.
@IMPROVE: find a better name
"""
function _Θedge(
    u::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Float64

    # which boundary are we close to?
    ul = (u > 0.) ? 1.0 : -1.0
    # compute the corresponding radius value
    rl = radius_from_anomaly(ul,a,e)
    # compute energy and angular momentum from the potential (allow for expansions)
    _, L = EL_from_ae(a, e, model, params)
    # compute the derivatives of the effective potential
    dψeffl, d2ψeffl = _dψeffdr(rl, L, model), _d2ψeffdr2(rl, L, model)
    # compute the derivatives of the Henon f function
    d2fl, d3fl, d4fl = _henond2f(ul), _henond3f(ul), _henond4f(ul)
    # define the prefactor
    pref = - ul * a * e
    # this denominator can be negative?
    combination = - a * e * dψeffl * d2fl

    # switch to safety: don't contribute anything at this point
    # In particular for radial orbits with ψ(r) = - Inf in r = 0
    # combination = - Inf
    if combination <= 0.
        # Switch to
        # extension of the function close to the border u ~ ul using the
        # linear interpolation between using the points ul ± tolu and ul ± 2*tolu
        u1 = ul * (1 - params.EDGE)
        u2 = ul * (1 - 2*params.EDGE)
        u3 = ul * (1 - 3*params.EDGE)

        Θ1 = Θ(u1, a, e, model, params)
        Θ2 = Θ(u2, a, e, model, params)
        Θ3 = Θ(u3, a, e, model, params)

        # It in fact is an extrapolation
        return _interpolate_order_2(u, u1, Θ1, u2, Θ2, u3, Θ3)
    end

    # if >0, sqrt is safe, proceed
    denom = sqrt(combination)

    zeroorder = d2fl
    firstorder = d3fl / 3
    secondorder = (
        (
            3 * (dψeffl * d2fl * d4fl - a * e * d2ψeffl * d2fl^3) 
            -  dψeffl * d3fl^2
        ) 
        / (24 * dψeffl * d2fl)
    )

    return pref / denom * (zeroorder + firstorder * (u - ul) + secondorder * (u - ul)^2)
end


########################################################################
#
# Θ(u) derivatives w.r.t. a and e
#
########################################################################

"""
    _Θ_derivatives_ae(u, a, e, modelparams)

numerical differentiation of Θ w.r.t. semimajor axis and eccentricity

@IMPROVE: right now the derivative procedure is hard-coded. Should incorparate it 
in the OrbitalParameters structure (not only the steps)
@IMPROVE: find a better name
@IMPROVE: no used for now as frequency derivatives are taken through finite difference.
"""
function _Θ_derivatives_ae(
    u::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters=OrbitalParameters()
)::Tuple{Float64,Float64}

    # Function to differentiate
    fun(atmp::Float64, etmp::Float64) = Θ(u, atmp, etmp, model, params)
    # Perform differentiation
    _, ∂Θ∂a, ∂Θ∂e = _derivatives_ae(fun, a, e, params)

    return ∂Θ∂a, ∂Θ∂e
end