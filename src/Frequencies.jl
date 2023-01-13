#=

wrapper to select the methodology for computing the frequencies

only selecting Henon anomaly mapping for now,
but this is where one could select for different anomalies

=#


# bring in the anomaly mapping (i.e. f(u))
include("Henon/Ufunc.jl")

# bring in the frequency mapping
include("Henon/Frequencies.jl")

# bring in the frequency inversion
include("Utils/NumericalInversion.jl")


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : Wrappers
#
########################################################################

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e,action,TOLECC,NINT,EDGE,TOLA)
wrapper to select which type of frequency computation to perform, from (a,e)
"""
function ComputeFrequenciesAE(ψ::F0,
                              dψ::F1,
                              d2ψ::F2,
                              d3ψ::F3,
                              d4ψ::F4,
                              a::Float64,
                              e::Float64,
                              TOLECC::Float64,
                              NINT::Int64,
                              EDGE::Float64,
                              TOLA::Float64)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return HenonΘFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE)

end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e[,action,TOLECC,NINT,EDGE,TOLA])
wrapper to select which type of frequency computation to perform, from (a,e)
with optional arguments
"""
function ComputeFrequenciesAE(ψ::F0,
                              dψ::F1,
                              d2ψ::F2,
                              d3ψ::F3,
                              d4ψ::F4,
                              a::Float64,
                              e::Float64;
                              TOLECC::Float64=0.001,
                              NINT::Int64=32,
                              EDGE::Float64=0.01,
                              TOLA::Float64=0.001)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE,TOLA)

end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT fourth derivative
"""
function ComputeFrequenciesAE(ψ::F0,
                              dψ::F1,
                              d2ψ::F2,
                              d3ψ::F3,
                              a::Float64,e::Float64;
                              TOLECC::Float64=0.001,
                              NINT::Int64=32,
                              EDGE::Float64=0.01,
                              FDIFF::Float64=1.e-8,
                              TOLA::Float64=0.001)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # define a numerical fourth derivative
    @inline d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE,TOLA)
end

"""ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e[,TOLECC,VERBOSE])
wrapper to select which type of frequency computation to perform, from (a,e)
EXCEPT third derivative
"""
function ComputeFrequenciesAE(ψ::F0,
                              dψ::F1,
                              d2ψ::F2,
                              a::Float64,e::Float64;
                              TOLECC::Float64=0.001,
                              NINT::Int64=32,
                              EDGE::Float64=0.01,
                              FDIFF::Float64=1.e-8,
                              TOLA::Float64=0.001)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    # define a numerical third derivative
    @inline d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF

    # Nul fourth derivative
    @inline d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE,TOLA)
end


########################################################################
#
# (a,e) -> (Ω1,Ω2) mapping : derivatives wrappers
#
########################################################################

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLECC,NINT,EDGE)
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

all parameters must be specified in this version
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,
                                       dψ::F1,
                                       d2ψ::F2,
                                       d3ψ::F3,
                                       d4ψ::F4,
                                       a::Float64,
                                       e::Float64,
                                       da::Float64,
                                       de::Float64,
                                       TOLECC::Float64,
                                       NINT::Int64,
                                       EDGE::Float64)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

        # first, check for values that need to be expanded

        # grid is structured like
        # (Ω1h,Ω2h) [+da]
        #    ^
        # (Ω1c,Ω2c)-> (Ω1r,Ω2r) [+de]

        # also need a TOLA: choose da
        # @IMPROVE make this a parameter
        TOLA = da/2

        # @IMPROVE watch out for close to TOLECC, will fail across boundary
        Ω1c,Ω2c = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE,TOLA)
        #Ω1c,Ω2c = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC,NINT,EDGE)

        # the offset in a
        a2 = a+da
        Ω1h,Ω2h = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a2,e,TOLECC,NINT,EDGE,TOLA)

        # the offset in e
        # if this is already a radial orbit, don't go to super radial
        e2 = e+de
        if e2 > 1.0
            de *= -1.0
            e2 = e+de
        end

        Ω1r,Ω2r = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e2,TOLECC,NINT,EDGE,TOLA)

        dΩ1da = (Ω1h-Ω1c)/da
        dΩ2da = (Ω2h-Ω2c)/da

        dΩ1de = (Ω1r-Ω1c)/de
        dΩ2de = (Ω2r-Ω2c)/de

        return Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLECC,NINT,EDGE)
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES

many parameters are optional in this version
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,
                                       dψ::F1,
                                       d2ψ::F2,
                                       d3ψ::F3,
                                       d4ψ::F4,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=ELTOLECC,
                                       NINT::Int64=32,
                                       EDGE::Float64=0.01)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLECC,NINT,EDGE)
end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e[,da,de,TOLECC,NINT,FDIFF,EDGE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT fourth derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,
                                       dψ::F1,
                                       d2ψ::F2,
                                       d3ψ::F3,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=ELTOLECC,
                                       NINT::Int64=32,
                                       FDIFF::Float64=1.e-8,
                                       EDGE::Float64=0.01)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # define a numerical fourth derivative
    @inline d4ψ(x::Float64) = (d3ψ(x+FDIFF)-d3ψ(x))/FDIFF

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;da,de,TOLECC,NINT,EDGE)

end

"""ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,a,e[,da,de,TOLECC,NINT,FDIFF,EDGE])
wrapper to select which type of frequency computation to perform, from (a,e), but DERIVATIVES
EXCEPT third derivative
"""
function ComputeFrequenciesAEWithDeriv(ψ::F0,
                                       dψ::F1,
                                       d2ψ::F2,
                                       a::Float64,
                                       e::Float64;
                                       da::Float64=0.0001,
                                       de::Float64=0.0001,
                                       TOLECC::Float64=ELTOLECC,
                                       NINT::Int64=32,
                                       FDIFF::Float64=1.e-8,
                                       EDGE::Float64=0.01)::Tuple{Float64,Float64,Float64,Float64,Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function}

    # define a numerical third derivative
    @inline d3ψ(x::Float64) = (d2ψ(x+FDIFF)-d2ψ(x))/FDIFF
    # Nul fourth derivative
    @inline d4ψ(x::Float64) = 0.

    return ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLECC,NINT,EDGE)
end


########################################################################
#
# (Ω1,Ω2) -> (a,e) mapping : Wrappers
#
########################################################################

"""ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,a,e[,eps,maxiter,TOLECC,TOLA])
wrapper to select which type of inversion to compute for (Omega1,Omega2)->(a,e)

all parameters must be specified
"""
function ComputeAEFromFrequencies(ψ::F0,
                                  dψ::F1,
                                  d2ψ::F2,
                                  d3ψ::F3,
                                  Ω1::Float64,Ω2::Float64,
                                  eps::Float64=1*10^(-12),
                                  maxiter::Int64=1000,
                                  TOLECC::Float64=ELTOLECC,
                                  TOLA::Float64=0.0001,
                                  da::Float64=0.0001,de::Float64=0.0001,
                                  VERBOSE::Int64=0)::Tuple{Float64,Float64} where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

        # @IMPROVE
        # use adaptive da, de branches
        # da max(0.0001,0.01a)
        # de min(max(0.0001,0.1a*e)

        a,e,_,_ = AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ;
                                                    eps=eps,ITERMAX=maxiter,
                                                    TOLECC=TOLECC,TOLA=TOLA,da=da,de=de)

        return a,e
end


########################################################################
#
# (E,L) -> (α,β) mapping : Jacobian
#
########################################################################

"""
compute the jacobian J = |d(E,L)/d(α,β)| = |d(E,L)/d(a,e)|/|d(α,β)/d(a,e)|
"""
function JacELToαβAE(ψ::F0,
                     dψ::F1,
                     d2ψ::F2,
                     d3ψ::F3,
                     d4ψ::F4,
                     a::Float64,
                     e::Float64;
                     NINT::Int64=64,
                     EDGE::Float64=0.02,
                     Ω₀::Float64=1.0,
                     TOLECC::Float64=ELTOLECC)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}


    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    #Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,TOLECC=TOLECC)
    JacELae = JacELToAE(ψ,dψ,d2ψ,a,e,TOLECC=TOLECC)

    # the (α,β) -> (a,e) Jacobian (below)
    Jacαβae = JacαβToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=NINT,EDGE=EDGE,Ω₀=Ω₀)

    # compute the Jacobian
    Jac = JacELae/Jacαβae

    # @IMPROVE: better adaptive checks here
    # do some cursory checks for quality
    if Jac < 0.0
        return 0.0
    end

    if isnan(Jac)
        return 0.0
    end

    return Jac

end

"""

@ATTENTION can use the isochrone-specific if you are using an isochrone. Otherwise this is a bit costly.


"""
function JacαβToAE(ψ::F0,
                   dψ::F1,
                   d2ψ::F2,
                   d3ψ::F3,
                   d4ψ::F4,
                   a::Float64,
                   e::Float64;
                   NINT::Int64=64,
                   EDGE::Float64=0.02,
                   Ω₀::Float64=1.0)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    # calculate the frequency derivatives
    α,β,∂α∂a,∂α∂e,∂β∂a,∂β∂e = OrbitalElements.DHenonΘFreqRatiosAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,NINT=NINT,EDGE=EDGE,Ω₀=Ω₀)

    # return the Jacobian
    Jacαβae = abs(∂α∂a*∂β∂e - ∂β∂a*∂α∂e)

end

"""

@ATTENTION this takes (a,e) as arguments.
@ATTENTION this combines several numerical derivatives; please take care!

@IMPROVE add massaging parameters for numerical derivatives
@IMPROVE fix boundary values when using limited development
@IMPROVE noisy at the boundaries

@IMPROVE, give this more derivatives!
"""
function JacELToαβAE(a::Float64,
                     e::Float64,
                     ψ::F0,
                     dψ::F1,
                     d2ψ::F2,
                     Ω₀::Float64=1.0;
                     nancheck::Bool=false,
                     NINT::Int64=64)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}

    tmpe = e
    # to be fixed for limited development...
    if e>0.99
        tmpe=0.99
    end

    if e<0.01
        #println("faking the eentricity...")
        tmpe=0.01
    end

    # get all numerical derivatives

    # these are dangerous, and break down fairly easily.
    Ω1c,Ω2c,dΩ1da,dΩ2da,dΩ1de,dΩ2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,a,tmpe,NINT=NINT)

    # @IMPROVE: use the version with derivatives specified
    # ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,da,de,TOLECC,NINT,EDGE)

    # this is nearly always safe
    # the (E,L) -> (a,e) Jacobian (in Utils/ComputeEL.jl)
    #Jac_EL_AE = JacELToAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e)
    JacELae = JacELToAE(ψ,dψ,d2ψ,a,tmpe)


    JΩ1Ω2ae = abs(dΩ1da*dΩ2de - dΩ1de*dΩ2da)

    # check for NaN or zero values
    if nancheck
        if isnan(JacELae )
            println("OrbitalElements.Frequencies.JacELToαβAE: J_EL_ae is NaN for a=$a,e=$e")
            return 0.0
        end

        if JacELae  <= 0.0
            println("OrbitalElements.Frequencies.JacELToαβAE: J_EL_ae is 0 for a=$a,e=$e")
            return 0.0
        end

        if isnan(JΩ1Ω2ae)
            println("OrbitalElements.Frequencies.JacELToαβAE: J_o12_ae is NaN for a=$a,e=$e")
            return 0.0
        end

        if JΩ1Ω2ae <= 0.0
            println("OrbitalElements.Frequencies.JacELToαβAE: J_o12_ae is 0 for a=$a,e=$e")
            return 0.0
        end
    end

    # combine and return
    return Ω1c*Ω₀*JacELae /JΩ1Ω2ae

end
