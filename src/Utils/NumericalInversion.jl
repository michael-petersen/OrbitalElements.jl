"""
numerical inversion of (Ω₁,Ω₂) -> (a,e)
by brute-forcing the derivative increments dΩ₁/da, dΩ₁/de, deΩ₂/da, dΩ₂/de

"""




"""AEFromOmega1Omega2Brute(Ω₁,Ω₂,ψ,dψ,d2ψ,d3ψ[,eps,maxiter,TOLECC,TOLA,da,de,verbose])

basic Newton-Raphson algorithm to find (a,e) from (Ω₁,Ω₂) brute force derivatives.

@IMPROVE add escape for circular orbits
"""
function AEFromOmega1Omega2Brute(Ω₁::Float64,Ω₂::Float64,
                                 ψ::Function,
                                 dψ::Function,
                                 d2ψ::Function,
                                 d3ψ::Function;
                                 eps::Float64=1*10^(-10),
                                 ITERMAX::Int64=1000,
                                 TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                 da::Float64=1.0e-5,de::Float64=1.0e-5,
                                 verbose::Int64=0,
                                 EDGE::Float64=0.03,
                                 NINT::Int64=64)
    #

    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = Omega1circ_to_radius(Ω₁,dψ,d2ψ)

    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5
    f1,f2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,aguess,eguess,NINT=NINT,EDGE=EDGE)

    # 2d Newton Raphson inversion and find new increments
    iter = 0
    while (((Ω₁ - f1)^2 + (Ω₂ - f2)^2) > eps^2)

        f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,d3ψ,aguess,eguess,da=da,de=de,TOLECC=TOLECC,verbose=verbose,NINT=NINT,EDGE=EDGE)

        jacobian = [df1da df1de ; df2da df2de]

        # this increment reports occasional failures; why?
        try
            increment = jacobian \ (-([f1 ; f2] - [Ω₁ ; Ω₂]))
            aguess,eguess = aguess + increment[1],eguess + increment[2]
        catch e # this catch appears to not work because LAPACK is doing something under the hood
            if verbose>0
                println("OrbitalElements.NumericalInversion.jl: bad division for Jacobian=$jacobian and (f1,f2)=($f1,$f2), (Ω₁,Ω₂)=($Ω₁,$Ω₂) at (a,e)=($aguess,$eguess).")
            end
            # are we just in some tiny bad patch? # reset to 'safe' values
            aguess,eguess = aguess + 10da,0.5
            increment = [0;0]

            if iter > ITERMAX
                finaltol = ((Ω₁ - f1)^2 + (Ω₂ - f2)^2)
                return aguess,eguess,-2,finaltol
            end
        end

        # the try...catch above is failing for some reason
        if (@isdefined increment) == false
            increment = [0;0]
        end

        # @WARNING: these appear to have broken something.
        # if bad guesses, needs to reset to a different part of space
        # can't go too small
        if eguess < TOLECC
            # go halfway between the previous guess and 0.
            try
                # reset eguess value
                eguess = eguess - increment[2]
                eguess = max(TOLECC,0.5eguess)
            catch e
                if verbose>0
                    println("OrbitalElements.NumericalInversion.jl: guessing close to ecc=0: ",eguess," (a=",aguess,")")
                end
                eguess = max(TOLECC,0.5eguess)
            end
        end

        if eguess >= (1-TOLECC)
            # go halfway between the previous guess and 1.
            try
                eguess = eguess - increment[2]
                eguess = min(1-TOLECC,eguess + 0.5*(1-eguess))
            catch e
                println("OrbitalElements.NumericalInversion.jl: guessing close to ecc=1: ",eguess," (a=",aguess,") for increment ",increment)
                eguess = min(1-TOLECC,eguess + 0.5*(1-eguess))
            end
        end

        iter += 1

        if iter > ITERMAX
            break
        end

    end


    if verbose > 0
        println("OrbitalElements.NumericalInversion.jl: niter=",iter)
    end

    finaltol = ((Ω₁ - f1)^2 + (Ω₂ - f2)^2)

    # check here to not allow bad values?
    if isnan(aguess) | isnan(eguess)
        if verbose>0
            println("OrbitalElements.NumericalInversion.jl: failed for inputs (Ω₁,Ω₂)=($Ω₁,$Ω₂).")
        end

        # return a semi-equivalent circular orbit, as the failure mode is mostly very small orbits
        return acirc,0.0,-1,finaltol

    else
        return aguess,eguess,iter,finaltol
    end
end




"""ae_from_omega1omega2_brute(Ω₁,Ω₂,ψ,dψ,d2ψ[,eps,maxiter])

basic Newton-Raphson algorithm to find (a,e) from (Ω₁,Ω₂) brute force derivatives.

@IMPROVE add escape for circular orbits
"""
function ae_from_omega1omega2_brute(Ω₁::Float64,Ω₂::Float64,
                                    ψ::Function,
                                    dψ::Function,
                                    d2ψ::Function,
                                    eps::Float64=1*10^(-10),
                                    maxiter::Int64=1000,
                                    TOLECC::Float64=0.001,TOLA::Float64=0.0001,
                                    da::Float64=0.0001,de::Float64=0.0001,
                                    verbose::Int64=0)
    #

    # get the circular orbit (maximum radius) for a given Ω₁,Ω₂. use the stronger constraint.
    acirc = Omega1circ_to_radius(Ω₁,dψ,d2ψ)


    # check to make sure we aren't very close to circular (radial isn't a problem)
    #f1circ,f2circ = ComputeFrequenciesAE(ψ,dψ,d2ψ,acirc,0.0)
    #if (((Ω₁ - f1circ)^2 + (Ω₂ - f2circ)^2) < eps^2)
    #    return acirc,0.0
    #end


    # then start from ecc=0.5 and take numerical derivatives
    aguess = acirc
    eguess = 0.5
    f1,f2 = ComputeFrequenciesAE(ψ,dψ,d2ψ,aguess,eguess)

    # 2d Newton Raphson inversion and find new increments
    iter = 0
    while (((Ω₁ - f1)^2 + (Ω₂ - f2)^2) > eps^2)

        f1,f2,df1da,df2da,df1de,df2de = ComputeFrequenciesAEWithDeriv(ψ,dψ,d2ψ,aguess,eguess,da=da,de=de,TOLECC=TOLECC,verbose=verbose)

        jacobian = [df1da df1de ; df2da df2de]

        # this increment reports occasional failures; why?
        try
            increment = jacobian \ (-([f1 ; f2] - [Ω₁ ; Ω₂]))
            aguess,eguess = aguess + increment[1],eguess + increment[2]
        catch e # this catch appears to not work because LAPACK is doing something under the hood
            if verbose>0
                println("OrbitalElements.NumericalInversion.jl: bad division for Jacobian=$jacobian and (f1,f2)=($f1,$f2), (Ω₁,Ω₂)=($Ω₁,$Ω₂).")
            end
            # are we just in some tiny bad patch? # reset to 'safe' values
            aguess,eguess = aguess + 1.e-6, 0.5
            increment = [0;0]
        end

        # the try...catch above is failing for some reason
        if (@isdefined increment) == false
            increment = [0;0]
        end




        # @WARNING: these appear to have broken something.
        # if bad guesses, needs to reset to a different part of space
        # can't go too small
        if eguess < TOLECC
            # go halfway between the previous guess and 0.
            try
                # reset eguess value
                eguess = eguess - increment[2]
                eguess = max(TOLECC,0.5eguess)
            catch e
                if verbose>0
                    println("OrbitalElements.NumericalInversion.jl: guessing close to ecc=0: ",eguess," (a=",aguess,")")
                end
                eguess = max(TOLECC,0.5eguess)
            end
        end

        if eguess >= (1-TOLECC)
            # go halfway between the previous guess and 1.
            try
                eguess = eguess - increment[2]
                eguess = min(1-TOLECC,eguess + 0.5*(1-eguess))
            catch e
                println("OrbitalElements.NumericalInversion.jl: guessing close to ecc=1: ",eguess," (a=",aguess,") for increment ",increment)
                eguess = min(1-TOLECC,eguess + 0.5*(1-eguess))
            end
        end

        if aguess < TOLA
            aguess = TOLA
            if verbose>1
                println("OrbitalElements.NumericalInversion.jl: guessing close to a=0: ",aguess," (e=",eguess,")")
            end
        end

        if (verbose>2)
            println("OrbitalElements/NumericalInversion.jl: iter=",iter," a=",aguess," e=",eguess)
        end

        iter += 1
        if iter > maxiter
            break
        end
    end

    if verbose > 0
        println("OrbitalElements/NumericalInversion.jl: niter=",iter)
    end

    finaltol = ((Ω₁ - f1)^2 + (Ω₂ - f2)^2)

    # check here to not allow bad values?
    if isnan(aguess) | isnan(eguess)
        if verbose>0
            println("OrbitalElements/NumericalInversion.jl: failed for inputs (Ω₁,Ω₂)=",Ω₁,Ω₂)
        end
        return acirc,0.5,-1,finaltol
    else
        return aguess,eguess,iter,finaltol
    end
end



"""ae_from_EL_brute(E,L,ψ,dψ,d2ψ[,eps,maxiter,TOLECC,verbose])
basic Newton-Raphson algorithm to find (a,e) from (E,L) brute force derivatives.
@IMPROVE add escape for circular orbits
"""
function ae_from_EL_brute(E::Float64,L::Float64,
                          ψ::Function,
                          dψ::Function,
                          d2ψ::Function,
                          eps::Float64=1*10^(-6),
                          maxiter::Int64=1000,
                          TOLECC::Float64=0.001,
                          verbose::Int64=0)
    #

    # get the circular orbit (maximum radius) for a given E. use the stronger constraint.
    #acirc = Omega1circ_to_radius(Ω₁,dψ,d2ψ)
    # is this the best launching eccentricity?
    aguess,eccguess = 1.,TOLECC

    rpguess,raguess = rpra_from_ae(aguess,eccguess)

    if (verbose>0)
      println("iter=",-1," aguess=",aguess," eguess=",eccguess)
    end

    Eguess,Lguess,dEda,dEde,dLda,dLde = dEdLFromRpRa(ψ,dψ,d2ψ,rpguess,raguess,da=0.0001,de=0.0001,TOLECC=TOLECC)



    # 2d Newton Raphson inversion and find new increments
    iter = 0
    while (((E - Eguess)^2 + (L - Lguess)^2) > eps^2)

        # convert to rp,ra for EL input
        rpguess,raguess = rpra_from_ae(aguess,eccguess)

        Eguess,Lguess,dEda,dEde,dLda,dLde = dEdLFromRpRa(ψ,dψ,d2ψ,rpguess,raguess,da=0.0001,de=0.0001,TOLECC=TOLECC)

        jacobian = [dEda dEde ; dLda dLde]
        increment = jacobian \ (-([Eguess;Lguess] - [E ; L]))

        aguess,eccguess = aguess + increment[1],eccguess + increment[2]

        # if bad guesses, needs to reset to a different part of space
        # can't go too small
        if eccguess < TOLECC
            # go halfway between the previous guess and 0.
            eccguess = eccguess - increment[2]
            eccguess = 0.5*eccguess
        end

        if eccguess >= 1.0-0.000001
            # go halfway between the previous guess and 1.
            eccguess = eccguess - increment[2]
            eccguess = eccguess + 0.5*(1-eccguess)
        end

        if aguess < 0.00000001
            aguess = 0.00000001
        end

        if (verbose>0)
            println("iter=",iter," aguess=",aguess," eguess=",eccguess)
        end

        iter += 1
        if iter > maxiter
            break
        end
    end

    return aguess,eccguess
end
