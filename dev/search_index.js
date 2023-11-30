var documenterSearchIndex = {"docs":
[{"location":"functions.html#Systems","page":"Functions","title":"Systems","text":"","category":"section"},{"location":"functions.html#Potentials","page":"Functions","title":"Potentials","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.dψMestel","category":"page"},{"location":"functions.html#OrbitalElements.dψMestel","page":"Functions","title":"OrbitalElements.dψMestel","text":"dψMestel(r[, R0, V0, epsilon])\n\nthe Mestel potential derivative.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#Distribution-functions","page":"Functions","title":"Distribution functions","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.Miyamoto_DF\nOrbitalElements.IsoDiscKal_DF\nOrbitalElements.IsoDiscPLB_ndDFdJ\nOrbitalElements.ZangOuterTaperdL","category":"page"},{"location":"functions.html#OrbitalElements.Miyamoto_DF","page":"Functions","title":"OrbitalElements.Miyamoto_DF","text":"Miyamoto_DF(E, L)\n\nMiyamoto distribution function for Kuzmin-Toomre disc.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.IsoDiscKal_DF","page":"Functions","title":"OrbitalElements.IsoDiscKal_DF","text":"IsoDiscKal_DF(E, L)\n\nKalnajs distribution function for isochrone disc.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.IsoDiscPLB_ndDFdJ","page":"Functions","title":"OrbitalElements.IsoDiscPLB_ndDFdJ","text":"IsoDiscPLB_ndDFdJ(n1, n2, E, L, ndotOmega)\n\nKalnajs DF derivative w.r.t. the actions J.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.ZangOuterTaperdL","page":"Functions","title":"OrbitalElements.ZangOuterTaperdL","text":"Zang_outer_tapering_dL(L[, Rout, V0, mu])\n\nZang outer tapering derivative.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#Changes-of-coordinates","page":"Functions","title":"Changes of coordinates","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.Ω1circular\nOrbitalElements.EFromAE\nOrbitalElements.ELFromAE\nOrbitalElements.ComputeActionsAE\nOrbitalElements.FrequenciesFromαβ\nOrbitalElements.Getϖ\nOrbitalElements.FindVbound\nOrbitalElements.RcircFromΩ2circ","category":"page"},{"location":"functions.html#OrbitalElements.Ω1circular","page":"Functions","title":"OrbitalElements.Ω1circular","text":"Ω1circular(dψ,d2ψ,a)\n\nradial frequency for circular orbits, from the epicyclic approximation a is the semi-major axis (equivalent to r for a circular orbit)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.EFromAE","page":"Functions","title":"OrbitalElements.EFromAE","text":"EFromAE(ψ,dψ,a,e,params)\n\nenergy as a function of (a,e) for a given potential ψ (and its derivatives)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.ELFromAE","page":"Functions","title":"OrbitalElements.ELFromAE","text":"ELFromAE(ψ,dψ,a,e,params)\n\ncombined energy + angular momentum as a function of (a,e) for a given potenial ψ (and its derivatives)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.ComputeActionsAE","page":"Functions","title":"OrbitalElements.ComputeActionsAE","text":"ComputeActionsAE(ψ,dψ,a,e,params)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.FrequenciesFromαβ","page":"Functions","title":"OrbitalElements.FrequenciesFromαβ","text":"FrequenciesFromαβ(α,β,Ω₀)\n\nconverts frequencies ratios to frequencies\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.Getϖ","page":"Functions","title":"OrbitalElements.Getϖ","text":"Getϖ(ω,n1,n2,dψ,d2ψ,params)\n\ntranslate a complex frequency into a rescaled frequency. maps omega to -11\n\nFouvry & Prunet B3\n\n@ASSUMPTION:     - ω is dimensionless, that is, rescaled by Ω₀ already.\n\n\n\n\n\nGetϖ(ω,ωmin,ωmax)\n\nϖ version with ωmin, ωmax\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.FindVbound","page":"Functions","title":"OrbitalElements.FindVbound","text":"FindVbound(n1,n2,dψ,d2ψ,Ω₀,rmin,rmax)\n\nfind any valid non- 0 or 1 v value at u=-1 or u=1\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.RcircFromΩ2circ","page":"Functions","title":"OrbitalElements.RcircFromΩ2circ","text":"RcircFromΩ2circ(Ω₂,dψ,d2ψ,rmin,rmax[,tolx,tolf])\n\nperform backwards mapping from Omega_2 for a circular orbit to radius\n\n@ASSUMPTIONS:     - Ω2circular is a decreasing function of radius     - d2ψ used for value at 0.\n\n\n\n\n\n","category":"function"},{"location":"functions.html#Isochrone-specific-functions","page":"Functions","title":"Isochrone specific functions","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.IsochroneβAE","category":"page"},{"location":"functions.html#OrbitalElements.IsochroneβAE","page":"Functions","title":"OrbitalElements.IsochroneβAE","text":"compute the dimensionless function for Omega2 from (a,e) (Fouvry 21 eq. G7)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#Plummer-specific-functions","page":"Functions","title":"Plummer specific functions","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.dΘRpRaPlummer\nOrbitalElements.PlummerELFromSpSa","category":"page"},{"location":"functions.html#OrbitalElements.dΘRpRaPlummer","page":"Functions","title":"OrbitalElements.dΘRpRaPlummer","text":"the wrapped Theta derivative function\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.PlummerELFromSpSa","page":"Functions","title":"OrbitalElements.PlummerELFromSpSa","text":"translate From (sp,sa) to (E,L)\n\n\n\n\n\n","category":"function"},{"location":"functions.html#Utility-functions","page":"Functions","title":"Utility functions","text":"","category":"section"},{"location":"functions.html","page":"Functions","title":"Functions","text":"OrbitalElements.NonAsciiHandle\nOrbitalElements.GetResLinesJL","category":"page"},{"location":"functions.html#OrbitalElements.NonAsciiHandle","page":"Functions","title":"OrbitalElements.NonAsciiHandle","text":"NonAsciiHandle(x)\n\nconvert some extra unicode characters to ascii\n\n\n\n\n\n","category":"function"},{"location":"functions.html#OrbitalElements.GetResLinesJL","page":"Functions","title":"OrbitalElements.GetResLinesJL","text":"GetResLinesJL(ω,ψ,dψ,d2ψ,Kv,tabresonances)\n\n\n\n\n\n","category":"function"},{"location":"index.html#OrbitalElements.jl","page":"Home","title":"OrbitalElements.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Galactic dynamics orbits in Julia","category":"page"}]
}
