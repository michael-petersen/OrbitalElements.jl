module OrbitalElements

# for zero-finding
using Optim
using Roots

# for interpolations
using Interpolations

# recommendations for various default parameters:
#   you can call these externally
const DEFAULT_TOLECC = 0.001
const DEFAULT_TOLA   = 0.001
const DEFAULT_EDGE   = 0.01
const DEFAULT_NINT   = 32
const DEFAULT_RMIN   = 1.e-6
const DEFAULT_RMAX   = 1.e6
const DEFAULT_TOL    = 1.e-12
const DEFAULT_AMIN   = 0.0
const DEFAULT_AMAX   = 1.0
const DEFAULT_DA     = 1.e-4
const DEFAULT_DE     = 1.e-4
const DEFAULT_ITERMAX= 100
const DEFAULT_O0     = 1.0

# utils: simple function extremisation, input/output, basic orbit coordinate transformations, simple integrators
include("Utils/ParameterStructure.jl")
include("Utils/Extremise.jl")
include("Utils/IO.jl")
include("Utils/OrbitDefinitions.jl")
include("Utils/Integrators.jl")
include("Utils/DerivationInterpolation.jl")

# bring in the test potentials (not strictly needed)
include("Potentials/isochrone.jl")
include("Potentials/plummer.jl")
include("Potentials/mestelzang.jl")
include("Potentials/kuzmintoomre.jl")

# bring in the test distribution functions (not strictly needed)
include("DistributionFunctions/isochrone.jl")
include("DistributionFunctions/plummer.jl")
include("DistributionFunctions/mestelzang.jl")
include("DistributionFunctions/miyamoto.jl")
include("DistributionFunctions/isochrone_discs.jl")

# enable energy and angular momentum computation (including expansions)
include("Utils/ComputeEL.jl")

# bring in the circular orbit frequencies
include("Circular/CircularFrequencies.jl")

# bring in the resonance mappings
include("Resonance/UVbounds.jl")
include("Resonance/ABtoUV.jl")

# the main wrapper for frequency and action calculations
# first, set the u (\in[-1,1]) integration limit to switch to expansions
const ULIMIT = 0.01
include("Frequencies.jl")

end # module
