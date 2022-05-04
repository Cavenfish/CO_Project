module MvH_Potential

include("MorsePotential.jl")
include("DispersionPotential.jl")
include("CoulombPotential.jl")
include("ExchangePotential.jl")

import .MorsePotential.MorseParameters
import .DispersionPotential.DispersionParameters
import .CoulombPotential.CoulombParameters
import .ExchangePotential.ExchangeParameters

export MorseParameters
export DispersionParameters
export CoulombParameters
export ExchangeParameters

end #module
