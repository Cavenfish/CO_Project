using NBodySimulator
import NBodySimulator.get_accelerating_function

struct MorseParameters <: PotentialParameters
  De = 11.23
  ρ0 = 2.3281*1.1282
  r0 = 1.1282
end

function get_accelerating_function(p::MorseParameters, sim::NBodySimulation)
  (dv, u, v, t, i) -> begin
    
  end #Morse acceleration
end
