module CustomSystem

using Unitful
using UnitfulAtomic
using StaticArrays

export AtomSystem

abstract type System
end

struct AtomSystem{bType <: Particle,
		  EType <: typeof(1.0u"eV")} <: System
  
  atoms::Vector{bType}
  KE::EType
  PE::EType
  TE::EType
end

end
