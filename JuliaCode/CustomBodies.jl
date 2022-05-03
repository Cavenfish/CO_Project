module CustomBodies

using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic


struct Nucleus{rType <: Unitful.Length,
	       vType <: typeof(1.0u"(eV*u)^0.5"),
	       mType <: typeof(1.0u"u"),
	       qType <: typeof(1.0u"e_au"),
	       sType <: SubString}

  r::rType
  v::vType
  m::mType
  q::qType
  s::sType
end
