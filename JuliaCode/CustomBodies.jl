module CustomBodies

using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

Unitful.register(@__MODULE__)

@unit Ang "Ang" Angstrom 0.1u"nm" true

export Nucleus

DiffEqBase.@def pos_vel begin
  r::SVector{3,rType}
  v::SVector{3,vType}
end

struct Nucleus{rType <: typeof(1.0u"Ang"),
	       vType <: typeof(1.0u"(eV/u)^0.5"),
	       mType <: typeof(1.0u"u"),
	       qType <: typeof(1.0u"e_au"),
	       sType <: AbstractString} <: NBodySimulator.Body

  @pos_vel
  m::mType
  q::qType
  s::sType
end

end 
