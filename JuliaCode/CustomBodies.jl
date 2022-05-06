module CustomBodies

using Unitful
using UnitfulAtomic
using StaticArrays

Unitful.register(@__MODULE__)

@unit Ang "Ang" Angstrom 0.1u"nm" true

export Atom

abstract type Particle 
end

DiffEqBase.@def pos_vel begin
  r::SVector{3,rType}
  v::SVector{3,vType}
end

struct Atom{rType <: typeof(1.0u"Ang"),
	    vType <: typeof(1.0u"(eV/u)^0.5"),
	    mType <: typeof(1.0u"u"),
	    qType <: typeof(1.0u"e_au"),
	    sType <: AbstractString} <: Particle

  @pos_vel
  m::mType
  q::qType
  s::sType
end

end 
