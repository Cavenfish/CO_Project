module CustomBodies

using Unitful
using UnitfulAtomic
using StaticArrays
using DiffEqBase
using NBodySimulator

Unitful.register(@__MODULE__)

@unit Ang "Ang" Angstrom 0.1u"nm" true

export Atom

DiffEqBase.@def pos_vel begin
  r::SVector{3,rType}
  v::SVector{3,vType}
end

struct Atom{rType <: Float64,
	    vType <: Float64,
	    mType <: Float64,
	    qType <: Float64,
	    bType <: UInt64,
	    sType <: AbstractString} <: NBodySimulator.Body

  @pos_vel
  m::mType
  q::qType
  b::bType
  s::sType
end

end 
