# Description of one XSPEC model parameter, carrying exactly the information
# that XSPEC's model.dat needs.

"""
    Parameter(name, default; unit = "", min, max, soft_min = min, soft_max = max,
              delta = 1% of the range, frozen = false, kind = :fit)

One fit parameter of an XSPEC model, as it will appear in `model.dat`.

- `name`: parameter name as shown in XSPEC (letters, digits and underscores).
- `default`: initial value.
- `min`, `max`: hard limits; XSPEC never evaluates the model outside them.
- `soft_min`, `soft_max`: soft limits (XSPEC's fit is cautious beyond them).
- `delta`: initial fit step. XSPEC freezes a parameter whose delta is negative,
  which is what `frozen = true` does.
- `kind`: `:fit` (ordinary), `:scale` (XSPEC `*` parameter: fixed scaling, not
  fitted) or `:switch` (XSPEC `\$` parameter: an integer-valued option).
"""
struct Parameter
    name::String
    unit::String
    default::Float64
    hard_min::Float64
    soft_min::Float64
    soft_max::Float64
    hard_max::Float64
    delta::Float64
    kind::Symbol
end

function Parameter(
    name::AbstractString,
    default::Real;
    unit::AbstractString = "",
    min::Real,
    max::Real,
    soft_min::Real = min,
    soft_max::Real = max,
    delta::Real = round((max - min) / 100; sigdigits = 2),
    frozen::Bool = false,
    kind::Symbol = :fit,
)
    occursin(r"^[A-Za-z][A-Za-z0-9_]*$", name) ||
        throw(ArgumentError("parameter name $(repr(name)) must be an identifier (letters, digits, underscores)"))
    min <= soft_min <= default <= soft_max <= max ||
        throw(ArgumentError("parameter $name needs min <= soft_min <= default <= soft_max <= max"))
    kind in (:fit, :scale, :switch) ||
        throw(ArgumentError("parameter kind must be :fit, :scale or :switch"))
    step = frozen ? -abs(delta) : abs(delta)
    return Parameter(String(name), String(unit), default, min, soft_min, soft_max, max, step, kind)
end

"""
    model_dat_line(parameter) -> String

The `model.dat` line describing this parameter:
`name unit default hard_min soft_min soft_max hard_max delta`.
Following XSPEC's own `model.dat`, scale parameters are written as
`*name unit value` and switch parameters as `\$name value`.
"""
function model_dat_line(p::Parameter)
    unit = isempty(p.unit) ? "\" \"" : p.unit
    if p.kind == :scale
        return "*$(p.name) $unit $(p.default)"
    elseif p.kind == :switch
        return "\$$(p.name) $(Int(round(p.default)))"
    end
    return join(string.([p.name, unit, p.default, p.hard_min, p.soft_min, p.soft_max, p.hard_max, p.delta]), " ")
end
