# Initial thermodynamic state from CICASS's RECFAST recombination table.
# For a cosmological simulation started at high redshift, the gas temperature and
# residual electron fraction must be set self-consistently (the gas is Compton-
# coupled to the CMB down to z ~ 150, then cools adiabatically). CICASS ships the
# RECFAST output `vbc_transfer/recfast/xeTrecfast.out` (columns: z, x_e, T_gas[K];
# z = 1630 -> 0); this exposes it so any starting redshift maps to an initial
# temperature and electron density for the MultiSpecies / Grackle fields.

const _RECFAST_CACHE = Ref{Union{Nothing,NTuple{3,Vector{Float64}}}}(nothing)

"Path to CICASS's RECFAST table (override with ENV[\"CICASS_RECFAST\"])."
recfast_path() = get(ENV, "CICASS_RECFAST",
    normpath(joinpath(cicass_root(), "vbc_transfer", "recfast", "xeTrecfast.out")))

# Read (z, x_e, T_gas) once; table is z-descending, first line is the row count.
function _recfast_table()
    _RECFAST_CACHE[] === nothing || return _RECFAST_CACHE[]
    p = recfast_path()
    isfile(p) || error("RECFAST table not found at $p (set ENV[\"CICASS_RECFAST\"]).")
    z = Float64[]; xe = Float64[]; tg = Float64[]
    for (i, line) in enumerate(eachline(p))
        i == 1 && continue                    # header: number of rows
        t = split(line)
        length(t) >= 3 || continue
        push!(z, parse(Float64, t[1]))
        push!(xe, parse(Float64, t[2]))
        push!(tg, parse(Float64, t[3]))
    end
    # store ascending in z for searchsortedfirst
    perm = sortperm(z)
    tab = (z[perm], xe[perm], tg[perm])
    _RECFAST_CACHE[] = tab
    return tab
end

_lininterp(x, xs, ys) = begin
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    j = searchsortedfirst(xs, x)
    w = (x - xs[j-1]) / (xs[j] - xs[j-1])
    ys[j-1] * (1 - w) + ys[j] * w
end

"""
    thermal_state(z) -> (; T_gas, x_e, T_cmb)

Initial gas temperature `T_gas` [K] and free-electron fraction
`x_e = n_e / n_H` at redshift `z`, interpolated from CICASS's RECFAST table
(valid 0 ≤ z ≤ 1630), plus the CMB temperature `T_cmb = 2.73·(1+z)`. At z ≈ 1000
`x_e ≈ 0.047`, `T_gas ≈ T_cmb`; the gas thermally decouples from the CMB near
z ≈ 150 and then `T_gas < T_cmb`.
"""
function thermal_state(z::Real)
    zs, xe, tg = _recfast_table()
    (0.0 <= z <= zs[end]) || @warn "thermal_state: z=$z outside RECFAST table [0, $(zs[end])]; extrapolating"
    return (; T_gas = _lininterp(Float64(z), zs, tg),
            x_e   = _lininterp(Float64(z), zs, xe),
            T_cmb = 2.73 * (1 + z))
end

"""
    multispecies_fractions(z; X_H=0.76) -> NamedTuple

Map CICASS's RECFAST state at redshift `z` to primordial MultiSpecies / Grackle
*number fractions by mass* (relative to total gas density), ready to load into
the species fields for a high-z start. Hydrogen carries the ionization
(`x_HII = x_e`); helium is taken neutral (it recombines at z > 1630):

  HI = (1−x_e)·X_H,  HII = x_e·X_H,  e = x_e·X_H·(m_e/m_H≈1 in Grackle's n_e·m_H),
  HeI = 1−X_H,  HeII = HeIII = H2I = H2II = HM ≈ 0.

Returns the fractions plus `T_gas` for setting the internal energy.
"""
function multispecies_fractions(z::Real; X_H::Real = 0.76)
    s = thermal_state(z)
    xe = s.x_e
    return (; HI = (1 - xe) * X_H, HII = xe * X_H, e = xe * X_H,
            HeI = 1 - X_H, HeII = 0.0, HeIII = 0.0,
            HM = 0.0, H2I = 0.0, H2II = 0.0,
            T_gas = s.T_gas, x_e = xe, z = Float64(z))
end
