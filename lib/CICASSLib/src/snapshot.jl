# Reader for the HDF5-free `.cicass` raw dump written by cicass/makeCosICs/capi_out.c.
# This is the hand-off vehicle into Enzo (gas grid + DM particles) and RAMSES
# (grafic gas fields + DM particles), carrying the baryon–DM streaming offset.

"""
    CICASSSnapshot

Realized streaming ICs loaded from a `.cicass` dump. Header scalars plus:

- `dm_pos :: Matrix{T}` — N³×3 DM positions, **box fraction [0,1)**
- `dm_vel :: Matrix{T}` — N³×3 DM velocities, **physical peculiar km/s**
- `gas_delta :: Vector{T}` — N³ baryon overdensity δ_b on the regular grid
- `gas_vel :: Matrix{T}` — N³×3 gas velocities, physical peculiar km/s
- `gas_temp :: Vector{T}` — N³ temperature field; T[K] = (gas_temp+1)·`tavg`

Grid ordering is CICASS-native C order `idx = i + j*N + k*N²` (i fastest); use
[`grid3d`](@ref) to reshape a field vector to an `(N,N,N)` array.

The streaming signature is `mean(gas_vel - dm_vel)` ≈ `vbc·(1+z)/1001` km/s along
a single axis (zero for `vbc=0`).
"""
struct CICASSSnapshot{T<:AbstractFloat}
    n::Int
    nspecies::Int
    box::Float64        # Mpc/h
    zinit::Float64
    omega_m::Float64
    omega_b::Float64
    omega_l::Float64
    hconst::Float64
    m_dm::Float64       # 1e10 Msun/h
    m_gas::Float64
    vbc::Float64        # km/s @ z=1000
    tavg::Float64       # K
    dm_pos::Matrix{T}
    dm_vel::Matrix{T}
    gas_delta::Vector{T}
    gas_vel::Matrix{T}
    gas_temp::Vector{T}
end

"Reshape a length-N³ CICASS field vector to an `(N,N,N)` array (C order, i fastest)."
grid3d(snap::CICASSSnapshot, field::AbstractVector) =
    reshape(field, snap.n, snap.n, snap.n)

function _read_col(io, ::Type{T}, n::Integer) where {T<:AbstractFloat}
    v = Vector{T}(undef, n)
    read!(io, v)
    return v
end

function _read_snapshot_body(io, ::Type{T}, n::Integer, nsp::Integer, hd::AbstractVector{Float64}) where {T<:AbstractFloat}
    box, zinit, omm, omb, oml, h, mdm, mgas, vbc, tavg = hd
    N3 = n * n * n
    rdcol() = _read_col(io, T, N3)
    # DM particles: pos x,y,z then vel x,y,z (axis-contiguous)
    px, py, pz = rdcol(), rdcol(), rdcol()
    vx, vy, vz = rdcol(), rdcol(), rdcol()
    # gas grid: delta, vel x,y,z, temp
    gd = rdcol()
    gvx, gvy, gvz = rdcol(), rdcol(), rdcol()
    gt = rdcol()

    dm_pos  = hcat(px, py, pz)
    dm_vel  = hcat(vx, vy, vz)
    gas_vel = hcat(gvx, gvy, gvz)
    return CICASSSnapshot(n, nsp, box, zinit, omm, omb, oml, h, mdm, mgas, vbc, tavg,
                          dm_pos, dm_vel, gd, gas_vel, gt)
end

"""
    read_snapshot(path::AbstractString) -> CICASSSnapshot

Load a `.cicass` raw dump. `CICASS01` stores f64 fields; `CICASS02` stores f32
fields. The temporary `CICASSF4` magic is accepted for local pre-merge f32 files.
Header metadata remains f64 in all formats.
"""
function read_snapshot(path::AbstractString)
    open(path, "r") do io
        magic_s = String(read(io, 8))
        field_type = magic_s == "CICASS01" ? Float64 :
                     (magic_s == "CICASS02" || magic_s == "CICASSF4") ? Float32 :
                     error("not a CICASS snapshot (bad magic $(repr(magic_s))): $path")
        n   = Int(read(io, Int32))
        nsp = Int(read(io, Int32))
        hd  = Vector{Float64}(undef, 10)
        read!(io, hd)
        return _read_snapshot_body(io, field_type, n, nsp, hd)
    end
end

"""
    streaming_velocity(snap) -> NTuple{3,Float64}

Mean `gas_vel - dm_vel` (physical peculiar km/s) — the coherent baryon–DM bulk
streaming offset realized in the ICs. ~`vbc·(1+zinit)/1001` along one axis.
"""
function streaming_velocity(snap::CICASSSnapshot)
    n = size(snap.dm_vel, 1)
    Tuple((sum(Float64, @view snap.gas_vel[:, d]) - sum(Float64, @view snap.dm_vel[:, d])) / n for d in 1:3)
end
