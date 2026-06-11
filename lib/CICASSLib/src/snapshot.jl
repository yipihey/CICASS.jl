# Reader for the HDF5-free `.cicass` raw dump written by cicass/makeCosICs/capi_out.c.
# This is the hand-off vehicle into Enzo (gas grid + DM particles) and RAMSES
# (grafic gas fields + DM particles), carrying the baryon–DM streaming offset.

"""
    CICASSSnapshot

Realized streaming ICs loaded from a `.cicass` dump. Header scalars plus:

- `dm_pos :: Matrix{Float64}` — N³×3 DM positions, **box fraction [0,1)**
- `dm_vel :: Matrix{Float64}` — N³×3 DM velocities, **physical peculiar km/s**
- `gas_delta :: Vector{Float64}` — N³ baryon overdensity δ_b on the regular grid
- `gas_vel :: Matrix{Float64}` — N³×3 gas velocities, physical peculiar km/s
- `gas_temp :: Vector{Float64}` — N³ temperature field; T[K] = (gas_temp+1)·`tavg`

Grid ordering is CICASS-native C order `idx = i + j*N + k*N²` (i fastest); use
[`grid3d`](@ref) to reshape a field vector to an `(N,N,N)` array.

The streaming signature is `mean(gas_vel - dm_vel)` ≈ `vbc·(1+z)/1001` km/s along
a single axis (zero for `vbc=0`).
"""
struct CICASSSnapshot
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
    dm_pos::Matrix{Float64}
    dm_vel::Matrix{Float64}
    gas_delta::Vector{Float64}
    gas_vel::Matrix{Float64}
    gas_temp::Vector{Float64}
end

"Reshape a length-N³ CICASS field vector to an `(N,N,N)` array (C order, i fastest)."
grid3d(snap::CICASSSnapshot, field::AbstractVector) =
    reshape(field, snap.n, snap.n, snap.n)

"""
    read_snapshot(path::AbstractString) -> CICASSSnapshot

Load a `.cicass` raw dump (little-endian f64; see `capi_out.c` for the layout).
"""
function read_snapshot(path::AbstractString)
    open(path, "r") do io
        magic = read(io, 8)
        String(magic) == "CICASS01" ||
            error("not a CICASS snapshot (bad magic $(repr(String(magic)))): $path")
        n   = Int(read(io, Int32))
        nsp = Int(read(io, Int32))
        hd  = Vector{Float64}(undef, 10)
        read!(io, hd)
        box, zinit, omm, omb, oml, h, mdm, mgas, vbc, tavg = hd
        N3 = n * n * n

        rdcol() = (v = Vector{Float64}(undef, N3); read!(io, v); v)
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
end

"""
    streaming_velocity(snap) -> NTuple{3,Float64}

Mean `gas_vel - dm_vel` (physical peculiar km/s) — the coherent baryon–DM bulk
streaming offset realized in the ICs. ~`vbc·(1+zinit)/1001` along one axis.
"""
function streaming_velocity(snap::CICASSSnapshot)
    n = size(snap.dm_vel, 1)
    Tuple((sum(@view snap.gas_vel[:, d]) - sum(@view snap.dm_vel[:, d])) / n for d in 1:3)
end
