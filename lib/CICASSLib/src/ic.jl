# Streaming-velocity cosmological IC "problem setup": a typed CICASSSpec →
# transfer.x (2D v_bc TF grid) → in-process cicass_generate → a `.cicass` raw dump
# ({DM particles + gas grid + bulk streaming offset}). This is the Julia-side
# surface the MultiCode framework calls to make streaming ICs for Enzo / RAMSES.

"""
    CICASSSpec

A streaming-velocity initial-conditions problem for CICASS, captured as a typed
Julia value (the streaming-aware analogue of `MusicSpec`).

The cosmology *shape* (σ₈, nₛ, the transfer function) is baked into the shipped
CAMB transfer functions in `cicass/vbc_transfer/TFs` — the default cosmology is
flat ΛCDM with h=0.71, Ωₘ=0.27, Ω_b=0.046, σ₈=0.8, n=0.95. Changing `Omega_m`,
`Omega_b` or `hconst` only rescales the realization; a *different* power-spectrum
shape requires regenerating the CAMB inputs.

Fields (all keyword):
- `boxlength` — comoving box size [Mpc/h]
- `zstart` — initial redshift (CICASS evolves the linear equations to here; ≳100 ok)
- `ngrid` — cells per axis N (gas grid is N³; DM is N³ particles)
- `glass_dim` — tiled glass-file dimension G (must divide `ngrid`)
- `vbc` — streaming velocity [km/s] at z=1000 (0 ⇒ no streaming, isotropic 1D TF)
- `Omega_m, Omega_b, hconst` — realization cosmology (`-O`, `-B`, `-H`)
- `species` — 1 (DM only) or 2 (DM + gas)
- `seed` — random-number seed (`-R`)
- `tf_mode` — transfer.x `-D` mode (1 = direct CAMB z=0,99,100,101 with v_bc applied;
  3 = full evolution from z≈1000)
- `tf_base` — base name of the CAMB transfer files (`-S`; shipped = `initSB_transfer_out`)
- `glass_file` — absolute path to the glass file (`-g`)
- `filename` — base output name (`<filename>.cicass`)
- `real_bytes` — field storage bytes in the `.cicass` snapshot: 8 (`CICASS01`) or 4 (`CICASS02`)
"""
function _default_real_bytes()
    n = tryparse(Int, get(ENV, "CICASS_REAL_BYTES", "4"))
    return n == 8 ? 8 : 4
end

Base.@kwdef struct CICASSSpec
    boxlength::Float64
    zstart::Float64    = 100.0
    ngrid::Int         = 128
    glass_dim::Int     = 128
    vbc::Float64       = 0.0
    Omega_m::Float64   = 0.27
    Omega_b::Float64   = 0.046
    hconst::Float64    = 0.71
    species::Int       = 2
    seed::Int          = 113334
    # Angulo & Pontzen (2016) variance-suppressed ICs. `fix_amplitude`: set every Fourier mode's
    # amplitude to exactly sqrt(P(k)) (random phases only) so the realized P(k) matches the input
    # mode-by-mode — kills box-to-box amplitude scatter. `flip_phase`: invert all phases for the
    # PAIRED run (average the fix_amplitude pair {flip_phase=false,true} to cancel leading
    # non-Gaussian variance). Both default off ⇒ standard Rayleigh realization.
    fix_amplitude::Bool = false
    flip_phase::Bool    = false
    tf_mode::Int       = 1
    tf_base::String    = "initSB_transfer_out"
    glass_file::String = ""
    filename::String   = "cicass_ics"
    real_bytes::Int    = _default_real_bytes()
end

"Default glass-file path from the sibling checkout."
_default_glass() = normpath(joinpath(cicass_root(), "glass", "glass_128_usethis"))

"CICASS's TF-grid / output filename convention: `initSimCartZI<Z>.1f_Vbc<V>.1f_<N>_<box>.1f.dat`."
_tf_gridname(spec::CICASSSpec) =
    string("initSimCartZI", _f1(spec.zstart), "_Vbc", _f1(spec.vbc),
           "_", spec.ngrid, "_", _f1(spec.boxlength), ".dat")

# CICASS uses C `%.1lf` formatting in the filename; reproduce it exactly.
_f1(x::Real) = string(round(Float64(x); digits = 1))

function _ic_threads()
    for key in ("CICASS_FFT_THREADS", "FFTW_NUM_THREADS", "OMP_NUM_THREADS")
        if haskey(ENV, key)
            n = tryparse(Int, ENV[key])
            n !== nothing && n > 0 && return n
        end
    end
    return Sys.CPU_THREADS
end

function _with_ic_thread_env(f; real_bytes::Union{Nothing,Integer} = nothing)
    n = string(_ic_threads())
    old = Dict(k => get(ENV, k, nothing) for k in ("CICASS_FFT_THREADS", "OMP_NUM_THREADS", "CICASS_REAL_BYTES"))
    ENV["CICASS_FFT_THREADS"] = n
    haskey(ENV, "OMP_NUM_THREADS") || (ENV["OMP_NUM_THREADS"] = n)
    real_bytes === nothing || (ENV["CICASS_REAL_BYTES"] = string(real_bytes == 8 ? 8 : 4))
    try
        return f()
    finally
        for (k, v) in old
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
    end
end

"""
    make_tf(spec::CICASSSpec; outdir=mktempdir(), force=false) -> String

Run the `transfer.x` helper to produce the 2D (k⊥, k∥) transfer-function grid that
the realizer consumes for `spec`, returning the absolute path to the grid file.
`transfer.x` runs with CWD = `cicass/vbc_transfer` (it reads `TFs/` and `recfast/`
relative), writing to `vbc_transfer/IC_outputs/`; the grid is then copied to
`outdir` under its canonical name. With `force=false` an existing grid is reused.
"""
function make_tf(spec::CICASSSpec; outdir::AbstractString = mktempdir(), force::Bool = false)
    mkpath(outdir)
    dest = joinpath(outdir, _tf_gridname(spec))
    (!force && isfile(dest)) && return dest

    tx = transfer_path()
    isfile(tx) || error("transfer.x not found at $tx — build with cicass/deps/build_cicass_darwin.sh " *
                        "or set ENV[\"CICASS_TRANSFER\"].")
    vbcdir = normpath(joinpath(cicass_root(), "vbc_transfer"))
    isdir(joinpath(vbcdir, "TFs")) ||
        error("CICASS TFs/ dir missing under $vbcdir (CAMB transfer inputs).")

    args = ["-B$(spec.boxlength)", "-N$(spec.ngrid)", "-V$(spec.vbc)",
            "-Z$(round(Int, spec.zstart))", "-D$(spec.tf_mode)", "-S$(spec.tf_base)"]
    src = joinpath(vbcdir, "IC_outputs", _tf_gridname(spec))
    isfile(src) && rm(src; force = true)            # don't mistake a stale grid for success
    # NB transfer.x returns a non-zero exit code even on success — gate on the
    # output file, not the status.
    cd(vbcdir) do
        _with_ic_thread_env() do
            run(pipeline(ignorestatus(`$tx $args`); stdout = devnull))
        end
    end
    isfile(src) || error("transfer.x did not produce $src (args: $(join(args, ' ')))")
    cp(src, dest; force = true)
    return dest
end

"""
    generate(spec::CICASSSpec; workdir=mktempdir(), tf_file=nothing, check=true) -> NamedTuple

Generate the streaming ICs for `spec`. Ensures the v_bc transfer grid is present
(via [`make_tf`](@ref) unless `tf_file` is supplied), stages it into `workdir`
under the name the realizer expects, runs [`cicass_generate`](@ref) in-process
(CWD = `workdir`), and returns `(; rc, dir, output, spec, tf)` where `output` is
the absolute path of the `.cicass` raw dump (load it with [`read_snapshot`](@ref)).
With `check=true`, a non-zero return raises with [`last_error`](@ref).
"""
function generate(spec::CICASSSpec; workdir::AbstractString = mktempdir(),
                  tf_file::Union{Nothing,AbstractString} = nothing, check::Bool = true)
    spec.real_bytes in (4, 8) || error("CICASSSpec.real_bytes must be 4 or 8 (got $(spec.real_bytes))")
    mkpath(workdir)
    glass = isempty(spec.glass_file) ? _default_glass() : spec.glass_file
    isfile(glass) || error("glass file not found: $glass")

    # 1) transfer grid in cwd under its canonical name
    tf = tf_file === nothing ? make_tf(spec; outdir = workdir) : String(tf_file)
    staged = joinpath(workdir, _tf_gridname(spec))
    (abspath(tf) == abspath(staged)) || cp(tf, staged; force = true)

    # 2) realize
    args = join(["-L$(spec.boxlength)", "-V$(spec.vbc)", "-N$(spec.ngrid)",
                 "-G$(spec.glass_dim)", "-Z$(round(Int, spec.zstart))",
                 "-H$(spec.hconst)", "-O$(spec.Omega_m)", "-B$(spec.Omega_b)",
                 "-S$(spec.species)", "-R$(spec.seed)",
                 "-F$(spec.fix_amplitude ? 1 : 0)", "-I$(spec.flip_phase ? 1 : 0)",
                 "-g$(glass)", "-o.", "-b$(spec.filename)"], " ")
    rc = cd(workdir) do
        _with_ic_thread_env(; real_bytes=spec.real_bytes) do
            cicass_generate(args)
        end
    end
    if check && rc != 0
        error("CICASS failed (rc=$rc) for vbc=$(spec.vbc): ", last_error())
    end
    return (; rc, dir = workdir, output = joinpath(workdir, spec.filename * ".cicass"),
            spec = spec, tf = staged)
end
