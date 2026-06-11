# Best-effort native build for CICASSLib: runs the CICASS C-ABI build script on a
# sibling `cicass` checkout. Mirrors MusicLib/deps/build.jl — never hard-fails the
# whole Pkg.build (so the package installs even without the checkout / a Mac
# toolchain); CICASSLib.available() gates every live call at runtime.

const CICASSLIB_SRC = @__DIR__
# CICASS.jl/lib/CICASSLib/deps → up 4 → Projects/, then cicass/
const CICASS = normpath(joinpath(CICASSLIB_SRC, "..", "..", "..", "..", "cicass"))
const SCRIPT = joinpath(CICASS, "deps", "build_cicass_darwin.sh")
const LIB = joinpath(CICASS, "build", Sys.isapple() ? "libcicass_capi.dylib" : "libcicass_capi.so")

if haskey(ENV, "CICASS_LIB")
    @info "CICASSLib: CICASS_LIB is set, skipping native build" lib = ENV["CICASS_LIB"]
elseif !isfile(SCRIPT)
    @warn "CICASSLib: no CICASS checkout / build script found; skipping native build. " *
          "Clone https://github.com/astromcquinn/CICASS to $CICASS, or set ENV[\"CICASS_LIB\"]." script = SCRIPT
elseif !Sys.isapple()
    @warn "CICASSLib: build_cicass_darwin.sh targets macOS; build libcicass_capi manually on this platform." script = SCRIPT
else
    try
        run(`bash $SCRIPT $CICASS`)
        isfile(LIB) || error("build script ran but $LIB is missing")
        @info "CICASSLib: built CICASS C-ABI library" lib = LIB
    catch err
        @warn "CICASSLib: native CICASS build failed; set ENV[\"CICASS_LIB\"] to a prebuilt library." exception = err
    end
end
