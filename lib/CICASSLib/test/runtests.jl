# CICASSLib tests. Pure-Julia spec/format tests always run; generation tests
# need the built libcicass_capi + transfer.x (gated on CICASSLib.available()).
# Run:  <julia> --project=test test/runtests.jl
#
# The v_bc=30 realization is ~60s (the 2D-TF interpolation); set
# ENV["CICASS_SKIP_STREAMING"]="1" to skip it.

using CICASSLib, Test
const CL = CICASSLib

@testset "CICASSLib" begin

    # ── spec + filename convention (no library needed) ──────────────────────
    @testset "spec / TF-grid name" begin
        spec = CICASSSpec(boxlength = 0.2, zstart = 100.0, ngrid = 128, vbc = 30.0)
        @test CL._tf_gridname(spec) == "initSimCartZI100.0_Vbc30.0_128_0.2.dat"
        @test spec.real_bytes in (4, 8)
        @test CICASSSpec(boxlength = 0.2, real_bytes = 4).real_bytes == 4
        s0 = CICASSSpec(boxlength = 1.0, zstart = 200.0, ngrid = 256, vbc = 0.0)
        @test CL._tf_gridname(s0) == "initSimCartZI200.0_Vbc0.0_256_1.0.dat"
    end

    @testset "CICASSF4 snapshot reader" begin
        path = tempname() * ".cicass"
        n = 2; n3 = n^3
        open(path, "w") do io
            write(io, "CICASSF4")
            write(io, Int32(n), Int32(2))
            write(io, Float64.([0.2, 100.0, 0.27, 0.046, 0.73, 0.71, 1.0, 0.1, 30.0, 100.0]))
            cols = [fill(Float32(0.1 * d), n3) for d in 1:11]
            cols[4] .= 0f0; cols[5] .= 0f0; cols[6] .= 0f0
            cols[8] .= 1f0; cols[9] .= 2f0; cols[10] .= 3f0
            foreach(c -> write(io, c), cols)
        end
        snap = read_snapshot(path)
        @test snap.n == n
        @test eltype(snap.dm_pos) === Float32
        @test eltype(snap.gas_delta) === Float32
        @test CL.streaming_velocity(snap) == (1.0, 2.0, 3.0)
    end

    # ── RECFAST thermal-state provider (no library needed, just the table) ──
    @testset "thermal_state / multispecies" begin
        if isfile(CICASSLib.recfast_path())
            s = thermal_state(1000.0)
            @test isapprox(s.x_e, 0.047; atol = 0.01)        # residual ionization at z=1000
            @test isapprox(s.T_gas, 2728.0; rtol = 0.05)     # gas ≈ T_CMB at z=1000
            @test isapprox(s.T_cmb, 2.73 * 1001; rtol = 1e-6)
            @test s.T_gas < s.T_cmb || isapprox(s.T_gas, s.T_cmb; rtol = 0.05)
            # decoupled at low z: gas colder than CMB
            @test thermal_state(50.0).T_gas < thermal_state(50.0).T_cmb
            mf = multispecies_fractions(1000.0)
            @test isapprox(mf.HI + mf.HII, 0.76; rtol = 1e-6) # hydrogen mass fraction
            @test isapprox(mf.HII, mf.e; rtol = 1e-9)         # charge balance (H-dominated)
            @test mf.H2I == 0.0
        else
            @info "CICASS RECFAST table not found; skipping thermal_state tests"
        end
    end

    if !CICASSLib.available()
        @info "CICASSLib: library not built; skipping generation tests" lib = CICASSLib.libpath()
    else
        @testset "version / precision" begin
            @test occursin("CICASS", CICASSLib.version())
            @test CICASSLib.real_bytes() == 8
        end

        # ── v_bc = 0: two-component IC, zero streaming offset ────────────────
        @testset "generate v_bc=0" begin
            spec = CICASSSpec(boxlength = 0.2, zstart = 100.0, ngrid = 128, vbc = 0.0,
                              filename = "test_v0")
            res = generate(spec)
            @test res.rc == 0
            @test isfile(res.output)
            snap = read_snapshot(res.output)
            @test snap.n == 128
            @test snap.vbc == 0.0
            @test 0.0 ≤ minimum(snap.dm_pos) && maximum(snap.dm_pos) < 1.0 + 1e-9
            # baryon (gas) particles lighter than DM by Ω_b/Ω_dm
            @test isapprox(snap.m_gas / snap.m_dm, 0.046 / (0.27 - 0.046); rtol = 1e-3)
            dv = CL.streaming_velocity(snap)
            @test maximum(abs, dv) < 1e-6          # no streaming
            @test abs(sum(snap.gas_delta)) / snap.n^3 < 1e-2   # mean δ_b ≈ 0
        end

        # ── v_bc = 30: coherent gas–DM bulk offset along one axis ───────────
        if get(ENV, "CICASS_SKIP_STREAMING", "0") != "1"
            @testset "generate v_bc=30 (streaming)" begin
                spec = CICASSSpec(boxlength = 0.2, zstart = 100.0, ngrid = 128, vbc = 30.0,
                                  filename = "test_v30")
                res = generate(spec)
                @test res.rc == 0
                snap = read_snapshot(res.output)
                dv = CL.streaming_velocity(snap)
                vexp = 30.0 * (1 + 100.0) / 1001.0          # ≈ 3.027 km/s phys
                # offset magnitude correct and concentrated on a single axis
                mag = sqrt(sum(abs2, dv))
                @test isapprox(mag, vexp; rtol = 5e-2)
                amax = argmax(abs.(dv))
                offaxis = sqrt(sum(abs2, dv) - dv[amax]^2)
                @test offaxis < 1e-3 * mag
            end
        end
    end
end
