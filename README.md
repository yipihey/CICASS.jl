# CICASS.jl

Julia bindings to **CICASS** (McQuinn & O'Leary's cosmological initial-conditions
generator that self-consistently includes the baryon–dark-matter **streaming
velocity**; [arXiv:1204.1344](https://arxiv.org/abs/1204.1344),
[arXiv:1204.1345](https://arxiv.org/abs/1204.1345)), in the same library style as
the other codes in the EnzoNG.jl unified multi-code framework
(Enzo / RAMSES / Arepo / MUSIC over `CodeBridge`).

CICASS is wrapped as a **C-ABI shared library** (`libcicass_capi.dylib`, built from
a fork of the CICASS source at [`yipihey/CICASS`](https://github.com/yipihey/CICASS),
a fork of [`astromcquinn/CICASS`](https://github.com/astromcquinn/CICASS)) plus an
HDF5-free `.cicass` raw dump, so a Julia-driven session can generate
streaming-velocity ICs in-process (or from a CodeBridge worker) and inject them
into live Enzo and RAMSES.

## Why CICASS, alongside MUSIC and DISCO-DJ

MUSIC's two-component (gas + dark matter) ICs share identical phases and differ
only in amplitude — no relative displacement, no bulk velocity offset. DISCO-DJ's
LPT is single-component. CICASS uniquely carries the Tseliakhovich–Hirata
**streaming velocity** `v_bc`: a 2D `(k⊥, k∥)` transfer function evolved forward
from recombination that gives the gas a coherent **bulk velocity offset** relative
to the dark matter (≈ `v_bc·(1+z)/1001` km/s, on one axis).

It is a **post-recombination** model — useful for `zstart ≲ 1000` (typically
z ≈ 100–200); it cannot initialize at pre-recombination redshifts.

## Layout

```
lib/CICASSLib   The CodeBridge wrapper: Bridge(:cicass) over libcicass_capi,
                CICASSSpec → make_tf (transfer.x 2D v_bc TF grid) → generate
                (.cicass raw dump) → read_snapshot, and a serve() worker loop.
```

The cross-code streaming gate (one `CICASSSpec` → live Enzo ≡ live RAMSES, both
carrying the streaming offset) lives in EnzoNG.jl's `MultiCodeCICASSExt`.

## Build

```
bash cicass/deps/build_cicass_darwin.sh   # from a sibling CICASS source checkout
                                          # → cicass/build/{transfer.x, libcicass_capi.dylib}
```

Needs Homebrew `gsl` + `fftw` and Apple clang. `CICASSLib.available()` gates every
live call; set `ENV["CICASS_LIB"]` to point at a prebuilt library.

## Test

```
<julia> --project=lib/CICASSLib/test lib/CICASSLib/test/runtests.jl
```
