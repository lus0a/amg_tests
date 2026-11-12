# UG4 AMG test scripts

The Shuai experiments are available as:

- `amgtest_elder_Shuai.lua`
- `amgtest_elder_Shuai_fresh_amg.lua`
- `amgtest_henry_Shuai.lua`
- `amgtest_laplace_Shuai.lua`

Use this repository as `ug4/apps/amg_tests`, with the
[AMG plugin](https://github.com/lus0a/plugin_amg), ConvectionDiffusion,
and the AMG solver-descriptor support in
[ugcore](https://github.com/lus0a/ugcore). The time-dependent scripts also
require the time-integration classes provided by the Limex plugin.

Run the scripts from `ug4/apps/amg_tests`:

```sh
../../bin/ugshell -ex amgtest_henry_Shuai.lua
../../bin/ugshell -ex amgtest_laplace_Shuai.lua
```

The default grids for these two scripts are included in `grids/`.
`henry_2x2.ugx` is copied from the local `apps/Examples/grids` directory;
`unit_square_01_tri_2x2.ugx` is copied from the local `apps/intro_2` directory.
Both scripts accept `-grid <path>` and `-maxBase <nodes>`; Laplace also accepts
`-numRefs <levels>`, `-RAepsilon <value>`, and `-RAalpha <degrees>`.
The bundled cases are two-dimensional.

VTK output is written to the ignored `amgtestvtk/` directory, created by
the Henry and Laplace scripts. Matrix dumps are disabled by default;
use `-matrixPath <existing-directory>/` to enable them.

The Henry-named script preserves the existing local experiment's 600 by 150
geometry and boundary conditions. Its numerical solution has not been
validated against a reference Henry benchmark. With the default refinement
and `maxBase`, its 81-node system uses the direct base solver. Laplace uses
scalar rotated-anisotropic diffusion and exercises multiple AMG levels.

The Elder scripts additionally need
`../Examples/grids/elder_quads_8x2_bnd.ugx`. Before running either Elder script,
create its output directory with `mkdir -p amgtestvtk`.
