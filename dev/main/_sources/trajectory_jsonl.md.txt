# OPTIMADE trajectory JSON Lines

`httk-trajectory-jsonl` is a streaming holding format. The public filename
convention is `.traj.jsonl`; the loader registers `.jsonl` because core strips
one compression suffix and dispatches on the remaining final suffix. Thus
`.traj.jsonl.gz` works through the normal text datastream/compression path.

The first line is an OPTIMADE 1.2.0 dense partial-data header with an
`x-httk-trajectory` description. Subsequent lines are frame objects containing
`index`, `fractional_site_positions`, and `observables`; variable-cell files
also contain `lattice_vectors`. See the `httk.io.optimade_jsonl` module
docstring for the normative schema. Values are float64 presentation values;
there is no exact-token channel.

Binary framing is intentionally not included. Its framing and random-access
trade-offs remain a separate design decision.

`TrajectoryJsonlFile.path` returns the source filename string used to construct
the lazy reader.
