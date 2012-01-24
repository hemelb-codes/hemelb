Configured by file {{CONFIG}} with a  {{SITES}} site geometry.
Recorded {{IMAGES}} images and {{SNAPSHOTS}} snapshots.
Ran with {{THREADS}} threads on {{MACHINES}} machines with {{DEPTHS}} deep topology.
Ran for {{STEPS}} steps over {{CYCLES}} cycles.
With {{STEPS_PER_SECOND}} time steps per second and {{STEPS_PER_CYCLE}} per cycle.
{{#DENSITIES}}
!! Maximum relative density difference allowed {{ALLOWED}} was violated: {{ACTUAL}} !!
{{/DENSITIES}}
{{#UNSTABLE}}
!! Simulation was unstable !!
{{/UNSTABLE}}

Sub-domains info:
{{#PROCESSOR}}
rank: {{RANK}}, fluid sites: {{SITES}}
{{/PROCESSOR}}

Timing data:
Name Local Min Mean Max Per
{{#TIMER}}
{{NAME}} {{LOCAL}} {{MIN}} {{MEAN}} {{MAX}} {{NORMALISATION}}
{{/TIMER}}
