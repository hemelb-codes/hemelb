Configured by file {{CONFIG}} with a {{SITES}} site geometry.
Recorded {{IMAGES}} images and {{SNAPSHOTS}} snapshots.
Ran with {{THREADS}} threads on {{MACHINES}} machines with {{DEPTHS}} deep topology.
Ran for {{STEPS}} steps over {{CYCLES}} cycles.
With {{STEPS_PER_CYCLE}} per cycle.
{{#DENSITIES}}
!! Maximum relative density difference allowed {{ALLOWED}} was violated: {{ACTUAL}} !!
{{/DENSITIES}}
{{#UNSTABLE}}
!! Simulation was unstable !!
{{/UNSTABLE}}

Sub-domains info:
{{#DOMAIN}}
rank: {{RANK}}, fluid sites: {{SITES}}
{{/DOMAIN}}

Timing data:
Name Local Min Mean Max
{{#TIMER}}
{{NAME}} {{LOCAL}} {{MIN}} {{MEAN}} {{MAX}}
{{/TIMER}}
