Configured by file {{CONFIG}} with a {{SITES}} site geometry.
There were {{BLOCKS}} blocks, each with {{SITESPERBLOCK}} sites (fluid and solid).
Recorded {{IMAGES}} images and {{SNAPSHOTS}} snapshots.
Ran with {{THREADS}} threads on {{MACHINES}} machines with {{DEPTHS}} deep topology.
Ran for {{STEPS}} steps of an intended {{TOTAL_TIME_STEPS}}.
With {{TIME_STEP_LENGTH}} seconds per time step.
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
Name Local Min Mean Max
{{#TIMER}}
{{NAME}} {{LOCAL}} {{MIN}} {{MEAN}} {{MAX}}
{{/TIMER}}

{{#BUILD}}
Revision number:{{REVISION}}
Steering mode: {{STEERING}}
Streaklines: {{STREAKLINES}}
Multimachine: {{MULTIMACHINE}}
Build type: {{TYPE}}
Optimisation level: {{OPTIMISATION}}
Built at: {{TIME}}
Reading group size: {{READING_GROUP_SIZE}}
Lattice: {{LATTICE_TYPE}}
Kernel: {{KERNEL_TYPE}}
Wall boundary condition: {{WALL_BOUNDARY_CONDITION}}

Communications options:
Point to point implementation: {{POINTPOINT_IMPLEMENTATION}}
All to all implementation: {{ALLTOALL_IMPLEMENTATION}}
Gathers implementation: {{GATHERS_IMPLEMENTATION}}
Separated concerns: {{SEPARATE_CONCERNS}}
{{/BUILD}}