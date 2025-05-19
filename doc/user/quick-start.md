
# Quickstart

This guide walks you through running a HemeLB simulation, and extracting simulation results on a Debian-based Linux machine. 

---

Assuming you already installed HemeLB and created HemeLB Geometry Tool environment and Python Tool environment. 

---

## 0. (Optional) Example File

If you want to use example files to get started, you can find it at `geometry-tool/tests/Model/data/`. 

---

## 1. Generate Geometry Files from STL

First of all, create your STL file of your intended geometry that you want to work on. Then you can use either GUI (preferable if you are using first time) or CLI (if you already have related files). 

Then activate the python environment where you installed the [geometry-tool](./geometry-tool.md).

### GUI:
```bash
hlb-gmy-gui --stl your_geometry.stl
```

You may face two errors in the command line (please refer [here](./geometry-tool.md) to know more), which you can ignore.

Then, 

- Mark inlets/outlets
- Set voxel size and seed point
- Click *Generate*
- [Optional] Save profile (`.pr2`)

### CLI:

If you have the profile (`.pr2`) file, and related STL file, you can use the command.

```bash
hlb-gmy-cli your_profile.pr2
```

Please ensure to keep the `.stl` file in the same directory or wherever it is referred inside the profile file.

If successfully executed, you will get something like this below:

```
Succesfully created closed polygon from input
The polyhedron has 160 facets 480 halfedges 0 border halfedges 82 vertices 
Preprocessing took: 0.000318 s 
Setup time: 0.198797 s
```

You should get your `.gmy` and `.xml` file after this. 

---

## 2. Prepare XML Configuration

Ensure your `input.xml` includes this below, to get your output file:

```xml
<properties>
  <propertyoutput file="whole.xtr" period="100">
    <geometry type="whole"/>
    <field type="velocity"/>
    <field type="pressure"/>
    <!-- You can also add more field types, please refer to XmlConfiguration.md for more. -->
  </propertyoutput>
</properties>
```

---

## 3. Run the Simulation

Now, run the command

```bash
mpirun -n number_of_processes hemelb -in input.xml -out ./output/
```

Here, you can put 1 (slower) or more (faster) as the number_of_processes. (For test file, 1 is enough)

After this, you should see result something like this below.

```
![1.3s]time step 2700
![1.4s]time step 2800
![1.4s]time step 2900
![1.5s]time step 3000
![1.5s]Finish running simulation.
```

If you do: Congratulations! You have successfully simulated using HemeLB. 

---

## 4. Post-Processing

You will get output files in `output` folder (based on previous command) and you can get `whole.xtr` file from there, which contains simulation output of your geometry. Afterwards, activate your python environment in which you installed the [python-tools](./python-tools.md). Then use the command below to convert `.xtr` output to a human readable `.csv` file:

```bash
hlb-dump-extracted-properties whole.xtr > whole.csv
```

After you get the CSV, you can view the results using your preferred tool (e.g., ParaView or Python).

---
