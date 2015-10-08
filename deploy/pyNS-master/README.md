pyNS
========

A modular solver framework for 0D/1D problems. Stable release, current version is 0.4.2

##Software architecture

pyNS has been designed with an object-oriented approach, which allows to abstracts the concept of element from the numerical solver itself. 

##Installation Requirements:

- [Python 2.6.x or 2.7.x](http://www.python.org/)
- [Numpy package](http://numpy.scipy.org/) (be careful to download the version related to your python version)
- [Scipy package](http://www.scipy.org)  ONLY FOR WINDOWS USERS

A browser for visualizing simulation results (internet connection is not required). Currently supported browser are Firefox, Chrome, Safari, Opera and Internet Explorer.

For additional features:
- [MEncoder library](http://www.mplayerhq.hu/design7/dload.html) (Only for velocity profile videos)
- [MatPlotLib library](http://matplotlib.sourceforge.net) (Post processing using png files instead of default browser visualization and for velocity profile videos.)
- [LXML package](http://lxml.de/) (Only for xml files validation feature).

### How to

Open a terminal, cd into pyNS directory and type:
[For Mac/Linux Users]
```bash
python pyNS.py --help 
```

[For Windows Users]
```bash
pyNS.py --help
```

Standard benchmark simulation:
By default pyNS runs a vascular network which represents the arterial vasculature of a right arm.
Open a terminal cd into pyNS directory and type:
[For Mac/Linux Users]
```bash
python pyNS.py --default
```

[For Windows Users]
```bash
pyNS.py --default
```
## Documentation
For a complete documentation please go to http://archtk.github.com

### The vascular network model

The vascular network model implemented in pyNS is a network graph consisting of segments and nodes serially connected on the basis of the anatomical configuration. Each segment is discretized in one or more elements, each of them modeled as an electronic circuit using the hydraulic analogy, with its own relation for pressure drop (p) and flow volume (q). Since elements are independent, we just sum across them in different ways for different circuit configurations. Thus, each kind of element represents a mathematical model which is implemented for the scope of modeling and studying pressure, volumetric flow rate and wall-shear stress distributions over the complete vascular system or a specific part of it. Each segment of the network only exposes its external nodes and each element has its own local degrees of freedom depending of the nature of the element itself. If an element is connected to one or more elements the connection node and its local degree of freedom are shared between elements. This structure allows to build complex vascular network, which can be potentially composed by different types of segments, ensuring mass and momentum conservation over the entire network.

## Guidelines for templates simulations

Guidelines for parameters.csv file, needed for arm templates simulations:

* Date of birth (dob) has to be in dd/mm/yyyy
* Date of surgery (dos) has to be in dd/mm/yyyy 
* Gender 0 Female, 1 Male
* Arm 0 Left, 1 Right
* ftype 0 End-to-end radio-cephalic fistula, 1 End-to-side radio-cephalic fistula, 3 End-to-side brachio-cephalic fistula, 5 End-to-side brachio-basilic fistula
* Diabetes 0 No, 1 Yes
* Hypertension 0 No, 1 Yes
* Height has to be in cm
* Weight has to be in Kg
* Systolic and diastolic pressures (sys and diap) have to be in mmHg
* Cardiac output is not mandatory but has to be in mL/min
* Period is the cardiac period, measured in seconds
* Brachial, radial and ulnar flows have to be in mL/min
* Hematocrit (ht) has to be in %
* Protein plasma concentration (cp) has to be in g/dL

If any parameter is not specified, pyNS assumes default values:

* dos = 27/07/2010
* dob = 27/07/1960
* gender Male
* Arm Right
* Fistula Type End-to-side radio-cephalic
* Height 175 cm
* Weight 70 Kg
* Systolic Pressure 110 mmHg
* Diastolic Pressure 70 mmHg
* Brachial Flow 130 mL/min
* Radial Flow 20 mL/min
* Ulnar Flow 30 mL/min
* Period 1 second
* Hematocrit 45 %
* Protein Plasma Concentration 7 g/dL
* Hypertension No
* Diabetes No

## License

Copyright (c) Simone Manini, Luca Antiga. 
Distributed under the BSD License.