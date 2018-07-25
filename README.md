# ADCS Simulation


## Features

### **Current features**

- [x] LEO Orbit propagation
- [x] 3-axis permanent magnet + hysteresis rods
- [x] Solar panel heat flux + surface temprature
- [x] Pertubation: Gravity gradient
- [x] Pertubation: Global magnetic model

### **Future implementation**

- [ ] Pertubation: Solar flux
- [ ] Pertubation: Aerodynamic drag
- [ ] Compare with SNAP simulator 
- [ ] GUI
- [ ] Automatically calculate period 

## Usage (Matlab Sims)

### **Inputs**
1. Sim duration [ s ]
2. Intial day 
3. Hysteresis rods
    * Rods per axix ( x + y )
    * Rod length + radius [ m ] 
    * Coercivity [ kA / m ]
    * Magnetization at saturation [ Tesla ] 
    * Magnetic remanence [Tesla ]
4. Permanent magnet
    * Magnetic moment [ A^2 * m ]
5. Orbit
    * Major Axis [ m ]
    * Eccentricity  
    * Inclination [ deg ]
    * Long. Of Ascending Node [ deg ]
    * Argument Of Periapsis [ deg ] 
    * Mean Anomaly at Epoch [ deg ]
    * Eclipse ratio
    * Period [ s ]
6. Satellite
    * Mass Moment of Intertia [ kg * m^4 ]

### **Output**
1. Angular Veocity plot
2. Error Angle plot
3. Torque Curves
4. Hysteresis Curves
5. Inertial Orbit graphic 
6. Earth Orbit graphic 
7. Temprature plot
8. Flux plot


## Authors

* Addy Bhatia (addy.bhatia@mail.utoronto.ca)

