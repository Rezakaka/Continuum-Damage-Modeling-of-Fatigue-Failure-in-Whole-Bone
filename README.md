# Density-Based Bone Damage Simulation in Abaqus (UMAT + UAMP)

This repository contains a Fortran implementation of Abaqus user subroutines for simulating fatigue-driven damage accumulation in bone under cyclic loading.

The model combines:
- A UMAT subroutine for material behavior and damage evolution
- A UAMP subroutine for global stiffness monitoring and termination control

---

## Overview

This code models bone as an anisotropic elastic material with progressive stiffness degradation due to fatigue damage.

Key features:
- Spatially varying material properties (based on element ranges)
- Strain-driven damage accumulation
- Progressive reduction of Young’s modulus
- Automatic termination when global stiffness drops below a threshold

---

## Subroutine Structure

### 1. UAMP – Global Stiffness Monitoring

The UAMP subroutine monitors structural stiffness using a displacement sensor:
NODE_TRACKED_HISTORY


A constant load (1260.0) is divided by the measured displacement:


stiffness = load / displacement


- Initial stiffness is stored in `svars(1)`
- Current stiffness is stored in `svars(2)`

The stiffness ratio is computed as:


threshold = svars(2) / svars(1)


#### Termination Criterion

If:


stiffness ratio < 0.75


Then the simulation step is terminated:


lFlagsDefine(iConcludeStep) = 1


This corresponds to approximately a 25% reduction in global stiffness.

---

### 2. UMAT – Material Model and Damage Evolution

The UMAT subroutine defines:
- Elastic material behavior
- Damage accumulation
- Stiffness degradation

---

## Damage Model

Damage evolves using a strain-based fatigue law:


deltD = deltN * 33.70 * e_m**3.2410


Where:
- `deltD` = damage increment
- `deltN` = cycle increment
- `e_m` = von Mises equivalent strain
- `33.70`, `3.2410` = empirical constants

Total damage is accumulated as:


STATEV(2) = STATEV(2) + deltD


---

## Stiffness Degradation

The updated Young’s modulus is computed as:


STATEV(1) = Y_m * (1 - STATEV(2))


A minimum Young’s modulus is enforced:


E >= 1 MPa


to maintain numerical stability.

---

## Element-Based Material Assignment

Material properties are assigned based on element number ranges using the `get_Y` subroutine.

Example:


Eset(30,1) = 18409
Eset(30,2) = 18957
mat(30) = 25.2824


This means elements numbered from 18409 to 18957 are assigned a Young’s modulus of 25.2824 MPa.

This approach allows representation of spatially varying bone stiffness, typically derived from density mapping.

---

## Anisotropic Elasticity

The material is modeled as anisotropic using directional relationships:


E1 = 0.574 * E3
E2 = 0.577 * E3
G12 = 0.195 * E3
G23 = 0.265 * E3
G31 = 0.216 * E3


The compliance matrix is constructed and inverted (via `INV6X6`) to obtain the stiffness matrix.

---

## State Variables

The UMAT uses the following state variables:

- `STATEV(1)` → Updated Young’s modulus  
- `STATEV(2)` → Accumulated damage  
- `STATEV(3)` → Von Mises equivalent strain  
- `STATEV(4)` → Damage increment  
- `STATEV(5)` → Damage (1 - E / Y_m)  
- `STATEV(9:14)` → Stiffness matrix diagonal terms  

