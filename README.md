# Data Driven Fluid Mechanics

### by
Amir Gholizad

---

## Project Overview

This project, submitted in partial fulfillment of the requirements for the course CMSC 6920, investigates the behavior of fluid flow using Proper Orthogonal Decomposition (POD) and Dynamic Mode Decomposition (DMD). The study focuses on the vorticity field of a two-dimensional flow past a cylinder at Reynolds number $Re = 200$.

### Instructor
Dr. Jahrul Alam

### Department
Department of Mathematics and Statistics, Memorial University

### Date
April 2023

### Location
St. John’s, Newfoundland and Labrador, Canada

---

## Abstract

The aim of this project is to analyze the fluid flow past a cylinder using POD and DMD methods. Data was collected on a $199×449$ grid at time intervals of $0.125$ seconds. Each snapshot is represented as a vector of length 89351, capturing the flow’s characteristics. The combined POD and DMD methods were used to predict the future modes of the flow.

---

## Table of Contents

1. [Introduction](#introduction)
    - Proper Orthogonal Decomposition (POD)
    - Dynamic Mode Decomposition (DMD)
2. [Methodology](#methodology)
    - Proper Orthogonal Decomposition (POD)
    - Dynamic Mode Decomposition (DMD)
3. [Prediction](#prediction)
4. [Simulation and Results](#simulation-and-results)
5. [Bibliography](#bibliography)

---

## Introduction

### Proper Orthogonal Decomposition (POD)

POD is used in fluid mechanics to analyze data that depends on both space and time, $U(x, t)$. It separates the variables into spatial modes $\phi_{j}(x)$ and temporal modes $a_{j}(t)$, allowing the representation of $U(x, t)$ as a sum of these modes. 

### Dynamic Mode Decomposition (DMD)

DMD identifies and analyzes dominant spatiotemporal patterns and dynamics within high-dimensional, time-varying datasets. It reconstructs and predicts future states of the system, often incorporating POD to identify dominant spatial modes.

---

## Methodology

### Proper Orthogonal Decomposition (POD)

The event $U(x, t)$ is separated into spatial and temporal modes. Data is assembled into a matrix with each column representing a snapshot at a given time and each row consisting of measurements across all times. The matrix undergoes Singular Value Decomposition (SVD) to obtain spatial modes and their corresponding temporal coefficients.

### Dynamic Mode Decomposition (DMD)

Each snapshot is related to the previous state using the equation $U_{k+1} = AU_{k}$, where $A$ is a linear operator. The matrix $A$ is solved using the reduced SVD of $U$, allowing the identification of eigenvalues and eigenvectors, leading to the calculation of Dynamic Modes.

---

## Prediction

The differential equation for fluid flow, $\frac{\partial U(x, t)}{\partial t} = LU$, is solved to predict future states. The eigenvalues of $L$ are determined, allowing the prediction of the system’s time evolution using $U(x, t) = \phi_{p} b e^{λt}$.

---

## Simulation and Results

Data is loaded into a Python script and reshaped. Random snapshots are visualized to understand the amplitude and importance of the modes. The energy of each mode is calculated using the squared $L_{2}$ norm of normalized singular values. Error levels for different truncation numbers are evaluated to determine the optimal number of modes. 

Dynamic Modes are calculated and compared with Projected modes to ensure convergence and stability. Predictions are made using the developed model, achieving high correlation with the actual data.

### Key Results
- Energy distribution of modes
- Error levels for different truncation numbers
- Dynamic Modes visualization
- Predicted vs. actual snapshots

---

## Bibliography

1. Miguel A. Mendez, Andrea Ianiro, Bernd R. Noack, and Steven L. Brunton, editors. *Data-Driven Fluid Mechanics*. Cambridge University Press, January 2023.
2. Xuan Dai, Da Xu, Mengqi Zhang, and Richard JAM Stevens. *A three-dimensional dynamic mode decomposition analysis of wind farm flow aerodynamics*. Renewable energy, 191:608–624, 2022.
3. Jonathan H Tu. *Dynamic mode decomposition: Theory and applications*. PhD thesis, Princeton University, 2013.