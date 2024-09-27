# Monte Carlo Program to Simulate Patchy Spheres with Crowding Agents
**Read** - J. Chem. Theory Comput. 2024, 20, 18, 7700–7707 https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c01118

## Overview

This program simulates the behavior of patchy particles with crowding agents using the Monte Carlo (MC) method. Patchy particles have specific interaction sites (patches) that control how they interact with each other, facilitating the study of self-assembly processes. The simulation also includes crowding agents, which influence the assembly and energy of the system. The program tracks particle positions, patch orientations, energies, and bonding interactions.

## Files

1. **input_mc.dat**: Contains input parameters such as the box size, number of MC cycles, number of particles, crowding agents, interaction parameters, and more.
2. **ini_simu.xyz**: Initial configuration file for the simulation (positions of particles and patches).
3. **initial_simu.xyz**: Output file to log the initial configuration of the simulation.
4. **output.dat**: Output file that logs the main simulation parameters and energy states.
5. **accep_prob.dat**: Records the acceptance probabilities of Monte Carlo moves.
6. **structure_dist.dat**: Stores the distribution of bond structures and the occurrence of assembly states.
7. **traj.xyz**: File that logs the trajectory of the system at specific intervals during the simulation.

## Input Parameters

- **boxl**: Box length for the simulation box.
- **nmc**: Number of Monte Carlo steps to run.
- **natom**: Number of patchy particles in the system.
- **ncrowd**: Number of crowding agents in the system.
- **d_m**: Maximum displacement for particle moves.
- **dpsi_max**: Maximum angle for patch rotations.
- **R**: Gas constant (Boltzmann constant).
- **T**: Temperature in Kelvin.
- **theta0, maxtheta**: Angular parameters that control patch-patch interaction alignment.
- **eps_patch**: Strength of the patch-patch interaction.
- **eps_drive**: Energy change associated with the external drive.
- **patch_width**: Width of the patch-patch interaction zone.
- **eps11, eps12, eps22**: Energy parameters for Lennard-Jones interactions between particles and crowding agents.
- **sgma, sgmatwelve, sgmahex**: Sigma parameters for Lennard-Jones interactions between particles.
- **sgma12, sgmatwelve12, sgmahex12**: Sigma parameters for Lennard-Jones interactions between particles and crowding agents.
- **sgma22, sgmatwelve22, sgmahex22**: Sigma parameters for interactions between crowding agents.

## Program Structure

### Initialization

1. **Reading Input**: The program reads input parameters from `input_mc.dat` and initializes variables.
2. **Random Particle Initialization**: Patchy particles and crowding agents are placed randomly within the simulation box.
3. **Patch Initialization**: Patches are placed based on the position of the patchy particles and the interaction angles.

### Monte Carlo Simulation

1. **Energy Calculation**: 
   - The energy between patchy particles and crowding agents is computed using the Lennard-Jones potential.
   - Patch-patch interactions are computed based on angular alignment.
   
2. **Monte Carlo Moves**:
   - Random displacements and rotations are applied to the particles.
   - Energy changes are computed for both particle-particle and particle-crowding agent interactions.
   - The Metropolis criterion is applied to accept or reject each move based on the energy difference.

3. **Bond Formation**: 
   - Bonds are formed between particles if they meet specific distance and angular criteria.
   - Patch orientations are updated based on rotational moves.
   - The system tracks the number of bonds formed and their structure.

4. **Trajectory Output**: At regular intervals (e.g., every 5000 steps), the trajectory of the system is written to the `traj.xyz` file.

### Final Output

1. **Energy and Bond Structure**: 
   - After each Monte Carlo cycle, the total energy and bond structure are computed and logged.
   - The program tracks probabilities of different bond configurations and the occurrence of fully assembled states.

2. **Acceptance Ratios**: 
   - The acceptance ratios of particle moves are recorded in `accep_prob.dat`.

3. **Bond Distribution**: 
   - The distribution of bonds and the probability of different structures are written to `structure_dist.dat`.

## Running the Program

1. Prepare the input files (`input_mc.dat`, `ini_simu.xyz`).
2. Compile and run the program:
   ```
   gfortran mc.f90 -o mc
   ./mc
   ```
3. The program will output the following files:
   - **initial_simu.xyz**: Logs the initial configuration of particles and patches.
   - **output.dat**: Contains the main simulation parameters and energy states.
   - **accep_prob.dat**: Acceptance ratios of Monte Carlo moves.
   - **structure_dist.dat**: Distribution of bond structures and assembly states.
   - **traj.xyz**: Trajectory of the system during the simulation.

## Key Outputs

- **Energy States**: The initial, intermediate, and final energies of the system are logged.
- **Bond Configurations**: The bond structure of the patchy particles is tracked and analyzed.
- **Assembly Metrics**: The program tracks when the system reaches a fully assembled state and logs the order parameter (proportion of formed bonds to total possible bonds).

## Modifications

To modify the simulation:
- Adjust the patch interaction angles (`theta0`, `maxtheta`) to explore different assembly patterns.
- Change the interaction strengths (`eps_patch`, `eps_drive`) to simulate different physical systems.
- Modify the number of particles and crowding agents by adjusting `natom` and `ncrowd` in `input_mc.dat`.

## Contact

For any issues, please contact via email: shubhadeepnag92@gmail.com.
