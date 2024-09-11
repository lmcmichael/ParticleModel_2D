# ParticleModel_2D
Particle model code building upon (PM-ABL v1.0) in the McMichael et al. 2024 (GMD) paper. Code is extended to two dimensions (x, y). 

This code generates random ship trajectories given a number of ships (set by ens_num) and the number of particles injected at each time step (part_num). To avoid numerical issues associated with periodic boundary conditions the ships randomly reemerge inside the domain when encountering the domain edge; however, the particles emitted from the ships are allowed to evolve with periodic boundary conditions. 

This code is set up to run a repeated diurnal cycle from the CGILS case. The number of days to be simulated is given by num_days. The injection schedule is controlled by injection_duration and total_cycle_duration. The default value is to inject continuously. This can be configured to only inject during the daytime hours to simulate a more realistic MCB injection scenario (injection_duration .ne. 24).

Aerosol deposition is parameterized through particle decay as a function of particle age, with the decay speed being controlled by the e-folding timescale (tau_decay), which has a default value of 2 days. This means that after 2 days, the original 1000 particle population will be reduced to ~370 particles with the deleted particle positions being chosen at random. 

"Super"-particle densities (# / m2) are converted to particle concentrations (# / mg). Particle bin size is set to 1 km2 and results in relatively smooth concentration PDFs. The 1-km binned concentrations are then coarse-grained to 5 km2 bins to further smooth the fields. 

Various quantities are plotted and saved as .mat files. Particle concentrations are categorized by age (0-1 days, 1-3 days, 3+ days) and plotted as a function of time. Total particle concentration and particle concentration PDFs are also plotted as a function of time.  
