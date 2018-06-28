MATLAB code for JAK/STAT case study. 

jakstat2.m calculated orthogonal basis vectors using Hermite polynomials. 

nonintrusiveorder2.m uses nonintrusive gPC to calculate time evolution of all 4 states with no noise, and ekf_jakstat.m propagates through the stochastic noise term of the system using an extended Kalman filter. jaktotalvar.m calculates the total variance when incorporating both nonintrusive gPC and the extended Kalman filter. This is compared to visualizations from jakstatbenchmark.m. Visualization data came from iPython notebook which is also in this GitHub. Everything else is saved .mat files that make computation and troubleshooting easier and other .m files that were used for method verification.
