=======================================================
Triangular Approximate SEmidefinite Relaxation (TASER) 
Massive MU-MIMO and JED in large SIMO Simulators
-------------------------------------------------------
(c) 2016 Christoph Studer & Oscar Castañeda
e-mail: studer@cornell.edu & oc66@cornell.edu
=======================================================

# Important information:

If you are using the simulator (or parts of it) for a publication, then you MUST cite our paper:

@article{CastanedaGoldsteinStuder:2016,
  author={O. Casta\~neda and T. Goldstein and C. Studer},
  journal={IEEE Transactions on Circuits and Systems I: Regular Papers},
  title={Data Detection in Large Multi-Antenna Wireless Systems via Approximate Semidefinite Relaxation},
  volume={63},
  number={12},
  pages={2334-2346},
  month={Dec.},
  year={2016}
}

and clearly mention this in your paper. More information about our research can be found at: http://vip.ece.cornell.edu

# How to start a simulation:

a. For large MU-MIMO wireless systems:

Simply run  

>> TASER_MIMO_sim

which starts a simulation in a 32 BS antenna, 32 user massive MU-MIMO system with QPSK modulation using MMSE and TASER detectors, as well as the SIMO lower bound. 

b. For JED in large SIMO wireless systems:

Simply run  

>> TASER_SIMO_JED_sim

which starts a simulation in a SIMO system with 16 BS antennas and transmission over 16 time slots with QPSK modulation using TASER detector. The simulator includes SIMO detection with both perfect receiver channel state information (CSI) and channel estimation (CHEST). 

In both simulators, and for smaller systems, you can also simulate ML and exact SDR detection. You can also do so for larger systems, but the simulation time is several hours! To use the exact SDR detector, you need to install CVX, which can be found here: 

http://cvxr.com/cvx/download/

Both the MIMO and SIMO simulators run with predefined simulation parameters. You can provide your own system and simulation parameters by passing your own "par"-structure (see the simulators for examples). 

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator. 

# Version 0.2 (November 27, 2016) - oc66@cornell.edu - re-expressing the code to match the algorithm as presented in the paper