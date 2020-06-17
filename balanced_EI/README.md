# BalancedEIMinimax
Code for "Minimax Dynamics of Optimally Balanced Spiking Networks of Excitatory and Inhibitory Neurons".
The four subfolders contains codes/data needed to reproduce each figure in the main text.
* Figure 1: 
  * To reproduce panels in figure 1 with data stored in Data folder
  ```
  plot_figure_1_panels.m
  ```
  * To simulate tightly balanced network with a given $\tau$
  ```
  [rate,re,ri,W_EE,W_EI,W_IE,W_II,F,NE,N,tau] = balanced_network(tau)
  ```
  * For firing rate prediction accuracy
  ```
  ropt = fr_prediction(tau,rate,re,ri,W_EE,W_EI,W_IE,W_II,F,NE,N)
  ```
* Figure 2:
  * To reproduce panels in figure 2 with data stored in Data folder
  ```
  plot_figure_2_panels.m
  ```
  * To simulate network for signal reconstruction of spike trains with constant firing rate, with NE excitatory neurons
  ```
  [] = signal_reconstruction_synthetic(NE,0)
  ```
  * To simulate network for signal reconstruction of time varying input
  ```
  [] = signal_reconstruction_synthetic(NE,1)
  ```
  * To simulate network for signal reconstruction of natural image patches n with NI inhibitory neurons
  ```
  [] = signal_reconstruction_image(n,NI)
  ```
* Figure 3:
  * To reproduce panels in figure 3 with data stored in Data folder
  ```
  plot_figure_3_panels.m
  ```
  * To simulate network with learned weights run
  ```
  [] = fixed_point_attractor(index, noise)
  ```
* Figure 4:
  * To reproduce panels in figure 4 with data stored in Data folder
  ```
  plot_figure_4_panels.m
  ```
  * To simulate network that stores ring attractor
  ```
  [] = spiking_ring_attractor(theta_1,theta_2,gamma)
  ```
  * Simulate network for grid attractor
  ```
  [] = spiking_grid_attractor()
  ```
For meanings of the parameters in the above example code type 
 ```
 help function_name
 ```
