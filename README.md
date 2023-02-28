# Avoiding Barren Plateaus in Variational Quantum Circuits

The barren plateau phenomenon is one of the main issues to solve in order to attain practical applications of variational quantum algorithms. This consists of the lack of convergence due to the vanishing of the cost function gradient. In this project, we implement several proposals to avoid barren plateaus in variational algorithms.


Variational quantum algorithms have the potential to solve real-world practical problems with quantum computers in the near term. They are hybrid quantum-classical algorithms that optimize an objective function encoded on a hamiltonian over a parametric quantum circuit. This objective function is efficiently evaluated in a quantum computer, while a classical computer is used to drive the optimization. The variational algorithms find applications in areas such as chemistry and finance.  

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig1_from_Noisy_intermediate_scale_quantum_algorithms.JPG?raw=true" width=500 style="display: block; margin: 0 auto" >

In order to reach real-word applications with variational algorithms, we have to overcome a big issue: the barren plateau phenomenon. This consists of the lack of convergence due to the vanishing of the cost function gradient. Pictorially, we can think that the landscape of the objective function is too flat, as shown in the figure, so any optimization method stuck. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig2_from_cost_function_induced_barren_plateau.JPG?raw=true" width=600 style="display: block; margin: 0 auto"  >

This problem can appear due to the following reasons:

1.- The cost function is global, i.e., involves measurements in many qubits. For example, consider the $n$-qubits Hamiltonian 
\begin{equation}
H = I^{\otimes n}-|0^{\otimes n}\rangle\langle 0^{\otimes n}|.
\end{equation}
This Hamiltonian require a measurement over $|0^{\otimes n}\rangle$, wich involve all the $n$ qubits. Thereby, the sample of shots $N_{shots}$ used to evaluate $H$ has to split into $2^n$ outputs. In the case $N_{shots}<2^n$, which always happens for large $n$, the number of shots is not enough large to perform a precise evaluation of $H$, so that its gradient vanish and the barren plateau appears. This occurs independently of the shape and the depth of the variational circuit.

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig3_from_cost_function_induced_barren_plateau.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

In [] the authors recommend avoiding this kind of Hamiltonian in the early iterations of the algorithm and instead employ a local hamiltonian $$H_L = I^{\otimes n}-\frac{1}{n}\sum_{j=1}^n I_{1}\otimes\cdots\otimes I_{j+1}\otimes|0_j\rangle\langle 0_j|\otimes I_{j-1}\otimes\cdots\otimes I_{n}.$$ Each term of this hamiltonian involves a single-qubits measurement, so that can be accurately evaluated with a small $N_{shots}$. The authors demonstrate that this kind of Hamiltonian can be used to train variational circuits with depth $\log(n)$.  


<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig4_from_cost_function_induced_barren_plateau.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

2.- The parametric quantum circuit is too expressive, that is equivalent to saying that the parametric circuit approximates a 2-design. When this condition is fulfilled the average gradient of the cost functions is null and its variance decay exponentially. $$ \mathbb{E}( \nabla \langle H\rangle )=0, \qquad {\rm Var}( \nabla \langle H\rangle ) \sim \frac{1}{2^n}.$$ Almost any variational circuit composed of layers of local gates intercalated with layers of entangling gates satisfies this condition. Reducing the expressibility of variational circuits is a way to avoid barren plateaus, as is suggested in []. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig5_from_mitigating_barren_plateaus_of_vqes.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

3.- The initial condition for the optimization is too far from the solution, which happens extremely often when is taken randomly. Having good initial conditions increase the probability to start the optimization with a non-zero gradient so that the protocol can converge to the optimum. For example, in [] is proposed to classically pre-train the variational circuit with tensor networks to have a good initial condition for the quantum optimization. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig6_from_MPS_pretraining.JPG?raw=true" width=400 style="display: block; margin: 0 auto"  >

4.- The problem of a lifetime: the noise []. The impact of the noise can be reduced by hardware-efficient circuits and local measurements. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig7_from_noise_induced_BP.JPG?raw=true" width=400 style="display: block; margin: 0 auto"  >




Install via `pip` as
``` sh
pip install git+https://github.com/jgidi/quantum-barren-plateaus
```

or, equivalently,
``` sh
python -m pip install git+https://github.com/jgidi/quantum-barren-plateaus
```
