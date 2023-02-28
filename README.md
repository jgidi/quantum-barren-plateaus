# __Avoiding Barren Plateaus in Variational Quantum Circuits__

A project by Felipe Quinteros (UdeC, Chile), Mariana Navarro (ICFO, Spain), Jorge Gidi (UdeC, Chile), and Luciano Pereira (IFF-CSIC, Spain) for the [Open Hackathon in the Qhack 2023](https://github.com/XanaduAI/QHack2023).

__Variational quantum algorithms__ have the potential to solve real-world practical problems with quantum computers in the near term [1]. They are hybrid quantum-classical algorithms that optimize an objective function encoded on a hamiltonian over a parametric quantum circuit. This objective function is efficiently evaluated in a quantum computer, while a classical computer is used to drive the optimization. The variational algorithms find applications in areas such as chemistry and finance.  

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig1_from_Noisy_intermediate_scale_quantum_algorithms.JPG?raw=true" width=500 style="display: block; margin: 0 auto" >

In order to reach real-word applications with variational algorithms, we have to overcome a big issue: the __barren plateau phenomenon__ [2,3]. This consists of the lack of convergence due to the vanishing of the cost function gradient. Pictorially, we can think that the landscape of the objective function is too flat, as shown in the figure, so any optimization method stuck. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig2_from_cost_function_induced_barren_plateau.JPG?raw=true" width=600 style="display: block; margin: 0 auto"  >

This problem can appear due to the following reasons:

1.- __The cost function is global__, i.e., involves measurements in many qubits [2]. For example, consider the $n$-qubits Hamiltonian 

$$H = I^{\otimes n}-|0^{\otimes n}\rangle\langle 0^{\otimes n}|.\qquad (1)$$
This Hamiltonian require a measurement over $|0^{\otimes n}\rangle$, wich involve all the $n$ qubits. Thereby, the sample of shots $N_{shots}$ used to evaluate $H$ has to split into $2^n$ outputs. In the case $N_{shots}<2^n$, which always happens for large $n$, the number of shots is not enough large to perform a precise evaluation of $H$, so that its gradient vanish and the barren plateau appears. This occurs independently of the shape and the depth of the variational circuit.

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig3_from_cost_function_induced_barren_plateau.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

In [2] the authors recommend avoiding this kind of Hamiltonian in the early iterations of the algorithm and instead employ a local hamiltonian $$H_L = I^{\otimes n}-\frac{1}{n}\sum_{j=1}^n I_{1}\otimes\cdots\otimes I_{j+1}\otimes|0_j\rangle\langle 0_j|\otimes I_{j-1}\otimes\cdots\otimes I_{n}.$$ Each term of this hamiltonian involves a single-qubits measurement, so that can be accurately evaluated with a small $N_{shots}$. The authors demonstrate that this kind of Hamiltonian can be used to train variational circuits with depth $\log(n)$.  

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig4_from_cost_function_induced_barren_plateau.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

2.- __The parametric quantum circuit is too expressive__, that is equivalent to saying that the parametric circuit approximates a 2-design [3]. When this condition is fulfilled the average gradient of the cost functions is null and its variance decay exponentially. $$\mathbb{E}( \nabla \langle H\rangle )=0, \qquad {\rm Var}( \nabla \langle H\rangle ) \sim \frac{1}{2^n}.$$ Almost any variational circuit composed of layers of local gates intercalated with layers of entangling gates satisfies this condition. Reducing the expressibility of variational circuits is a way to avoid barren plateaus, as is suggested in [4]. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig5_from_mitigating_barren_plateaus_of_vqes.JPG?raw=true" width=500 style="display: block; margin: 0 auto"  >

3.- __The initial condition for the optimization is too far from the solution__, which happens extremely often when is taken randomly. Having good initial conditions increase the probability to start the optimization with a non-zero gradient so that the protocol can converge to the optimum. For example, in [5] is proposed to classically pre-train the variational circuit with tensor networks to have a good initial condition for the quantum optimization. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig6_from_MPS_pretraining.JPG?raw=true" width=400 style="display: block; margin: 0 auto"  >

4.- __The problem of a lifetime: the noise__ [6]. The impact of the noise can be reduced by hardware-efficient circuits and local measurements. 

<img src="https://github.com/jgidi/quantum-barren-plateaus/blob/main/imgs/fig7_from_noise_induced_BP.JPG?raw=true" width=400 style="display: block; margin: 0 auto"  >

In this project, we develop a __Qiskit module__ that includes several proposals to reduce the impact of barren plateaus in variational quantum algorithms. We also provide an early implementation on __Pennylane__. For more details visit the introductory notebook or each individual tutorial notebook.

[1] Noisy intermediate-scale quantum algorithms, Rev. Mod. Phys. 94, 015004 (2002).

[2] Cost function dependent barren plateaus in shallow parametrized quantum circuits, Nat Commun 12, 1791 (2021).

[3] Barren plateaus in quantum neural network training landscapes, Nat Commun 9, 4812 (2018).

[4] Mitigating barren plateaus of variational quantum eigensolvers, arXiv:2205.13539v2 (2022).

[5] Matrix product state pre-training for quantum machine learning, Quantum Sci. Technol. 7, 035014 (2022).

[6] Noise-induced barren plateaus in variational quantum algorithms, Nat Commun 12, 6961 (2021).

Tutorials:

1.- [Introduction to QuBaP](https://github.com/jgidi/quantum-barren-plateaus/blob/d89ac2b072019616ebb313357f6346980e6d42ca/Tutorials/Introduction_to_QuBaP.ipynb)

2.- [Avoiding Cost-function induced Barren Plateaus](https://github.com/jgidi/quantum-barren-plateaus/blob/d89ac2b072019616ebb313357f6346980e6d42ca/Tutorials/CostFunction_BarrenPlateaus.ipynb)

3.- [State Efficient Ansatz](https://github.com/jgidi/quantum-barren-plateaus/blob/932915a8d6a35c4e6cbca76322e7f5a825a4e000/Tutorials/State_Efficient_Ansatz_BeH2.ipynb)

4.- [Classical Pretraining](https://github.com/jgidi/quantum-barren-plateaus/blob/3ab5be7d7c68b66296c6460ba1fd1bb7c8f4d8c8/Tutorials/MPS_pretraining.ipynb)

5.- [Combining the methods](https://github.com/jgidi/quantum-barren-plateaus/blob/d89ac2b072019616ebb313357f6346980e6d42ca/Tutorials/Combining_methods.ipynb)

6.- [Early implementation on Pennylane](https://github.com/jgidi/quantum-barren-plateaus/blob/dceb5d06908d3c97c6324d5bc29d48592db81e09/Tutorials/Tutorial_early_pennylane.ipynb)

Install via `pip` as
``` sh
pip install git+https://github.com/jgidi/quantum-barren-plateaus
```

or, equivalently,
``` sh
python -m pip install git+https://github.com/jgidi/quantum-barren-plateaus
```
