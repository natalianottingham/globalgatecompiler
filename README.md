# Global Gate Compiler
The code in this repository corresponds to our paper *"Circuit decompositions and scheduling for neutral atom architectures with limited local addressibility"* which will appear in QCE 2024.

We provide a compiler framework to convert an arbitrary input circuit into a native gate set that involves global gates, with a focus on the decomposition and scheduling steps. Our approach specifically minimizes the number of global gates and total global gate rotation amount in the final circuit. Our work is motivated by the fact that many current neutral atom architectures do not natively support local addressing of single-qubit rotations about an axis in the xy-plane of the Bloch sphere. We assume a native gate set of local Rz($\lambda$) rotations of arbitrary angle $\lambda$, local CZ and CCZ entangilng gates, and global GR($\theta,\phi$) gates that simultaneously rotate all qubits in the circuit by an angle $\theta$ about an axis parameterized by $\phi$. However, our code can easily be adjusted to work with other native gate sets that involve global gates. 

### <ins>INSTALLATION</ins>: 
To install the most recent version:

```bash
python3 -m pip install git+https://github.com/natalianottingham/globalgatecompiler@master
```

The following packages are required to run our code:
```bash
qiskit 1.1.1
numpy 2.0.0
networkx 3.3
```

### <ins>GETTING STARTED</ins>:
See the provided `tutorial_notebook.ipynb` for detailed examples and explanations of how to use our code.

### <ins>CITING OUR CODE</ins>:
If you use our code in your work, please cite the corresponding paper:<br>
[will add citation info once publication comes out]

For further questions/comments about our code: nottingham@uchicago.edu