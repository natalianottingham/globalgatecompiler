import qiskit

def transpile_to_u3_cz_ccz(circuit,optimization_level=3):  
    '''
    Transpile to gate set {u3,CZ,CCZ}. Assumes 3-q. gates in input circuit are CCZ, CCX, or CSWAP.
    (Qiskit transpilation will transpile to {u3,CZ}, even if CCZ is include in basis_gates)
    
    3-qubit gates are converted to a combination of CCZs + 1q gates in the following way:
        CCZs: remains the same
        CCX(c0,c1,t): H(t)-CCZ(c0,c1,t)-H(t)
        CSWAP(c,t1,t0): H(t1)-CCZ(c,t1,t0)-H(t0)-H(t1)-CCZ(c,t1,t0)-H(t0)-H(t1)-CCZ(c,t1,t0)-H(t1)
    Note that H is equivalent to U3(pi/2,0,pi)
    
    Blocks of 1-q & 2-q gates in between the 3-qubit gates are transpiled to {u3,CZ} using Qiskit.

    *** NOTE: 3-qubit gates will not be taken into account with Qiskit optimization passes to, e.g., 
        commute and cancel gates and/or merge single-qubit gates, as these passes are applied only to 
        the blocks of 1-qubit and 2-qubit gates in between; however, if there are opportunities in the 
        circuit to commute/merge gates, this is taken into account in our Theta-Opt scheduling code. ***
        
    '''
    new_circ = qiskit.QuantumCircuit(len(circuit.qubits))

    curr_block = qiskit.QuantumCircuit(len(circuit.qubits))

    for op in circuit:
        if len(op[1])>2:
            assert isinstance(op[0],qiskit.circuit.library.standard_gates.z.CCZGate) or \
            isinstance(op[0],qiskit.circuit.library.standard_gates.x.CCXGate) or \
            isinstance(op[0],qiskit.qiskit.circuit.library.standard_gates.swap.CSwapGate)

            c0,c1,target = [qubit.index for qubit in op[1]]

            if isinstance(op[0],qiskit.circuit.library.standard_gates.z.CCZGate): # CCZ
                transpile_and_append_curr_block(new_circ,curr_block,optimization_level)
                new_circ.ccz(c0,c1,target)
                curr_block = qiskit.QuantumCircuit(len(circuit.qubits))

            elif isinstance(op[0],qiskit.circuit.library.standard_gates.x.CCXGate): # CCX
                curr_block.h(target)
                transpile_and_append_curr_block(new_circ,curr_block,optimization_level)
                new_circ.ccz(c0,c1,target)
                curr_block = qiskit.QuantumCircuit(len(circuit.qubits))
                curr_block.h(target)
                
            else: # CSWAP
                curr_block.h(target)
                transpile_and_append_curr_block(new_circ,curr_block,optimization_level)
                new_circ.ccz(c0,c1,target)
                new_circ.u(np.pi/2,0,np.pi,c1)
                new_circ.u(np.pi/2,0,np.pi,target)
                new_circ.ccz(c0,c1,target)
                new_circ.u(np.pi/2,0,np.pi,c1)
                new_circ.u(np.pi/2,0,np.pi,target)
                new_circ.ccz(c0,c1,target)
                curr_block = qiskit.QuantumCircuit(len(circuit.qubits))
                curr_block.h(target)

        else: # two-qubit or single-qubit gate
            curr_block.append(op[0],[qubit.index for qubit in op[1]])

    if len(curr_block)>0:
        transpile_and_append_curr_block(new_circ,curr_block,optimization_level)
        
    return new_circ