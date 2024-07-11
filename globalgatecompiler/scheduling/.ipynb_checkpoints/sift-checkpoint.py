import qiskit
import networkx as nx

def circuit_to_dag(circuit):
    '''
    Inputs:
        - circuit: Qiskit circuit to convert to DAG
        - draw: if True, draws DAG labeled with node labels and gate type
    Outputs:
        - dag: networkx DAG representation of the input circuit
        - node_label_to_gate: Dict with keys given by DAG node labels (int) and values given by
                              gate (CircuitInstruction) corresponding to the node in the DAG. 
    '''
    
    dag = nx.DiGraph()
    
    qubit_paths = {}
    node_labels = {}
    node_label_to_gate = {}
    
    label = 0
    for op in circuit:
        dag.add_node(label)
        node_labels[label] = op[0].name
        node_label_to_gate[label] = op
        for qubit in op[1]:
            if qubit in qubit_paths:
                qubit_paths[qubit].append(label)
            else:
                qubit_paths[qubit] = [label]
        label += 1
                
    for qubit in qubit_paths:
        path = qubit_paths[qubit]
        for i in range(len(path)-1):
            if (path[i],path[i+1]) not in dag.edges:
                dag.add_edge(path[i],path[i+1])
        
    return dag,node_label_to_gate

def sift_commute_cz(c,node_label_to_gate):
    '''
    Specialized version of the general sift function given below, which takes into account the fact that
    CZ and CCZ gates can commute with each other. Assumes parallelism is maximized for single qubit gates,
    so that v_caught contains U3 gates and v_passed contains CZ and CCZ. If a given CZ or CCZ gate is not 
    disjoint from the set of qubits in v_remaining, but if the only gate it conflicts with is also a CZ or
    CCZ gate (and not a U3 gate), it can commute past and be placed in v_passed.
    
    Note: in theory, there may exist a circuit in which this produces a better schedule than the general 
    version of sift given below (which doesn't account for commutability of CZs and CCZs). In practice, 
    however, we find it almost always gives the same schedule. Thus, in most cases, it makes sense to use 
    the version below, which has slightly faster compile time complexity.
    '''
    
    v_passed = []
    v_caught = []
    v_remaining = []
    
    u3_qubits = set() # qubits from u3 ops in v_caught and v_rem
    cz_qubits = set() # qubits from cz ops in v_rem
    
    for node_label in c:
        op = node_label_to_gate[node_label]
        op_qubits = op.qubits
        
        if len(op[1])==1: # single-qubit gate
            v_rem_caught_qubits = u3_qubits.union(cz_qubits)
            if v_rem_caught_qubits.isdisjoint(op_qubits):
                v_caught.append(node_label)
            else:
                v_remaining.append(node_label)
            u3_qubits.update(op_qubits)
            
        else: # multi-qubit gate
            if u3_qubits.isdisjoint(op_qubits):
                v_passed.append(node_label)
            else:
                v_remaining.append(node_label)
                cz_qubits.update(op_qubits)
                
    return v_passed, v_caught, v_remaining

def sift(c,node_label_to_gate,num_qubits,indicator_fn=lambda op: len(op[1])==1):
    '''
    Given the remaining part of circuit that has not been scheduled yet, find the next groups of gates
    that can be scheduled. Gates with indicator function returning "True" that can be scheduled during 
    this round of sifting are "caught" by the sifter (added to the set v_caught). Gates of the other gate
    type that can be scheduled during this round are "passed" by the sifter (added to the set v_passed). 
    Gates in v_caught are scheduled to the same moment, and gates in v_passed are scheduled to the same 
    moment. All other gates are added to the set v_remaining and are scheduled during later sifting rounds.
    
    When in DAG form, C_caught, C_passed, and C_remaining form mutually disjoint subsets of the circuit C 
    input into the sift function. Concatenating C_caught, C_passed, and C_remaining gives the original 
    circuit C. The sets v_caught, v_passed, and v_remaining are comprised of the gates in the subcircuits
    C_caught, C_passed, and C_remaining, respectively. 
    '''
    
    v_passed = []
    v_caught = []
    v_remaining = []
    
    v_caught_rem_qubits = set()
    
    for node_label in c:
        op = node_label_to_gate[node_label]
        op_qubits = op.qubits
        if len(v_caught_rem_qubits)<num_qubits and v_caught_rem_qubits.isdisjoint(op_qubits):
            if indicator_fn(op):
                v_caught.append(node_label)
                v_caught_rem_qubits.update(op_qubits)
            else:
                v_passed.append(node_label)
        else: 
            v_remaining.append(node_label)
            v_caught_rem_qubits.update(op_qubits)
    
    return v_passed, v_caught, v_remaining

def get_sifted_schedule(circuit,indicator_fn=lambda op: len(op[1])==1):
    '''
    Given circuit, convert to DAG form and repeatedly call the sift function (which maximizes parallelism 
    of the specified gate type) until all gates have been scheduled. We assume the circuit (corresponding to 
    the input DAG) has previously been transpiled to the gate set {u3,CZ}, and that each moment in the output 
    schedule contains only one type of gate (i.e., a u3 gate and CZ gate cannot be scheduled into the same moment). 
    This is necessary to later decompose the circuit into the neutral atom gate set {Rz,GR,CZ}. 
    [Note, however, this code can be easily adjusted to account for other gate types if necessary.]
    
    Maximizing parallelism of u3 gates will minimize the total number of GR gates when the circuit is 
    decomposed to {Rz,GR,CZ}. This is motivated by the fact that GR gates often dominate gate durations, 
    especially on architectures where GR gates are implemented using microwave beams. However, maximizing
    u3 gate parallelism may result in less CZ gate parallelism. Thus, especially on architectures where CZ
    gate times are dominant, maximizing parallelism of u3 gates may not lead to the optimal circuit duration.
    
    Inputs:
        - circuit: Qiskit circuit to be scheduled.
        - indicator_fn: determines which gate type has parallelism maximized when scheduling.

    Outputs: 
        - schedule: list of moments in the order they appear in the schedule, where each moment is
                    a list of CircuitInstruction objects that are scheduled to in execute in parallel.
    '''
    n = len(circuit.qubits)
    
    schedule = []
    dag,node_label_to_gate = circuit_to_dag(circuit)
    v_remaining = list(nx.topological_sort(dag))
    
    while len(v_remaining)!=0:
        v_passed,v_caught,v_remaining = sift(v_remaining,node_label_to_gate,n,indicator_fn=indicator_fn)   
        schedule.append([node_label_to_gate[v] for v in v_passed])
        schedule.append([node_label_to_gate[v] for v in v_caught])
        
    if len(schedule[0])==0: # scheduling sometimes occurs in a way such that first moment is empty
        return schedule[1:] # return only the non-empty moments (empty moment irrelevant)
    return schedule