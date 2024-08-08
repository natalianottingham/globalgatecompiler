import qiskit
import networkx as nx

from .schedule_class import Schedule

def circuit_to_dag(circuit):
    '''
    Inputs:
        - circuit: Qiskit circuit to convert to a directed acyclic graph (DAG).
    Outputs:
        - dag: networkx DAG representation of the input circuit.
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

def sift(c,node_label_to_gate,num_qubits,indicator_fn=lambda op: len(op[1])==1):
    '''
    Given the remaining part of circuit that has not been scheduled yet, find the next two moments of gates
    following the steps outlined in Alg.1 and Sec. V-B of our paper. 

    Inputs:
        - c: the remaining part of the circuit that has yet to be scheduled.
        - node_label_to_gate: Dict with keys given by DAG node labels (int) and values given by
                              gate (CircuitInstruction) corresponding to the node in the DAG. 
        - num_qubits: the number of qubits in the circuit.
        - indicator_fn: defines the categories of gates, where gates in different categories cannot 
                        be scheduled to the same moment, and parallelism is maximized for the 
                        gate type with indicator function value of 1 (True). 
    Outputs:
        - v_passed: moment of gates such that 0 (False) is returned when input into the indicator function
                    and all dependencies have already been scheduled; v_passed is added as the next moment
                    in the circuit.
        - v_caught: moment of gates such that 1 (True) is returned when input into the indicator function
                    and all dependencies have already been scheduled (either in previous sift calls or in
                    the v_passed moment returned by this sift call); v_caught is added as the next moment
                    after v_passed.
        - v_remaining: all remaining gates which will be scheduled by later calls to sift.
    
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
    Schedules a circuit using repeated calls to sift, as described in Sec. V-B of our paper.
    If indicator_fn is set to the default, single-qubit gate parallelism is maximized; assuming the
    circuit has previously been transpiled to {u3,CZ,CCZ}, maximizing u3 parallelism during scheduling
    will result in a final circuit with the minimum possible number of GR gates, once decomposed to the
    gate set {Rz,GR,CZ,CCZ}. 
    
    Inputs:
        - circuit: Qiskit circuit to be scheduled.
        - indicator_fn: defines the categories of gates, where gates in different categories cannot 
                        be scheduled to the same moment, and parallelism is maximized for the 
                        gate type with indicator function value of 1 (True). 

    Outputs: 
        - schedule: list of moments in the order they appear in the schedule.
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
        return Schedule(schedule[1:],circuit) # return only the non-empty moments (empty moment irrelevant)
    return Schedule(schedule,circuit)