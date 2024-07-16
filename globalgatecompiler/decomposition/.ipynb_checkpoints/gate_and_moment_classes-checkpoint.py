import qiskit
from qiskit.circuit import Gate
from typing import Dict
from itertools import product
import numpy as np

class GRGate(Gate):
    '''
    Global rotation gate, applied to all qubits in the circuit.
        - num_qubits: number of qubits in the circuit
        - theta: defines the rotation amount
        - phi: defines axis of rotation within the xy-plane
    '''
    
    def __init__(self, num_qubits: int, theta: float, phi: float)->None:
        name = f"GR"#({theta:.2f}, {phi:.2f})"
        super().__init__(name,num_qubits,[theta,phi])
        
    def _define(self):
        theta = self.params[0]
        phi = self.params[1]
        
        qr = qiskit.QuantumRegister(self._num_qubits)
        qc = qiskit.QuantumCircuit(qr, name=self.name)
        rules = [
            (qiskit.circuit.library.RGate(theta,phi), [qr[i]], []) for i in range(len(qr))
        ]
        for instr, qargs, cargs in rules:
            qc._append(instr, qargs, cargs)
        
        self.definition = qc
        
class GRMoment:
    '''
    A global gate moment. Contains only one GR gate. 
    '''
    
    def __init__(self,theta:float,phi:float,num_qubits_in_circuit:int) -> None:
        self.theta = theta
        self.phi = phi
        
        self.qubits = set([q for q in range(num_qubits_in_circuit)])
        self.duration = None
        self.matrix = None
        
    def set_duration(self,gate_durations,connectivity_graph=None):
        self.duration = gate_durations['gr']*abs(self.theta)/np.pi
        
class RzMoment:
    '''
    A group of Rz gates that can all be executed in parallel. Duration of the entire moment 
    is determined by the longest gate (largest angle of rotation) within the moment. 
    '''
    
    def __init__(self,angles:Dict[int,int]) -> None:
        self.angles = angles
        
        self.qubits = set([q for q in list(self.angles.keys())])
        self.duration = None
        self.matrix = None
    
    def set_duration(self,gate_durations,connectivity_graph=None):
        self.duration = gate_durations['rz']*max([abs(angle) for angle in list(self.angles.values())])/np.pi

class MultiQubitGateMoment:
    '''
    A group of CZ and CCZ gates that occur consecutively between other Rz and GR moments. 
    Note: this group of gates may or may not all be executable in parallel;
    parallelism constraints between multi-qubit gates are enforced in set_duration.
    '''
    
    def __init__(self,qubit_tuples:list) -> None:
        '''
        - self.qubit_tuples: list of tuples, where each tuple is a subset of qubits acted upon by a gate
        - self.submoments: list of submoments, where each submoment is a group of parallel CZs & CCZs
        '''
        self.qubit_tuples = qubit_tuples
        self.submoments = None
        
        self.duration = None
        self.cz_duration = None
        self.ccz_duration = None
        
        self.qubits = set(list(sum(qubit_tuples, ())))
        self.matrix = None
        
    def set_duration(self,gate_durations,connectivity_graph):
        '''
        Determine which CZ & CCZ gates can be executed in parallel, and based on that, calculate duration.
        
        More than one multi-qubit gate can be executed in parallel iff both conditions are satisfied:
            - subsets of the qubits acted on by the gates are disjoint. 
            - blockade radii of operands of *different* gates do not overlap, i.e.,
               
        E.g., CZ(qa,qb) and CZ(qc,qd) can execute in parallel if {qa,qb} and {qc,qd} are disjoint and 
              if (qa,qc), (qa,qd), (qb,qc), (qb,qd) are *not* in the edge set of the connectivity graph.
        '''
        
        initial_set = set()
        initial_set.add(self.qubit_tuples[0])
        self.submoments = [initial_set]
        
        for q_tuple in self.qubit_tuples[1:]:
            not_scheduled = True
            i = 0
            
            while not_scheduled and i < len(self.submoments):
                
                curr_submoment = self.submoments[i]
                
                disjoint_qubit_subsets = all([[q not in curr_gate for q in q_tuple] 
                                              for curr_gate in curr_submoment])
                nonoverlapping_blockades = all([[(p[0],p[1]) not in connectivity_graph 
                                                for p in product(q_tuple,curr_gate)] 
                                                for curr_gate in curr_submoment])
                
                if disjoint_qubit_subsets and nonoverlapping_blockades:
                    self.submoments[i].add(q_tuple)
                    not_scheduled = False
            
                i+=1
                
            if not_scheduled:
                new_submoment = set()
                new_submoment.add()
                self.submoments.append(new_submoment)
                
        cz_duration = 0
        ccz_duration = 0
        for submoment in self.submoments:
            if any([len(q_tuple)==3 for q_tuple in submoment]):
                ccz_duration+=gate_durations['ccz']
            else:
                cz_duration+=gate_durations['cz']
                
        self.duration = cz_duration + ccz_duration 
        self.cz_duration = cz_duration
        self.ccz_duration = ccz_duration