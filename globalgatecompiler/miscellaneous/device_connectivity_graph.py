import qiskit
from qiskit import QuantumCircuit
import networkx as nx
import math

def device_connectivity_graph(circuit,blockade_radius:int,spacing:int,grid_height:int=None,grid_width:int=None) -> nx.Graph:
    '''
    Returns graph where nodes are physical/hardware qubits, i.e. atom sites, and an edge exists between 
    two nodes if a two-qubit gate can be executed between the two corresponding atom sites (i.e., if the 
    atom sites are within a blockade radius of each other). Edge weights give distance between atoms.
    Atoms are assumed to be arranged in a 2D grid, but code can be easily adjusted for other geometries.
    
    Inputs:
        - blockade_radius: atoms must be within this distance of each other to execute a two-qubit gate.
                           (sometimes referred to as "Rydberg radius"). Units of micrometers.
        - spacing: distance between neighboring atom sites. Units of micrometers.
        - n: number of program qubits in the circuit.
        - grid_height,grid_width: number of rows and columns, respectively, in the 2D grid. 
                                  NOTE: need grid_height*grid_width>=n, i.e., number of hardware qubits must 
                                  be greater than or equal to the nubmer of program qubits in the circuit.
    Outputs:
        - cg: networkx Graph object defining the connectivity graph.
    '''
    
    n = len(circuit.qubits)
    
    assert blockade_radius >= spacing, 'To execute two-qubit gates, need blockade radius >= spacing' 
    if grid_width==None and grid_height==None:
        grid_width = math.ceil(math.sqrt(n))
        grid_height = math.ceil(n/grid_width) 
    elif grid_width is not None:
        grid_height = math.ceil(n/grid_width)
    elif grid_height is not None:
        grid_width = math.ceil(n/grid_height)
    else: # both grid_width and grid_height are defined
        assert grid_width*grid_height >= n, 'Number of hardware qubits must be >= number of program qubits'     
    
    num_phys_qubits = grid_height*grid_width
    cg = nx.Graph()
    cg.add_nodes_from([i for i in range(num_phys_qubits)]) # nodes = hardware qubits
    
    hardware_location = {} # defines row and column index of hardware qubit in atom array
    for i in range(grid_height):
        for j in range(grid_width):
            phys_qubit = i*grid_width+j
            hardware_location[phys_qubit] = (i,j)
    
    for k in range(num_phys_qubits):
        ki, kj = hardware_location[k]
        for m in range(k+1,num_phys_qubits):
            mi, mj = hardware_location[m]
            distance = math.sqrt((ki-mi)**2 + (kj-mj)**2)*spacing
            if distance <= blockade_radius:
                cg.add_edge(k,m,weight=distance)
                
    return cg