import qiskit

def get_schedule_GR_cost(schedule):
    final_cost = 0
    for moment in schedule:
        if len(moment)!=0 and len(moment[0][1])==1:
            final_cost += max([float(op[0].params[0]) for op in moment])
    return round(final_cost,3)

def print_schedule(schedule, with_angles=False):
    i = 0
    for moment in schedule:
        tab = '     ' if i<10 else '    '
        string = '['
        if all([isinstance(op[0],qiskit.circuit.library.standard_gates.u3.U3Gate) for op in moment]):
            for op in moment:
                if with_angles:
                    theta, phi, lam = [round(float(x),2) for x in op[0].params]
                    string += f'U3({theta},{phi},{lam},q{op[1][0].index}), '
                else:
                    string += f'U3(q{op[1][0].index}), '
            #print(f'Moment {i}:'+ tab + string[:-2] + ']')
            print(string[:-2] + ']')
        else:
            for op in moment:
                string += f'CZ(q{op[1][0].index},q{op[1][1].index}), ' if len(op[1])==2 \
                          else f'CCZ(q{op[1][0].index},q{op[1][1].index},q{op[1][2].index}), '
            #print(f'Moment {i}:'+ tab + string[:-2] + ']')
            print(string[:-2] + ']')

        i+=1
        
def get_circuit_from_schedule_with_barriers(schedule,n):
    circuit = qiskit.QuantumCircuit(n)
    for moment in schedule:
        for op in moment:
            circuit.append(op[0],[q.index for q in op[1]])
        circuit.barrier([q for q in range(n)])
    return circuit

def get_circuit_from_schedule(schedule,n):
    circuit = qiskit.QuantumCircuit(n)
    for moment in schedule:
        for op in moment:
            circuit.append(op[0],[q.index for q in op[1]])
    return circuit