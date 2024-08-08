import qiskit
import numpy as np
import math
from .gate_and_moment_classes import *

def convert_angle(angle):
    '''
    If abs(angle)>pi, convert to equivalent rotation amount with abs value <= pi. 
    '''
    angle = angle % (2*np.pi) # gives value >= 0
    if angle > np.pi:
        angle = angle - 2*np.pi
    assert abs(angle) <= np.pi
    return angle

def get_axial_decomposition_angles(moment,circuit,eta=0):
    '''
    For a given moment in the circuit containing only u3 gates, return the Rz and GR angles
    necessary to decompose into the neutral atom gate set using the axial decomposition method.
    Each moment is decomposed into 3 columns of Rz gates, with columns separated by GR gates.
    Rz gates appear on qubits that had u3 gates in the original moment; GRs act on all qubits.
    Eta can be changed to optimize Rz rotation amount in the final circuit.
    
    Input: 
        - moment: list of CircuitInstruction objects, each of which specifies a u3 gate;
                  we assume that there is at most one gate on each qubit in the input moment.
    Output:
        - dict specifying Rz and GR angles in the decomposition:
            angles_all['rz']['first'][q]: rotation angle for Rz gate in first column on qubit q
                                          (same idea for 'middle' and 'last' columns of Rz gates)
            angles_all['gr']['first']: theta, phi for GR gate between first and middle Rz columns
            angles_all['gr']['last']: theta, phi for GR gate between middle and last Rz columnns
    '''
    
    assert all([isinstance(op[0],qiskit.circuit.library.standard_gates.u3.U3Gate) or
                isinstance(op[0],qiskit.circuit.library.UGate) for
                op in moment]), 'Input should be a moment of only u3 gates'
    
    angles_all = {}
    angles_all['rz'] = {'first': {}, 'middle': {}, 'last': {}}
    angles_all['gr'] = {'first': {'theta': np.pi/2, 'phi': eta}, 
                    'last': {'theta': -np.pi/2, 'phi': eta}} 

    for op in moment:
        euler_angles = op[0].params
        theta_ = float(euler_angles[0])
        phi_ = float(euler_angles[1])
        lambda_ = float(euler_angles[2])
        rz_angles_this_op = {'first': lambda_ + eta, 'middle': theta_, 'last': phi_ - eta}

        qubit = circuit.find_bit(op[1][0]).index
        assert all([qubit not in angles_all['rz'][x] for x in 
                    ['first','middle','last']]), 'Maximum one u3 gate on each qubit per moment'
        for x in rz_angles_this_op:
            angles_all['rz'][x][qubit] = rz_angles_this_op[x] 
                
    return angles_all

def get_transverse_decomposition_angles(moment,circuit,eta=0,sign_theta_max=1,sign_sigma_j=1):
    '''
    Similar to get_axial_decomposition_angles, except applying the transverse decomposition method. 
    For a given moment of u3 gates, the transverse decomposition minimizes GR rotation amount.
    Eta, sign_theta_max, sigma_j can be changed to optimize Rz rotation amount in the final circuit.
    '''
    
    assert sign_sigma_j in [-1,1] and sign_theta_max in [-1,1]
    assert all([isinstance(op[0],qiskit.circuit.library.standard_gates.u3.U3Gate) or
                isinstance(op[0],qiskit.circuit.library.UGate) for
                op in moment]), 'Input should be a moment of only u3 gates'
    
    theta_max = max([abs(float(op[0].params[0])) for op in moment])
    
    angles_all = {}
    angles_all['rz'] = {'first': {}, 'middle': {}, 'last': {}}
    angles_all['gr'] = {'first': {'theta': -theta_max/2, 'phi': np.pi/2 + eta}, 
                        'last': {'theta': theta_max/2, 'phi': np.pi/2 + eta}}

    for op in moment:
        euler_angles = op[0].params
        theta = float(euler_angles[0])
        phi_minus = float(euler_angles[1])
        phi_plus = float(euler_angles[2])
        
        if abs(round(theta,5))!=abs(round(theta_max,5)):
            kappa = math.sqrt((math.sin(theta/2)**2)/(math.sin(theta_max/2)**2-math.sin(theta/2)**2))
            alpha = math.atan(math.cos(theta_max/2)*kappa)
            chi = -sign_sigma_j*2*math.atan(kappa)
        else: # |theta|=|theta_max| --> kappa = inf, arctan(inf) = pi/2
            alpha = np.pi/2
            chi = -sign_sigma_j*np.pi
        beta = theta/abs(theta)*sign_theta_max*np.pi/2 if theta!=0 else 0
    
        gamma_plus = phi_plus+sign_sigma_j*(alpha+beta)
        gamma_minus = phi_minus+sign_sigma_j*(alpha-beta)
        
        rz_angles_this_op = {'first': gamma_plus+eta, 'middle': chi, 'last': gamma_minus-eta}
        
        qubit = circuit.find_bit(op[1][0]).index
        assert all([qubit not in angles_all['rz'][x] for x in 
                    ['first','middle','last']]), 'Maximum one u3 gate on each qubit per moment'
        for x in rz_angles_this_op:
            angles_all['rz'][x][qubit] = rz_angles_this_op[x]   
                
    return angles_all

def decompose_to_neutral_atom_gate_set(schedule,decomposition_type='transverse',use_backlog=True,
                                       eta=0,sign_theta_max=1,sign_sigma_j=1):
    '''
    Given a circuit in terms of the gate set {u3,CZ,CCZ}, decompose to the gate set {Rz,GR,CZ,CCZ}.
    
    Inputs:
        - schedule: list of moments, in the order they appear in the circuit.
                    Each moment is a list of CircuitInstruction objects and contains only u3 or only CZ/CCZ gates
                    (i.e., a CZ/CCZ gate and u3 gate cannot appear within the same moment). 
        - decomposition_type: 'axial' or 'transverse'
        - use_backlog: if True, for each u3 moment, the last column of Rz gates is combined with the first column 
                       of Rz gates in the next u3 moment. This reduces total Rz gate cost in the circuit. 
        - eta: can be any number in [0,2pi); used to optimize Rz costs.
        - sign_theta_max: either +1 or -1; used to optimize Rz costs. Only relevant with transverse decomposition.
        - sign_sigma_j: either +1 or -1; used to optimize Rz costs. Only relevant with transverse decomposition.
    
    Outputs: 
        - c: the decomposed circuit.
        - decomposed_moments: list of RzMoment, GRMoment, and MultiQubitGateMoment objects, in the order they
                              appear in the circuit. This is used as input when calculating fidelity and duration.
    '''
    
    assert decomposition_type in ['axial','transverse'], 'Supported decomposition types: axial, transverse'

    n = len(schedule.circuit.qubits)
    qr = qiskit.QuantumRegister(n)
    c = qiskit.QuantumCircuit(qr)
    
    decomposed_moments = []
    
    if use_backlog:
        backlog = np.zeros(n)
        
    for moment in schedule.schedule:
        if len(moment)==0:
            continue

        if all([len(op[1])>=2 for op in moment]): # Multi-Qubit Gate Moment (MQGM) --> keep same
            assert all([isinstance(op[0],qiskit.circuit.library.standard_gates.z.CZGate) or
                        isinstance(op[0],qiskit.circuit.library.standard_gates.z.CCZGate) 
                        for op in moment]), 'Input circuit should be in gate set {u3,cz,ccz}'
            for op in moment:
                c.append(op[0],qargs=[qr[schedule.circuit.find_bit(qubit).index] for qubit in op[1]])
            decomposed_moments.append(MultiQubitGateMoment([tuple(schedule.circuit.find_bit(qubit).index 
                                                                  for qubit in op[1]) for op in moment]))

        else: # Single-Qubit Gate Moment (SQGM) --> need to decompose
            assert all([len(op[1])==1 for op in moment]),'Cannot have 1-qubit and 2-qubit gates in same moment'

            # get angles to use for Rz and GR gates in the decomposition
            if decomposition_type=='axial':
                decomposition_angles = get_axial_decomposition_angles(moment,schedule.circuit,eta=eta)
            else:
                assert decomposition_type=='transverse', 'Supported decomposition types: axial, transverse'
                decomposition_angles = get_transverse_decomposition_angles(moment,schedule.circuit,eta=eta,
                                                                           sign_theta_max=sign_theta_max,
                                                                           sign_sigma_j=sign_sigma_j)
                if round(decomposition_angles['gr']['first']['theta'],5)==0:
                    for qubit in decomposition_angles['rz']['first']:
                        backlog[qubit] += sum([decomposition_angles['rz'][x][qubit] 
                                               for x in ['first','middle','last']])
                    continue

            # add gates to the circuit based on the angles from above
            moment_angles = {}
            for qubit in decomposition_angles['rz']['first']:
                if use_backlog:
                    angle = decomposition_angles['rz']['first'][qubit] + backlog[qubit]
                else:
                    angle = decomposition_angles['rz']['first'][qubit]
                if abs(angle) > np.pi:
                    angle = convert_angle(angle)
                if round(angle,5)!=0:
                    c.rz(angle,qr[qubit])
                    moment_angles[qubit] = angle
            if len(moment_angles) > 0:
                decomposed_moments.append(RzMoment(moment_angles))
            
            theta = decomposition_angles['gr']['first']['theta']
            phi = decomposition_angles['gr']['first']['phi']
            c.append(GRGate(n,theta,phi),qr)
            decomposed_moments.append(GRMoment(theta,phi,n))
                
            assert len(decomposition_angles['rz']['middle'])!=0 # if len==0, GR gates cancel
            moment_angles = {}
            for qubit in decomposition_angles['rz']['middle']:
                angle = decomposition_angles['rz']['middle'][qubit]
                if abs(angle) > np.pi:
                    angle = convert_angle(angle)
                if round(angle,5)!=0:
                    c.rz(angle,qr[qubit])
                    moment_angles[qubit] = angle
            if len(moment_angles) > 0:
                decomposed_moments.append(RzMoment(moment_angles))

            theta = decomposition_angles['gr']['last']['theta']
            phi = decomposition_angles['gr']['last']['phi']
            c.append(GRGate(n,theta,phi),qr)
            decomposed_moments.append(GRMoment(theta,phi,n))
                
            if use_backlog:
                for qubit in decomposition_angles['rz']['last']:
                    backlog[qubit] = decomposition_angles['rz']['last'][qubit]
            else:
                moment_angles = {}
                for qubit in decomposition_angles['rz']['last']:
                    angle = decomposition_angles['rz']['last'][qubit]
                    if abs(angle) > np.pi:
                        angle = convert_angle(angle)
                    if round(angle,5)!=0:
                        c.rz(angle,qr[qubit])
                        moment_angles[qubit] = angle
                if len(moment_angles) > 0:
                    decomposed_moments.append(RzMoment(moment_angles))
                    
    # clear remaining backlog
    if use_backlog:
        for q in range(len(backlog)):
            if backlog[q]!=0:
                c.rz(backlog[q],q)
            
    return c, decomposed_moments