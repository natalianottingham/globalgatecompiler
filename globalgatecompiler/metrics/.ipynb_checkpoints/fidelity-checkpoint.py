from ..decomposition.gate_and_moment_classes import RzMoment, GRMoment, MultiQubitGateMoment
import numpy as np

def get_fidelity_cost_data(decomposed_moments,rz_constant,gr_constant,cz_fidelity,ccz_fidelity,t2_star,
                           rz_pi_duration=None,gr_pi_duration=None,cz_duration=None,ccz_duration=None,connectivity_graph=None):
    
    if decomposed_moments[0].duration==None:
        assert rz_pi_duration!=None
        assert gr_pi_duration!=None
        assert cz_duration!=None
        assert ccz_duration!=None
        assert connectivity_graph!=None
        
        gate_durations = {'rz': rz_pi_duration, 'gr': gr_pi_duration, 'cz': cz_duration, 'ccz': ccz_duration}
        for moment in decomposed_moments:
            moment.set_duration(gate_durations,connectivity_graph) 
    
    data = {'gate': 1, 'idle': 1}
    
    for moment in decomposed_moments:
        if isinstance(moment,RzMoment):
            for angle in moment.angles.values():
                data['gate'] *= 1 - rz_constant*abs(angle)    
        
        elif isinstance(moment,GRMoment):
            angle = moment.theta
            data['gate'] *= 1 - gr_constant*(angle**2)   

        elif isinstance(moment,MultiQubitGateMoment):
            for qubit_tuple in moment.qubit_tuples:
                if len(qubit_tuple)==2: # CZ gate
                    data['gate'] *= cz_fidelity
                else:
                    assert len(qubit_tuple)==3 # CCZ gate
                    data['gate'] *= ccz_fidelity  
                        
        data['idle'] *= np.e**(-moment.duration/t2_star)
        
    data['total'] = data['gate']*data['idle']
    for error_type in data:
        data[error_type] = float(round(data[error_type],3))
    return data