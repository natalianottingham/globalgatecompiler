from ..decomposition.gate_and_moment_classes import RzMoment, GRMoment, MultiQubitGateMoment

def get_time_cost_data(decomposed_moments,rz_pi_duration,gr_pi_duration,cz_duration,ccz_duration,connectivity_graph=None):
    '''
    Inputs: 
        - decomposed_moments: list of moments in the circuit, in the order that they appear in the circuit; 
                              each moment is an object of class RzMoment, GRMoment, or MultiQubitGateMoment.
        - rz_pi_duration: time required to execute an Rz gate with rotation angle of pi. 
                          (duration of each Rz gate in the circuit is scaled by rotation amount)
        - gr_pi_duration: time required to execute an GR gate with rotation angle of pi.
                          (duration of each GR gate in the circuit is scaled by rotation amount)
        - cz_duration: time required to execute a CZ gate.
        - ccz_duration: time required to execute a CCZ gate. 
                          
    Output: 
        - data: Dict giving the amount of time spent executing gates in each category ('rz','gr','multi')
                over the entire circuit that was specified by the input moments. Overall duration for the
                entire circuit is given by data['total']. The category 'multi' includes CZ and CCZ gates.
    '''
    
    if decomposed_moments[0].duration==None:
        assert connectivity_graph!=None
        
        gate_durations = {'rz': rz_pi_duration, 'gr': gr_pi_duration, 'cz': cz_duration, 'ccz': ccz_duration}
        for moment in decomposed_moments:
            moment.set_duration(gate_durations,connectivity_graph) 
    
    data = {'rz': 0, 'gr': 0, 'multi': 0}
    
    for moment in decomposed_moments:
        if isinstance(moment,RzMoment):
            data['rz'] += moment.duration

        elif isinstance(moment,GRMoment):
            data['gr'] += moment.duration

        elif isinstance(moment,MultiQubitGateMoment):
            data['multi'] += moment.cz_duration + moment.ccz_duration
            
    data['total'] = data['rz'] + data['gr'] + data['multi']               
    for gate_type in data:
        data[gate_type] = float(round(data[gate_type],3))
    return data