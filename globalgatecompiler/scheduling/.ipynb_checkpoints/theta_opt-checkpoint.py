import qiskit
import itertools
import networkx as nx

from .sift import sift, circuit_to_dag
from .schedule_class import Schedule

def get_next_dp_args_list(c,v_c0,v_p1,v_c1,v_rem,dag_node_to_gate,check_last_cond):
        
    # initializations
    v_p1_star = [node for node in v_p1]
    v_p1_square = []
    v_c1_star = [node for node in v_c1]
    v_c1_square = []
        
    v_c0_cost = max([float(dag_node_to_gate[node][0].params[0]) for node in v_c0])
    v_c1_cost = max([float(dag_node_to_gate[node][0].params[0]) for node in v_c1]) if v_c1 else 0
    orig_vc0_plus_vc1_cost = v_c0_cost + v_c1_cost
    new_vc1_cost = max(v_c0_cost,v_c1_cost)
    
    args = [(v_c0,v_p1,v_c1,v_rem,v_c0_cost)] # args for case where M_k=v_c0 (always considered)


    v_c0.sort(key=lambda x: float(dag_node_to_gate[x][0].params[0]), reverse=True)
    potential_arg_info = []
    prev_theta_value = round(abs(float(dag_node_to_gate[v_c0[0]][0].params[0])),5)
    mk_qubits_just_added = [c.find_bit(dag_node_to_gate[v_c0[0]][1][0]).index]
    for i in range(1,len(v_c0)):
        this_op_theta_value = round(float(dag_node_to_gate[v_c0[i]][0].params[0]),5)
        if this_op_theta_value!=prev_theta_value:
            potential_arg_info.append((v_c0[i:],v_c0[:i],mk_qubits_just_added,this_op_theta_value))
            prev_theta_value = this_op_theta_value
            mk_qubits_just_added = []
        mk_qubits_just_added.append(c.find_bit(dag_node_to_gate[v_c0[i]][1][0]).index)

    
    for M_k,m_k,mk_qubits_just_added,max_theta in potential_arg_info: 
    # iterate through each unique theta value (this guarantees condition 1 is satisfied)
    # starts with largest theta (sorted in reverse order above)
            
        # check condition 2 (MQGM in between is non-empty)
        nodes_to_push_back = []
        vp1_square_qubits_just_added = set()
        for node in v_p1_star:
            node_qubits = set([c.find_bit(q).index for q in dag_node_to_gate[node][1]])
            if not node_qubits.isdisjoint(mk_qubits_just_added):
                nodes_to_push_back.append(node)
                vp1_square_qubits_just_added.update(node_qubits)
        for node in nodes_to_push_back:
            v_p1_star.remove(node)
            v_p1_square.append(node)
            
        if len(v_p1_star)==0:
            break # if fails for this k, will also fail for all later k
            
        # check condition 3 (m_k doesn't make up its own SQGM)
        nodes_to_push_back = []
        for node in v_c1_star:
            node_qubits = set([c.find_bit(q).index for q in dag_node_to_gate[node][1]])
            if not node_qubits.isdisjoint(vp1_square_qubits_just_added.union(mk_qubits_just_added)):
                nodes_to_push_back.append(node)
        for node in nodes_to_push_back:
            v_c1_star.remove(node)
            v_c1_square.append(node)
            
        if len(v_c1_star)==0:
            break # if fails for this k, will also fail for all later k
            
        # check condition 4 (prev cost vs updated cost) if applicable
        if (not check_last_cond) or (check_last_cond and max_theta+new_vc1_cost<orig_vc0_plus_vc1_cost):
            args.append((M_k,
                        [node for node in v_p1_star],
                        m_k+[node for node in v_c1_star],
                        [node for node in v_p1_square]+[node for node in v_c1_square]+v_rem,
                        max_theta))

    return args

def dp(c,v_p0,v_c0,v_rem,dag_node_to_gate,n,subproblem_graph,subproblem_node_to_moments,
       prev_node,node_count,memo,check_last_cond):

    # check if solution already computed
    id_ = ','.join(str(i) for i in v_p0)+';'+','.join(str(i) for i in v_c0)+';'+','.join(str(i) for i in v_rem)
    if id_ in memo:
        cost,subproblem_node,edge_weight = memo[id_]
        subproblem_graph.add_edge(prev_node,subproblem_node,weight=edge_weight)
        return cost
    
    # base case
    if len(v_rem)==0:
        base_case_cost = float(dag_node_to_gate[v_c0[0]][0].params[0]) if v_c0 else 0 # max theta in v_c0 (v_c0 sorted in reverse)
        _ = add_node_to_subproblem_graph(v_p0,v_c0,base_case_cost,
                                         subproblem_graph,subproblem_node_to_moments,prev_node,node_count)
        return base_case_cost
    
    # recursive/dp case
    v_p1,v_c1,v_rem = sift(v_rem,dag_node_to_gate,n)
    args = get_next_dp_args_list(c,v_c0,v_p1,v_c1,v_rem,dag_node_to_gate,check_last_cond)
    cost_list = []
    for new_v_c0,new_v_p1,new_v_c1,new_v_rem,new_vc0_cost in args:
        new_node = add_node_to_subproblem_graph(v_p0,new_v_c0,new_vc0_cost,
                                                subproblem_graph,subproblem_node_to_moments,prev_node,node_count)
        cost_list.append(dp(c,new_v_p1,new_v_c1,new_v_rem,dag_node_to_gate,n,subproblem_graph,
                            subproblem_node_to_moments,new_node,node_count,memo,check_last_cond) + new_vc0_cost)   
    cost = min(cost_list)
    memo[id_] = (cost, new_node, new_vc0_cost)
    return cost

def add_node_to_subproblem_graph(vp,vc,vc_cost,subproblem_graph,subproblem_node_to_moments,prev_node,node_count):
    new_node = node_count[0]
    node_count[0] += 1
    subproblem_graph.add_node(new_node)
    subproblem_graph.add_edge(prev_node, new_node, weight=vc_cost)
    subproblem_node_to_moments[new_node] = [vp,vc]
    return new_node

def build_schedule_from_subproblem_graph(subproblem_graph,subproblem_node_to_moments,
                                         dag_node_to_gate,return_node_label_schedule):
    
    # find path in subproblem graph with the lowest weight (i.e. path with best cost)
    leaf_nodes = [node for node in subproblem_graph.nodes() if subproblem_graph.out_degree(node)==0]
    shortest_path_to_each_leaf = [nx.shortest_path(subproblem_graph,source=-1,target=leaf_node,weight='weight') 
                                  for leaf_node in leaf_nodes]
    best_path = min(shortest_path_to_each_leaf, key=lambda x:nx.path_weight(subproblem_graph,x,'weight'))
    
    # build schedule by traversing the path with best cost
    schedule = []
    for node in best_path[1:]:
        v_passed, v_caught = subproblem_node_to_moments[node]
        if return_node_label_schedule:
            schedule.append(v_passed)
            schedule.append(v_caught)
        else:
            schedule.append([dag_node_to_gate[label] for label in v_passed])
            schedule.append([dag_node_to_gate[label] for label in v_caught])
        
    if schedule[-1] == []: # if last iteration of DP had empty v_caught
        schedule = schedule[:-1]
    if schedule[0] == []: # if first iteration of sift returned empty v_passed
        schedule = schedule[1:]
        
    return schedule

def get_theta_opt_schedule(circuit,check_last_cond=False,return_node_label_schedule=False):
    n = len(circuit.qubits)
    
    dag, dag_node_to_gate = circuit_to_dag(circuit)
    g = list(nx.topological_sort(dag)) # list of all gates in topological order
    v_p0,v_c0,v_rem = sift(g,dag_node_to_gate,n)
    
    subproblem_graph = nx.DiGraph()
    subproblem_graph.add_node(-1)
    subproblem_node_to_moments = {}
    
    _ = dp(circuit,v_p0,v_c0,v_rem,dag_node_to_gate,n,subproblem_graph,subproblem_node_to_moments,-1,[0],{},check_last_cond)
    schedule = build_schedule_from_subproblem_graph(subproblem_graph,subproblem_node_to_moments,dag_node_to_gate,
                                                    return_node_label_schedule=return_node_label_schedule)
    return Schedule(schedule,circuit)