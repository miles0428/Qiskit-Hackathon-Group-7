from qiskit.transpiler.basepasses import TransformationPass,AnalysisPass
from qiskit.dagcircuit import DAGCircuit
from qiskit.circuit.library import XGate, YGate
from qiskit import transpiler
from qiskit.dagcircuit.dagnode import DAGNode, DAGOpNode, DAGInNode, DAGOutNode
from qiskit.circuit.gate import Gate
from qiskit.circuit.quantumregister import QuantumRegister, Qubit
from collections import defaultdict
from qiskit.transpiler.passes.synthesis.unitary_synthesis import UnitarySynthesis
from qiskit.circuit import Parameter
from qiskit.quantum_info import Operator, Statevector, DensityMatrix
from qiskit_algorithms.optimizers import COBYLA
import numpy as np 
from qiskit.circuit.library import ZZFeatureMap, RealAmplitudes, Initialize
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime.fake_provider import FakeKyoto, FakeNairobi
from qiskit.circuit.random import random_circuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Optimize1qGates
from qiskit.transpiler.passes.optimization.consolidate_blocks import ConsolidateBlocks





class Collect_3q_blocks(AnalysisPass):
    def run(self, dag: DAGCircuit) -> DAGCircuit:
        three_q_blocks , centers= self.collect_3q_blocks(dag)
        self.property_set["commutation_set"] = defaultdict(list)
        self.property_set["block_list"] = three_q_blocks
        self.property_set["centers"] = centers
        # print(centers)
        # dag.collect_2q_runs()
        return dag
    
    def collect_3q_blocks(self,dag: DAGCircuit):
        two_q_blocks = dag.collect_2q_runs()
        q_args_list=[]
        
        for i in two_q_blocks:
            counter = 0
            q_args = [] 
            for node in i :
                for qubit in node.qargs:
                    if qubit not in q_args:
                        q_args.append(qubit)
                        counter += 1
                    if counter == 2:
                        break
                if counter == 2:
                    break
            q_args_list.append(q_args)

        # print(q_args_list)
        group = []
        groups = []
        group_qubit = []
        centers = []
        for i,q_args in enumerate(q_args_list):
            Pass = True
            for elem in q_args:
                if elem not in group_qubit:
                    if len(group_qubit)>=3:
                        Pass = False
                        break
                    else :
                        group_qubit.append(elem)
                else:
                    Pass = True
            if Pass :
                group.append(i)                    
            else:
                groups.append(group)             
                group = [i]
                group_qubit = [i for i in q_args]
        groups.append(group)
        three_q_blocks =[]
        # pack each group into a block
        # print(groups)
        for group in groups:
            three_q_block = []
            for index in group:
                assert type(two_q_blocks[index]) == list
                for node in two_q_blocks[index]:
                    three_q_block.append(node)
            three_q_blocks.append(three_q_block)
        
        return three_q_blocks, centers

class MLUnitarySynthesis(TransformationPass):
    def __init__(self, basis_gates=["ecr", 'id', 'rz', 'sx', 'x'], coupling_map=None):
        super().__init__()
        self.basis_gates = basis_gates
        self.coupling_map = coupling_map
    def run(self, dag: DAGCircuit) -> DAGCircuit:
        # print(self.property_set["commutation_set"])
        centers = self.property_set["centers"]
        nodes = dag.op_nodes()
        index = len(centers)-1
        for node in dag.op_nodes():
            if node.op.name =="unitary":
                # print(index)
                center_index = 3
                # center_index = node.qargs.index(centers[index])
                for i,o in enumerate(node.qargs):
                    done = False
                    for j,p in enumerate(node.qargs):
                        if i == j :
                            pass
                        else:
                            indice = [o._index, p._index]
                            indice2 = [p._index, o._index]
                            # print(indice)
                            if (indice not in self.coupling_map) and (indice2 not in self.coupling_map):
                                center_index -= i
                                center_index -= j
                                done = True
                                break
                            # else:
                            #     print("hello")
                    if done:
                        break
                # print(node.qargs)
                matrix = node.matrix
                qc = QuantumCircuit(3)
                qc_target = QuantumCircuit(3)
                # print(matrix)
                qc_target.unitary(matrix, [0, 1, 2])
                state_vector_target = Statevector(qc_target)
                p_index = 0
                for i in range(3):
                    qc.sx(i)
                    qc.rz(Parameter(f'p{p_index}'), i)
                    p_index += 1
                    qc.sx(i)
                    qc.rz(Parameter(f'p{p_index}'), i)
                    p_index += 1
                for i in range(3):
                    if i == center_index:
                        pass
                    else:
                        qc.ecr(center_index, i) 
                for i in range(3):
                    qc.sx(i)
                    qc.rz(Parameter(f'p{p_index}'), i)
                    p_index += 1
                    qc.sx(i)
                    qc.rz(Parameter(f'p{p_index}'), i)
                    p_index += 1

                def cost_func(params):
                    qc_test = qc.assign_parameters(params)
                    state_vector = Statevector(qc_test)
                    fidelity = state_vector.inner(state_vector_target)
                    return 1- fidelity
                
                optimizer = COBYLA(maxiter=500)
                initial_point = np.random.rand(p_index)*2*np.pi
                res = optimizer.minimize(cost_func, initial_point)
                qc = qc.assign_parameters(res.x)

                minidag = DAGCircuit()
                Register = QuantumRegister(3)
                minidag.add_qreg(Register)
                minidag.apply_operation_back(qc.to_instruction(), [Register[i] for i in range(3)])
                index -=1
                dag.substitute_node_with_dag(node, minidag)
            else :
                pass
        return dag

from qiskit import QuantumCircuit


# qc = ZZFeatureMap(5,2)
qc = random_circuit(5,20)
# state_vector_actual = Statevector(qc)
backend = FakeNairobi()
coupling_map = backend.configuration().coupling_map
# print(coupling_map)
pass_manager = generate_preset_pass_manager(basis_gates=['ecr', 'id', 'rz', 'sx', 'x'], backend=backend,optimization_level=3)






pass_manager_custom = PassManager()
pass_manager_custom.append(Collect_3q_blocks())
pass_manager_custom.append(ConsolidateBlocks())
pass_manager_custom.append(MLUnitarySynthesis(coupling_map=coupling_map))

transpiled_qc = pass_manager.run(qc)

print(f'Depth: {transpiled_qc.depth()}')

# print(transpiled_qc.draw(idle_wires=False))


transpiled_qc = pass_manager_custom.run(transpiled_qc).decompose()
state_vector_actual = Statevector(transpiled_qc)
# print(transpiled_qc.draw(idle_wires=False))
print(f'Depth: {transpiled_qc.depth()}')

state_vector_transpiled = Statevector(transpiled_qc)

print(f'fidelity: {state_vector_actual.inner(state_vector_transpiled)}')


