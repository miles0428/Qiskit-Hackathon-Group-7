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
from qiskit_algorithms.optimizers import SLSQP, COBYLA, SPSA, ADAM
import numpy as np 
from qiskit.circuit.library import ZZFeatureMap, RealAmplitudes, Initialize
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime.fake_provider import FakeKyoto, FakeNairobi
from qiskit.circuit.random import random_circuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Optimize1qGates, CommutationAnalysis,CommutativeCancellation,Collect2qBlocks
from qiskit.transpiler.passes.optimization.consolidate_blocks import ConsolidateBlocks
from qiskit import QuantumCircuit
from qiskit.transpiler.passes.synthesis.unitary_synthesis import UnitarySynthesis
import pyscf
import pyscf.cc
import pyscf.data.elements
import ffsim
from qiskit_aer import QasmSimulator,UnitarySimulator
import qiskit.visualization as visualization
import matplotlib.pyplot as plt

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
        for group in groups:
            three_q_block = []
            for index in group:
                assert type(two_q_blocks[index]) == list
                for node in two_q_blocks[index]:
                    three_q_block.append(node)
                
            # if len(two_q_blocks[index]) >20:
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
        norm_distance_record = []
        # index = len(centers)-1
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
                    if done:
                        break
                # print(node.qargs)
                matrix = node.matrix
                if matrix.shape != (8, 8):
                    # print(matrix.shape)
                    continue
                qc = QuantumCircuit(3)
                qc_target = QuantumCircuit(3)
                # print(matrix)
                qc_target.unitary(matrix, [0, 1, 2])
                target_unitary = UnitarySimulator().run(circuits=qc_target).result().get_unitary()
                state_vector_target = np.array(Statevector(qc_target)).reshape(-1)
                # state_vector_target = Statevector(qc_target)
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
                for i in range(3):
                    if i == center_index:
                        pass
                    else:
                        qc.ecr(i, center_index) 
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
                    state_vector = np.array(state_vector).reshape(-1)
                    backend = UnitarySimulator()
                    r = np.array(backend.run(circuits=qc_test).result().get_unitary())
                    # print(r)
                    norm_distance = np.linalg.norm(target_unitary - r)
                    # print(norm_distance)
                    return norm_distance
        
                optimizer = COBYLA(maxiter=600)
                optimizer2 = SPSA(maxiter=600)
                initial_point = np.random.rand(p_index)*2*np.pi
                # res = optimizer.minimize(cost_func, initial_point)
                # res = optimizer2.minimize(cost_func, res.x)
                # norm_distance_record.append(cost_func(res.x))
                qc = qc.assign_parameters(initial_point)
                # qc = qc.assign_parameters(res.x)

                # statevector = np.array(Statevector(qc)).reshape(-1)
                # print(f'fidelity: {np.abs(np.dot(np.conj(statevector), state_vector_target))}')

                minidag = DAGCircuit()
                Register = QuantumRegister(3)
                minidag.add_qreg(Register)
                minidag.apply_operation_back(qc.to_instruction(), [Register[i] for i in range(3)])
                # index -=1
                dag.substitute_node_with_dag(node, minidag)
            else :
                pass
        return dag


# class MLUnitarySynthesis2(TransformationPass):
#     def __init__(self, basis_gates=["ecr", 'id', 'rz', 'sx', 'x'], coupling_map=None):
#         super().__init__()
#         self.basis_gates = basis_gates
#         self.coupling_map = coupling_map
#     def run(self, dag: DAGCircuit) -> DAGCircuit:
#         for node in dag.op_nodes():
#             if node.op.name =="unitary":
#                 # print(index)
#                 center_index = 3
#                 # center_index = node.qargs.index(centers[index])
#                 for i,o in enumerate(node.qargs):
#                     done = False
#                     for j,p in enumerate(node.qargs):
#                         if i == j :
#                             pass
#                         else:
#                             indice = [o._index, p._index]
#                             indice2 = [p._index, o._index]
#                             # print(indice)
#                             if (indice not in self.coupling_map) and (indice2 not in self.coupling_map):
#                                 center_index -= i
#                                 center_index -= j
#                                 done = True
#                                 break
#                     if done:
#                         break
#                 # print(node.qargs)
#                 matrix = node.matrix
#                 if matrix.shape != (8, 8):
#                     print(matrix.shape)
#                     continue
#                 qc = QuantumCircuit(3)
#                 qc.unitary(matrix, [0, 1, 2])
#                 coupling_map = []
#                 for i in range(3):
#                     if i == center_index:
#                         pass
#                     else:
#                         coupling_map.append([i, center_index])
#                         coupling_map.append([center_index, i])

#                 pass_manager = generate_preset_pass_manager(basis_gates=['ecr', 'id', 'rz', 'sx', 'x'],optimization_level=3, coupling_map=coupling_map)
#                 pass_manager_c = PassManager()
#                 pass_manager_c.append(Collect2qBlocks())
#                 pass_manager_c.append(ConsolidateBlocks())
#                 pass_manager_c.append(UnitarySynthesis())
#                 pass_manager_c.append(Optimize1qGates())
#                 pass_manager_c.append(CommutativeCancellation())
#                 pass_manager_c.append(Optimize1qGates())
#                 pass_manager_c.append(CommutativeCancellation())
#                 qc = pass_manager.run(qc)
#                 qc = pass_manager_c.run(qc)

#                 print(qc.draw())
#             else :
#                 pass
#         return dag



if __name__ == "__main__":
    
    # qc = random_circuit(6, 3)
    # qc = QuantumCircuit(5)
    # qc.cx(2, 3)
    # qc.ry(2.2, 4)
    # qc.cx(3,4)
    # qc.ry(2.2, 3)
    # qc.cx(1, 2)
    # qc.x(0)
    # qc.ry(2.2, 1)
    # qc.cx(0, 1)
    benchmark_num_qubits = [i%20+5 for i in range(100)]
    benchmark_depth_default = []
    benchmark_depth_custom = []

    # state_vector_actual = Statevector(qc)
    backend = FakeKyoto()
    coupling_map = backend.configuration().coupling_map
    # print(coupling_map)
    pass_manager = generate_preset_pass_manager(basis_gates=['ecr', 'id', 'rz', 'sx', 'x'], backend=backend,optimization_level=3)
    pass_manager_custom = PassManager()
    pass_manager_custom.append(Collect_3q_blocks())
    pass_manager_custom.append(ConsolidateBlocks())
    pass_manager_custom_2 = PassManager()
    pass_manager_custom_2.append(MLUnitarySynthesis(coupling_map=coupling_map))
    pass_manager_custom_2.append(CommutationAnalysis())
    pass_manager_custom_2.append(CommutativeCancellation())
    pass_manager_custom_2.append(Optimize1qGates())

    for i in benchmark_num_qubits:
        qc = random_circuit(i, 10)
        qc = pass_manager.run(qc)
        benchmark_depth_default.append(qc.depth())
        qc = pass_manager_custom.run(qc)
        qc = pass_manager_custom_2.run(qc).decompose()
        qc = pass_manager_custom.run(qc)
        qc = pass_manager_custom_2.run(qc).decompose()
        benchmark_depth_custom.append(qc.depth())
    
    plt.scatter(benchmark_num_qubits, benchmark_depth_default, label='Qiskit default')
    plt.scatter(benchmark_num_qubits, benchmark_depth_custom, label='Collect_3q_blocks')
    plt.xlabel('number of qubits')
    plt.ylabel('depth')
    plt.legend()
    # plt.show()
    plt.savefig('benchmark_qubit.png')
    plt.show()
    plt.clf()
    plt.scatter(benchmark_num_qubits, np.array(benchmark_depth_default)-np.array(benchmark_depth_custom), )
    plt.xlabel('number of qubits')
    plt.ylabel('depth')
    plt.title('difference between default and custom')
    plt.show()
    plt.savefig('benchmark_qubit_diff.png')



    # transpiled_qc = pass_manager.run(qc)
    # # state_vector_actual = Statevector(transpiled_qc)

    # print(f'Depth: {transpiled_qc.depth()}')

    # # print(transpiled_qc.draw(idle_wires=False))


    # transpiled_qc = pass_manager_custom.run(transpiled_qc)
    # # state_vector_actual = Statevector(transpiled_qc)
    # print(transpiled_qc.draw(idle_wires=False))

    # transpiled_qc = pass_manager_custom_2.run(transpiled_qc).decompose()
    # print(f'Depth: {transpiled_qc.depth()}')
    # print(transpiled_qc.draw(idle_wires=False))
    # statevector = Statevector(transpiled_qc)

    # fidelity = np.abs(np.dot(np.conj(np.array(statevector).reshape(-1)), np.array(state_vector_actual).reshape(-1)))
    # print(f'fidelity: {fidelity}')



