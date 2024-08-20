from qiskit.transpiler.basepasses import TransformationPass,AnalysisPass
from qiskit.dagcircuit import DAGCircuit
from qiskit.circuit.library import XGate, YGate
from qiskit import transpiler
from qiskit.dagcircuit.dagnode import DAGNode, DAGOpNode, DAGInNode, DAGOutNode
from qiskit.circuit.gate import Gate
from qiskit.circuit.quantumregister import QuantumRegister, Qubit
from collections import defaultdict
from qiskit.transpiler.passes.synthesis.unitary_synthesis import UnitarySynthesis




class Collect_3q_blocks(AnalysisPass):
    def run(self, dag: DAGCircuit) -> DAGCircuit:
        three_q_blocks , centers= self.collect_3q_blocks(dag)
        self.property_set["commutation_set"] = defaultdict(list)
        self.property_set["block_list"] = three_q_blocks
        self.property_set["centers"] = centers
        dag.collect_2q_runs()
        return dag
    
    def collect_3q_blocks(self,dag: DAGCircuit):
        two_q_blocks = dag.collect_2q_runs()
        print(type(two_q_blocks))
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
                    else :
                        group_qubit.append(elem)
                else:
                    Pass = True
            if Pass :
                group.append(i)                    
            else:
                groups.append(group)
                # find the center
                for elem in group_qubit:
                    is_center = True
                    for index in group: 
                        if elem in two_q_blocks[index][0].qargs or elem in two_q_blocks[index][1].qargs:
                            is_center = True
                        else:
                            is_center = False
                            break
                    if is_center:
                        centers.append(elem)
                        break                      
                group = [i]
                group_qubit = [i for i in q_args]
        groups.append(group)
        for elem in group_qubit:
            is_center = True
            for index in group: 
                if elem in two_q_blocks[index][0].qargs or elem in two_q_blocks[index][1].qargs:
                    is_center = True
                else:
                    is_center = False
                    break
            if is_center:
                centers.append(elem)
                break  
        print(groups)
        three_q_blocks =[]
        # pack each group into a block
        for group in groups:
            three_q_block = []
            for index in group:
                for node in two_q_blocks[index]:
                    three_q_block.append(node)
            three_q_blocks.append(three_q_block)
        
        return three_q_blocks, centers

    
    # def collect_2q_runs(self, dag: DAGCircuit):
    #     """Return a set of non-conditional runs of 3q "op" nodes."""

    #     to_qid = {}
    #     for i, qubit in enumerate(dag.qubits):
    #         to_qid[qubit] = i

    #     def filter_fn(node):
    #         if isinstance(node, DAGOpNode):
    #             return (
    #                 isinstance(node.op, Gate)
    #                 and len(node.qargs) <= 2
    #                 and not getattr(node.op, "condition", None)
    #                 and not node.op.is_parameterized()
    #             )
    #         else:
    #             return None

    #     def color_fn(edge):
    #         if isinstance(edge, Qubit):
    #             return to_qid[edge]
    #         else:
    #             return None

    #     return rx.collect_bicolor_runs(dag._multi_graph, filter_fn, color_fn)


from qiskit import QuantumCircuit

qc = QuantumCircuit(4)
qc.h(0)
qc.cx(0, 1)
qc.h(1)
qc.cx(1, 2)
qc.h(2)
qc.h(0)
qc.cx(0, 1)
qc.h(1)
qc.cx(1, 2)
qc.h(2)
qc.cx(2, 3)
qc.cx(2, 1)
qc.measure_all()
print(qc.draw())


from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Optimize1qGates
from qiskit.transpiler.passes.optimization.consolidate_blocks import ConsolidateBlocks


pass_manager = PassManager()
pass_manager.append(Collect_3q_blocks())
pass_manager.append(ConsolidateBlocks())

transpiled_qc = pass_manager.run(qc)
print(transpiled_qc)
