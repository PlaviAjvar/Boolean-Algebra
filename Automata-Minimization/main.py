import os
import collections

# tests if all partitions are reducible
# in other words if the superstates we have formed are all valid and we can reduce the automaton
def proper_partitions(partitions, transition_labels):
    for label, partition in partitions.items():
        partition_transit_labels = [transition_labels[state] for state in partition]
        # check if all elements are equal
        if partition_transit_labels[1:] != partition_transit_labels[:-1]:
            return False
    return True

# construct reduced automaton from properly partitioned automaton
# works for mealy and moore
def reduced_automaton(partitions, transition_labels, num_transitions, is_mealy):
    new_state_labels = {}
    new_automaton = []
    new_outputs = []
    state_counter = 0

    for label, partition in partitions.items():
        new_state_labels[label] = state_counter
        state_counter += 1

    for label, partition in partitions.items():
        representative = partition[0]
        new_automaton.append([new_state_labels[transition_label] for transition_label in transition_labels[representative]])
        if is_mealy:
            new_outputs.append(label[:num_transitions])
        else:
            new_outputs.append([label[0]])

    return new_automaton, new_outputs

# function that implements minimization for Mealy automaton
def minimize_automaton(automaton, outputs, is_mealy):
    num_transitions = len(automaton[0])
    # labels of states and inverse dictionary label -> state
    state_label = [tuple(output) for output in outputs]
    partitions = collections.defaultdict(list)
    for state, label in enumerate(state_label):
        partitions[label].append(state)
    # labels of neighbors of a state, in order (of output indices)
    transition_labels = [[state_label[neighbor] for neighbor in state] for state in automaton]
    while True:
        if proper_partitions(partitions, transition_labels):
            return reduced_automaton(partitions, transition_labels, num_transitions, is_mealy)

        # update state labels
        for label, partition in partitions.items():
            for state in partition:
                state_label[state] = state_label[state] + tuple(transition_labels[state])

        # update other two tables
        partitions = collections.defaultdict(list)
        for state, label in enumerate(state_label):
            partitions[label].append(state)
        transition_labels = [[state_label[neighbor] for neighbor in state] for state in automaton]

# helper function for inputing Moore automaton
def input_moore(lines, num_states, num_transitions):
    # empty lines list corresponds to standard input
    if lines == None:
        lines = [list(map(int, input().split())) for _ in range(num_states * 2)]
    else:
        lines = [list(map(int, line.split())) for line in lines]
    # input all state transitions
    automaton = lines[:num_states]
    for state in automaton:
        if len(state) != num_transitions:
            raise Exception("Invalid state transitions.")
    # input all transition outputs
    outputs = lines[num_states:]
    for state in outputs:
        if len(state) != 1:
            raise Exception("Invalid state output.")
    return automaton, outputs

# helper function for inputing Mealy automaton
def input_mealy(lines, num_states, num_transitions):
    # empty lines list corresponds to standard input
    if lines == None:
        lines = [list(map(int, input().split())) for _ in range(num_states*2)]
    else:
        lines = [list(map(int, line.split())) for line in lines]
    # input all state transitions
    automaton = lines[:num_states]
    for state in automaton:
        if len(state) != num_transitions:
            raise Exception("Invalid state transitions.")
    # input all transition outputs
    outputs = lines[num_states:]
    for state in outputs:
        if len(state) != num_transitions:
            raise Exception("Invalid transition outputs.")
    return automaton, outputs

# the input is a Mealy or Moore automaton
# the states are asumed to be zero indexed and inputed in order
# same holds for state transitions
if __name__ == "__main__":
    input_source = input("Input source(stdin, file, test): ")
    if input_source == "stdin":
        automaton_type = input("Input automaton type(Mealy, Moore): ")
        num_states, num_transitions = tuple(map(int, input("Number of states an number of transitions: ").split()))
        print("Input the automaton:")
        if automaton_type == "Mealy":
            automaton, outputs = input_mealy(None, num_states, num_transitions)
        elif automaton_type == "Moore":
            automaton, outputs = input_moore(None, num_states, num_transitions)
        else:
            raise Exception("Invalid automaton type.")
    elif input_source == "file":
        with open("input.txt","r") as file:
            file_lines = file.readlines()

        automaton_type = file_lines[0][:-1]
        num_states, num_transitions = tuple(map(int, file_lines[1].split()))
        if automaton_type == "Mealy":
            automaton, outputs = input_mealy(file_lines[2:], num_states, num_transitions)
        elif automaton_type == "Moore":
            automaton, outputs = input_moore(file_lines[2:], num_states, num_transitions)
        else:
            raise Exception("Invalid automaton type.")
    elif input_source == "test":
        test_num = input("Test number: ")
        # define file path for tests folder
        script_dir = os.path.dirname(__file__)
        rel_path = "tests/test" + test_num + ".txt"
        abs_file_path = os.path.join(script_dir, rel_path)

        with open(abs_file_path,"r") as file:
            file_lines = file.readlines()

        automaton_type = file_lines[0][:-1]
        num_states, num_transitions = tuple(map(int, file_lines[1].split()))
        if automaton_type == "Mealy":
            automaton, outputs = input_mealy(file_lines[2:], num_states, num_transitions)
        elif automaton_type == "Moore":
            automaton, outputs = input_moore(file_lines[2:], num_states, num_transitions)
        else:
            raise Exception("Invalid automaton type.")
    else:
        raise("Invalid input source.")

    minimized_automaton, minimized_outputs = minimize_automaton(automaton, outputs, (automaton_type == "Mealy"))
    # the output goes to a file output.txt
    with open("output.txt","w+") as file:
        file.write(str(len(minimized_automaton)) + " " + str(num_transitions) + "\n")
        file.writelines([" ".join([str(transition) for transition in state]) + "\n" for state in minimized_automaton])
        file.writelines([" ".join([str(output) for output in state]) + "\n" for state in minimized_outputs])