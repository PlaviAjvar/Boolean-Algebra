import os

# reduce two implicants which we assume are reducible
def reduce(first_implicant, second_implicant):
    form, mask, _ = first_implicant
    # add new reduced bit to mask
    new_mask = mask | (second_implicant[0] - first_implicant[0])
    is_dont_care = (first_implicant[2] and second_implicant[2])
    return (form, new_mask, is_dont_care)

# test if two implicants are reducible
def can_reduce(first_implicant, second_implicant):
    # if the input labels are different we can't reduce
    if first_implicant[1] != second_implicant[1]:
        return False
    # it's necessary and sufficient that second_implicant-first_implicant is a power of 2
    if second_implicant[0] <= first_implicant[0]:
        return False
    difference = second_implicant[0] - first_implicant[0]
    # n is a power of two if and only if n & (n-1) is zero
    # it's easy to see why, looking at the binary form of the number
    return not(difference & (difference-1))

# calculate bitcount of number
def bitcount(n):
    return bin(n).count("1")

def get_prime_implicants(constituents, num_var, dont_care):
    # sort one constituents by their bitcount
    implicant_buckets = [set() for _ in range(num_var+1)]
    # the third parameter is True if it's a don't care combination
    for constituent in constituents:
        implicant_buckets[bitcount(constituent)].add((constituent, 0, False))
    # add don't care combinations
    for dc_comb in dont_care:
        implicant_buckets[bitcount(dc_comb)].add((dc_comb, 0, True))

    # find all prime implicants
    # if a prime implicant is obtained by reducing only don't care combinations, don't list it
    # it's not hard to see it will never be in the optimal solution to the problem
    prime_implicants = []
    while True:
        next_buckets = [set() for _ in range(num_var)]
        any_reductions = False
        # tells us which implicants have been used in the iteration
        label = [[False] * len(bucket) for bucket in implicant_buckets]

        # try to reduce neighboring implicants (bitcount off by one)
        for bucket_idx, bucket in enumerate(implicant_buckets[1:], start=1):
            last_bucket = implicant_buckets[bucket_idx-1]

            for implicant_idx, implicant in enumerate(bucket):
                for old_implicant_idx, old_implicant in enumerate(last_bucket):
                    if can_reduce(old_implicant, implicant):
                        # label them as used
                        label[bucket_idx-1][old_implicant_idx] = True
                        label[bucket_idx][implicant_idx] = True
                        # add reduced implicant to next step of algorithm
                        # the bitcount will be the lesser of the two bitcounts
                        next_buckets[bitcount(old_implicant[0])].add(reduce(old_implicant, implicant))
                        any_reductions = True

        # all unlabeled implicants are therefore prime implicants
        for bucket_idx, bucket in enumerate(implicant_buckets):
            for implicant_idx, implicant in enumerate(bucket):
                if not label[bucket_idx][implicant_idx] and implicant[2] != True:
                    prime_implicants.append(implicant[:-1])

        # we terminate if none of the implicants were able to reduce
        if not any_reductions:
            break
        implicant_buckets = next_buckets

    return prime_implicants

# returns True if constituent is implied by prime implicant being true
def implies(implicant, constituent):
    # the only differing bits should be the reduced bits
    # bitwise xor gives all bits which are different
    pattern, mask = implicant
    difference = pattern ^ constituent
    # all differing bits should be submask of mask
    # Ergo, difference & mask = difference (bitwise AND)
    return (difference & mask) == difference

# return bitmask which represents which constituents are implied by the prime implicant
def cover(implicant, constituents):
    mask = 0
    for constituent_idx, constituent in enumerate(constituents):
        if implies(implicant, constituent):
            mask |= (1 << constituent_idx)
    return mask


# core method for performing Quine-McCluskey algorithm
def quine_mccluskey(constituents, num_var, dont_care):
        # get number of letters used in implicant
        # conjuctions, disjunctions and negations don't count
        def word_length(implicant):
            return num_var - bitcount(implicant[1])

        # get the prime implicants
        prime_implicants = get_prime_implicants(constituents, num_var, dont_care)

        # find minimal set of prime implicants
        # we will use a bitmask dynamic programming approach
        # the number of prime implicants could be upto 3^n ln(n)
        # the number of one constituents is max 2^n
        # this means that the it's better in the worst case that the exponential is over the constituents
        inf = 10**9
        num_implicants = len(prime_implicants)
        num_constituents = len(constituents)
        # minimal_form[i][mask] is the length of the minimal form which uses the first i-1 implicants
        # and implies the constituents covered by mask
        # for certain implementation reasons, minimal_form[i][mask] is larger than it should be by exactly one
        # for all non-trivial cases
        minimal_form = [[inf] * 2**num_constituents for _ in range(num_implicants+1)]
        backtrack = [[None] * 2**num_constituents for _ in range(num_implicants+1)]
        # base case
        minimal_form[0][0] = 0
        backtrack[0][0] = 0
        # bottom up DP, also forward (better for bitwise ones like this)
        for implicant_idx, implicant in enumerate(prime_implicants):
            for mask in range(2**num_constituents):
                # if the state is unreachable (infinite) we skip it
                if minimal_form[implicant_idx][mask] >= inf:
                    continue
                    
                # first case - we don't use the current implicant
                # the mask stays the same
                if minimal_form[implicant_idx+1][mask] > minimal_form[implicant_idx][mask]:
                    minimal_form[implicant_idx+1][mask] = minimal_form[implicant_idx][mask]
                    backtrack[implicant_idx+1][mask] = mask

                # second case - we use the current implicant
                # we do a union of the mask and the cover of the current implicant (bitwise OR)
                union = mask | cover(implicant, constituents)
                if minimal_form[implicant_idx + 1][union] > minimal_form[implicant_idx][mask] + word_length(implicant):
                    minimal_form[implicant_idx + 1][union] = minimal_form[implicant_idx][mask] + word_length(implicant)
                    backtrack[implicant_idx + 1][union] = mask

        # now we can reconstruct the optimal solution from the state (num_implicants, 11...1)
        solution = []
        cur_mask = 2**num_constituents - 1

        for implicant_idx in reversed(range(num_implicants)):
            new_mask = backtrack[implicant_idx+1][cur_mask]
            if new_mask != cur_mask:
                solution.append(prime_implicants[implicant_idx])
            cur_mask = new_mask

        return solution, (minimal_form[num_implicants][2**num_constituents-1])

# utility function for inverting inputs in implicant
def invert_input(implicant, num_var):
    return ((implicant[0] ^ (2**num_var-1)), implicant[1])

# inverts inputs for negated disjunctive form
# this gets us the nonnegated conjuctive form inputs
def invert_inputs(implicants, num_var):
    return [invert_input(implicant, num_var) for implicant in implicants]

# invert the one constituents of the expression
# if K was in the expression, it no longer is
# if K wasn't in the expression, it now is
# dont cares stay the same
def invert_constituents(one_constituents, num_var, dont_care):
    in_expr = [False for _ in range(2**num_var)]
    for one_constituent in one_constituents:
        in_expr[one_constituent] = True
    for dc_comb in dont_care:
        in_expr[dc_comb] = True
    # the remaining False entries were zeros in the original function
    # fill the list with these values
    inverted_list = []
    for one_constituent in range(2**num_var):
        if not in_expr[one_constituent]:
            inverted_list.append(one_constituent)
    return inverted_list

# utility function for getting string from disjunctive/conjuctive form
# uses uppercase letters of the english alphabet as input labels
def string_from_expression(constituents, num_var, is_dnf):
    expr_str = ""

    for index, (pattern, mask) in enumerate(constituents):
        # add conjuction between constituents
        if index > 0:
            if is_dnf:
                expr_str += "v"
        # if it's in conjuctive form put parenthesis before implicant
        if not is_dnf:
            expr_str += "("
        # iterate through bits of implicant, output active (nonreduced) bits
        # first factor only applies if is_dnf = False (conjuctive form)
        # we need to see if we've put the first OR in the current conjuction, yet
        first_factor = True
        for bit in reversed(range(num_var)):
            if mask & (1 << bit) == 0:
                # add disjuntion after every input label starting from the second one
                if not first_factor and not is_dnf:
                    expr_str += "v"
                expr_str += chr(ord("A") + num_var-bit-1)
                first_factor = False
                if not (pattern & (1 << bit)):
                    expr_str += "'"

        # at this point if first_factor = True, means there were no factors
        if first_factor:
            expr_str += "1"
        # if it's in conjuctive form put parenthesis after implicant
        if not is_dnf:
            expr_str += ")"

    return expr_str

# input logic expression in the form of thee set of one constituents
# a one constituent is a conjuction of all input terms(possibly negated)
# for the logic function y = ABC v AB'C'
# the one constituents are {ABC, AB'C'} for the input terms A,B,C

# if we define an ordering of the input variables(A,B,C), we can define a binary number for each constituent
# for example ABC would be 111 and AB'C' would be 100
# therefore the input is a set of binary integers from 0 to 2^#inputs

# we need to minimize this form, by combining and reducing one constituents
# for example ABC v ABC' = AB
# the form with the minimal number of letters (A or A' are one letter, conjuction and disjunction are one)

# another form we can write expressions in is the conjuctive form, for example y = (AvBvC)(AvB'vC')
# we will also find the minimal conjuctive form for the expression
# the minimum conjuctive form is equivalent to the minimum disjunctive form of the y' (inverted function)

if __name__ == "__main__":
    # input the constituents
    print("Select input method (standard, file, test):")
    input_method = input()
    if input_method == "standard":
        num_var = int(input("Number of variables:"))
        one_constituents = list(map(int, input("Input constituent ones:").split()))
        dont_care = list(map(int, input("Input don't care combinations:").split()))
    elif input_method == "file":
        with open("input.txt", "r") as file:
            file_lines = file.readlines()
        num_var = int(file_lines[0])
        one_constituents = list(map(int, file_lines[1].split()))
        if len(file_lines) > 2:
            dont_care = list(map(int, file_lines[2].split()))
        else:
            dont_care = []
    elif input_method == "test":
        test_number = int(input("Input test number: "))
        # define file path for tests folder
        script_dir = os.path.dirname(__file__)
        rel_path = "tests/test" + str(test_number) + ".txt"
        abs_file_path = os.path.join(script_dir, rel_path)

        with open(abs_file_path, "r+") as file:
            file_lines = file.readlines()
        num_var = int(file_lines[0])
        one_constituents = list(map(int, file_lines[1].split()))
        if len(file_lines) > 2:
            dont_care = list(map(int, file_lines[2].split()))
        else:
            dont_care = []

        print("The results of testing on test" + str(test_number) + ":\n")
    else:
        raise Exception("Invalid input method.")

    # find the minimal forms using the Quine-McCluskey algoritm
    minimal_disjunctive, dis_length = quine_mccluskey(one_constituents, num_var, dont_care)
    # sort the implicants for clarity
    minimal_disjunctive = sorted(minimal_disjunctive)

    minimal_conjuctive, con_length = quine_mccluskey(invert_constituents(one_constituents, num_var, dont_care), num_var, dont_care)
    minimal_conjuctive = invert_inputs(minimal_conjuctive, num_var)
    minimal_conjuctive = sorted(minimal_conjuctive)

    # print the results
    print("\nThe minimal disjunctive form is:", string_from_expression(minimal_disjunctive, num_var, is_dnf = True))
    print("The corresponding prime implicants are", minimal_disjunctive)
    print("It's length is", dis_length)

    print("\nThe minimal conjuctive form is:", string_from_expression(minimal_conjuctive, num_var, is_dnf = False))
    print("The corresponding prime implicants are", minimal_conjuctive)
    print("It's length is", con_length)