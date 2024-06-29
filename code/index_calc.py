from sage.all import *
from helpers import find_prime_base, gen_random_powers, modpow, find_smooth_powers, prime_factor, solve_strt_system

def solveDLP(g, h, p, B):
    """
    Given prime p, g primitive root mod p, h in F_p, and hyperparamter B, use the index calculus method to solve the DLP:
    g^x cong. h (mod p)
    Outputs x in (0, p)
    """

    # STEP 1: Manually find the factor base of primes bounded by B
    factor_base = find_prime_base(B)

    # STEP 2: Take random powers of g modulo p and pick pi(B) ones that  are B-smooth. 
    pi_B = len(factor_base)
    
    # STEP 3: Convert these into linear equations with variables equal to the log_g(l) for l in the factor base 
    coefficients, powers_i = gen_random_powers(g, p, pi_B, factor_base)

    # factor p - 1 
    prime_factorization = prime_factor(p - 1)

    successful_random_pows = False
    while not successful_random_pows:
        successful_random_pows = True
        sun_tzu_consts = []
        for q_i, e_i in prime_factorization:
            # NOTE: Simplifying assumption that p-1 is square-free
            if e_i > 1: raise Exception("Out of scope: p-1 is not square-free")
            # STEP 4: Solve the congruencies modulo the prime factors of p - 1
            coeffs_matrix = Matrix(GF(q_i), coefficients)
            consts_matrix = vector(GF(q_i), powers_i)
            aug = coeffs_matrix.augment(consts_matrix)
            R = aug.rref()
            # extract the roots from the matrix
            system_solns = []
            for i in range(pi_B):
                if R[i, i] != 1: # extract from pivot rows
                    continue
                system_solns.append(R[i, pi_B])

            if len(system_solns) != pi_B: 
                # Oh no, this choice of powers actually doesn't give us a solution. Try again with a new set of random powers!
                coefficients, powers_i = gen_random_powers(g, p, pi_B, factor_base)
                successful_random_pows = False

            sun_tzu_consts.append((q_i, system_solns))

    # STEP 6: Combine the solutions using STRT
    sun_tzu_solved = solve_strt_system(sun_tzu_consts, pi_B)

    # STEP 7: Choose random values of k until we find a value of hg^{-k} that is B-smooth 
    while True:
        k = randint(1, p-1)
        hg_minus_k = (h * modpow(g, -k, p)) % p
        smooth_powers = find_smooth_powers(hg_minus_k, factor_base)
        if smooth_powers != None:
            # STEP 8: Calculate the DLP by taking the linear combo of log_g(l)
            x = k
            for sm_pwr, log_g in zip(smooth_powers, sun_tzu_solved):
                x += sm_pwr * log_g
            return x % (p - 1)

if __name__ == "__main__":
    g = int(input("enter a value for g: "))
    h = int(input("enter a value for h: "))
    p = int(input("enter a value for p: "))
    B = int(input("enter a value for B: "))
    print("The solution to the DLP is:", solveDLP(g, h, p, B))

