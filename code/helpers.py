from sage.all import *
def modpow(a,n,m):
    """
    Given a, n, and m, return a^n mod m using L-to-R binary fast powering
    """
    # if n is negative, take the inverse of a
    if n < 0:
        # since au - mv = 1, we see that au is 1 mod m, hence u is the inverse of a
        a, _ = ext_gcd(a, m)
        n = -1 * n
    # represent n in binary 
    n_bin = bin(n)[2:]
    prod = 1
    for digit in n_bin:
        # square the product
        prod = (prod ** 2) % m
        if digit == '1':
            prod = (prod * a) % m
    return prod 

def ext_gcd(a,b):
    """
    Use extended gcd to compute a inverse mod m (the first return value)
    """
    x, y, z, w = 1, 0, 0, 1
    while b!= 0:
        x, y, z, w = z, w, x - (a//b)*z, y - (a//b)*w
        a, b = b, a%b
    return x, y

def find_prime_base(B):
    """
    Given a number, find all of the prime numbers p such that 2 <= p <= B
    Uses the Sieve of Eratosthenes
    """
    # Keeps track of which numbers are prime in the range
    primes = [True for i in range(B + 1)]
    # The first prime is 2 :D
    p = 2
    while (p ** 2 <= B):
        if primes[p]: 
            for i in range (p ** 2, B + 1, p):
                primes[i] = False
        p += 1

    output_primes = []
    for p in range(2, B + 1):
        if(primes[p]):
            output_primes.append(p)

    return output_primes


def find_smooth_powers(n, factor_base):
    """
    Given a factor base for a number B, determine whether a number n is B smooth.
    If so return a list of exponents corresponding to the factor base. If not, return None
    """
    output = []
    for b in factor_base:
        # for each of the factor base numbers, factor out as much as you can 
        power = 0
        while n % b == 0:
            power += 1
            n = n // b
        output.append(power)
    # if we perfectly factored, success! return the powers. 
    if n == 1: 
        return output
    else:
        return None


def gen_random_powers(g, p, pi_B, factor_base):
    coefficients = []
    powers_i = []
    
    while len(coefficients) <  pi_B: 
        rand_pwr_i = randint(1, p-1)
        u_l = find_smooth_powers(modpow(g, rand_pwr_i, p), factor_base)
        if u_l != None:
            coefficients.append(u_l)
            powers_i.append(rand_pwr_i)

    return (coefficients, powers_i)

def prime_factor(n):
    """
    Factor a number n. Returns a list of tuples, the first elt being the prime factor and the second being its power
    """
    # find all of the possible prime factors of n
    primes = find_prime_base(n)
    # now find how they factorize n
    pows = find_smooth_powers(n, primes)
    output = [] 
    for prime, power in zip(primes, pows):
        if power != 0:
            output.append((prime, power))
    return output

def solve_strt_system(sun_tzu_consts, pi_B):
    """
     Given a Sun Tzu system of the form 
     [(modulus, [list of constants])] 
     Solve the system and return (modulus, [list of solutions])
    """
    solns = []
    # for each base
    for i in range(pi_B):
        # for each prime 
        strt_in = []
        m = 1
        for j in range(len(sun_tzu_consts)):
            m_j = sun_tzu_consts[j][0]
            strt_in.append((m_j, int(sun_tzu_consts[j][1][i])))
            m *= m_j
        solns.append(solve_single_strt(m, strt_in))
    return solns


def solve_single_strt(m, strt_in):
    """
    Given input of the form [(base_1, const_1), (base_2, const_2), ...] return 
    the solution to the STRT system, where m is the product of all of the bases 
    """
    soln = 0
    for m_i, b_i in strt_in:
        n_i = m // m_i
        # know that r_i * m_i + n_i * s_i = 1 
        _, s_i = ext_gcd(m_i, n_i)
        e_i = n_i * s_i
        soln += b_i * e_i
    return soln % m





    












