import random
from functools import reduce
from fractions import Fraction


# Divisibility
def divides(a, b):
    """Returns True if a divides b (i.e. a | b)"""
    return a != 0 and b % a == 0

def gcd(a, b):
    return gcd(b, a % b) if b != 0 else a


# Primality
def are_coprime(a, b):
    return True if gcd(a, b) == 1 else False

def primes_leq(n):
    ns = list(range(2, n + 1))
    primes = set()
    while ns:
        p = ns.pop(0)
        ns = [a for a in ns if not divides(p, a)]
        primes.add(p)
    return primes

def prime_divisors_of_(n):
    return {p for p in primes_leq(n) if divides(p, n)}

def _miller_rabin(n):
    """Returns True if num is a prime number."""
    s = n - 1
    t = 0
    while s % 2 == 0:
        # keep halving s while it is even (and use t
        # to count how many times we halve s)
        s = s // 2
        t += 1

    for _ in range(5): # try to falsify num's primality 5 times
        a = random.randrange(2, n - 1)
        v = pow(a, s, n)
        if v != 1: # this test does not apply if v is 1.
            i = 0
            while v != (n - 1):
                if i == t - 1:
                    return False
                else:
                    i = i + 1
                    v = (v ** 2) % n
    return True

def is_prime(n):
    """
    Return True if num is a prime number. This function does a quicker
    prime number check before calling rabinMiller().
    """
    if (n < 2):
        return False # 0, 1, and negative numbers are not prime

    # About 1/3 of the time we can quickly determine if num is not prime
    # by dividing by the first few dozen prime numbers. This is quicker
    # than rabinMiller(), but unlike rabinMiller() is not guaranteed to
    # prove that a number is prime.
    lowPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
        53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
        269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
        349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
        431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
        503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593,
        599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
        673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757,
        761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853,
        857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
        947, 953, 967, 971, 977, 983, 991, 997]

    if n in lowPrimes:
        return True

    # See if any of the low prime numbers can divide num
    for prime in lowPrimes:
        if (n % prime == 0):
            return False

    # If all else fails, call rabinMiller() to determine if num is a prime.
    return _miller_rabin(n)

def generate_large_prime(keysize=1024):
    """Return a random prime number of keysize bits in size."""
    while True:
        num = random.randrange(2**(keysize-1), 2**(keysize))
        if is_prime(num):
            return num

def is_odd_prime(n):
    return is_prime(n) and n != 2

def phi(n):
    """Euler totient/phi function.

    Returns the number of positive integers up to a given n that
    are relatively prime to n. The function uses Euler's product
    formula. See
    https://en.wikipedia.org/wiki/Euler%27s_totient_function
    for more info.
    """
    if n < 1:
        return None
    product = 1
    for p in prime_divisors_of_(n):
        product *= 1 - Fraction(1, p)
    return int(n * product)


# Group of Units
def Z_mod_(n):
    return {a for a in range(n)}

def U(n):
    """Returns the group of units of Z mod n"""
    return {a for a in Z_mod_(n) if are_coprime(a, n)}

def is_primitive_root(a, n):
    """Returns True if a is a primitive root modulo n"""
    phi_n = phi(n)
    for p in prime_divisors_of_(phi_n):
        if a ** (phi_n / p) % n == 1:
            return False
    return True

def primitive_roots(n):
    return {a for a in U(n) if is_primitive_root(a, n)}


# Quadratic residues
def Q(n):
    """Returns the quadratic residues modulo n"""
    return {pow(a, 2, n) for a in U(n)}

def is_quadratic_residue(a, n):
    """Returns True if a is a quadratic residue modulo n."""
    if is_odd_prime(n) and are_coprime(a, n):
        # We can use Euler's Criterion
        return True if pow(a, int((n - 1) / 2), n) == 1 else False
    return a in Q(n)

def Legendre_symbol(a, p):
    assert is_odd_prime(p)
    if divides(p, a):
        return 0
    return 1 if is_quadratic_residue(a, p) else -1
