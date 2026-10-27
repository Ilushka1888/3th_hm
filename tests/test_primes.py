import math
import pytest
from hypothesis import given, strategies as st
from crypto_math import random_exact_bits
from crypto_math.primes import (is_prime_fast, primes, prime_range, next_prime, previous_prime, isprime_bytes, prime_pi_exact)

def reference(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
@given(st.integers(min_value=-1000,max_value=100000))
def test_prime_reference(n):
    assert is_prime_fast(n)==reference(n)
@given(st.integers(min_value=0,max_value=2000))
def test_sieve(n):
    expected=[p for p in range(n+1) if reference(p)]
    assert primes(n)==expected
    assert [p for p,f in enumerate(isprime_bytes(n)) if f]==expected
    assert prime_pi_exact(n)==len(expected)
@given(st.integers(min_value=-50,max_value=1000),st.integers(min_value=-50,max_value=1000))
def test_range(a,b):
    assert prime_range(a,b)==[p for p in range(a,b) if reference(p)]
@given(st.integers(min_value=-100,max_value=10000))
def test_next(n):
    p=next_prime(n)
    assert p>n and reference(p)
    assert not any(reference(k) for k in range(max(2,n+1),p))
@given(st.integers(min_value=3,max_value=10000))
def test_previous(n):
    p=previous_prime(n)
    assert p<n and reference(p)
    assert not any(reference(k) for k in range(p+1,n))
@given(st.integers(min_value=1,max_value=4096))
def test_bits(bits):
    assert random_exact_bits(bits).bit_length()==bits
def test_edges():
    for n in (-100,0,1,2):
        with pytest.raises(ValueError): previous_prime(n)
    for f in (primes,isprime_bytes):
        with pytest.raises(ValueError): f(-1)
    with pytest.raises(ValueError): random_exact_bits(0)
    assert prime_range(2,3)==[2] and prime_range(3,3)==[]
