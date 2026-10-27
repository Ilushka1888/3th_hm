import math
import pytest
from hypothesis import given, strategies as st, settings
from crypto_math.number_theory import *

pair = st.tuples(st.integers(),st.integers()).filter(lambda p:p!=(0,0))
@given(st.integers(),st.integers())
def test_divides(a,b):
    assert divides(a,b) == (b==0 if a==0 else b%a==0)
@given(pair)
def test_gcd(pair):
    a,b=pair
    assert gcd(a,b)==math.gcd(a,b)
@given(pair)
def test_bezout(pair):
    a,b=pair; g,x,y=extended_gcd(a,b)
    assert g==math.gcd(a,b) and a*x+b*y==g
def test_zero_pair():
    for f,args in [(gcd,(0,0)),(extended_gcd,(0,0)),(solve_linear_diophantine,(0,0,1))]:
        with pytest.raises(ValueError): f(*args)
@given(st.integers(),st.integers(min_value=2))
def test_inverse(a,n):
    if math.gcd(a,n)==1:
        u=mod_inverse(a,n)
        assert 0<=u<n and a*u%n==1 and u==pow(a,-1,n)
    else:
        with pytest.raises(ValueError): mod_inverse(a,n)
@given(st.integers(),st.integers(min_value=0,max_value=2**64),st.integers(min_value=1))
def test_pow(a,e,n):
    assert mod_pow(a,e,n)==pow(a,e,n)
@given(st.integers(),st.sampled_from([[3,5,7],[8,9,5],[89,97]]))
def test_crt(x,ms):
    rs=[x%m for m in ms]; N=math.prod(ms)
    assert crt(rs,ms)==garner_crt(rs,ms)==x%N
@given(pair,st.integers(),st.integers())
def test_diophantine(pair,c,t):
    a,b=pair; g=math.gcd(a,b); sol=solve_linear_diophantine(a,b,c)
    if c%g: assert sol is None
    else:
        x,y,dx,dy=sol
        assert a*x+b*y==c and (dx,dy)==(b//g,-a//g)
        assert a*(x+dx*t)+b*(y+dy*t)==c
        assert math.gcd(dx,dy)==1
@given(st.integers(),st.integers(),st.integers(min_value=2,max_value=150))
def test_linear(a,b,n):
    xs=solve_linear(a,b,n); g=math.gcd(a,n)
    assert xs==sorted(set(xs))
    assert len(xs)==(g if b%g==0 else 0)
    assert all(0<=x<n and (a*x-b)%n==0 for x in xs)
    assert xs==[x for x in range(n) if (a*x-b)%n==0]
def test_invalid_contracts():
    for f,args in [(mod_inverse,(1,1)),(solve_linear,(1,1,1)),(mod_pow,(2,-1,5)),(mod_pow,(2,2,0)),(crt,([],[])),(garner_crt,([1,2],[6,9]))]:
        with pytest.raises(ValueError): f(*args)


# Свойства Notebook 3.
def _divisors_small(n):
    return [d for d in range(1,n+1) if n%d==0]

def test_notebook3_prime_interfaces():
    flags=[False]*10001
    for p in primes(10000): flags[p]=True
    assert all(is_prime_fast(n)==flags[n] for n in range(10001))

@given(st.integers(min_value=1,max_value=1000000))
def test_factorization_roundtrip(n):
    D=decompose(n)
    assert recompose(D)==n
    assert all(is_prime_fast(p) and e>=1 for p,e in D.items())
    if n>1: assert smallest_factor(n)==min(D)
    else: assert D=={}

@given(st.integers(min_value=2,max_value=500),st.integers(min_value=-10000,max_value=10000))
def test_euler_and_order(n,a):
    if math.gcd(a,n)!=1:
        with pytest.raises(ValueError): multiplicative_order(a,n)
        return
    phi=euler_phi(n)
    assert mod_pow(a,phi,n)==1
    t=multiplicative_order(a,n)
    assert phi%t==0 and mod_pow(a,t,n)==1
    assert all(mod_pow(a,d,n)!=1 for d in _divisors_small(t) if d<t)
    assert multiplicative_order(a,n,phi*2)==t

@given(st.integers(min_value=1,max_value=1000))
def test_totient_divisor_identity(n):
    assert sum(euler_phi(d) for d in _divisors_small(n))==n

@given(st.integers(min_value=2,max_value=100))
def test_primitive_root_minimal(n):
    g=primitive_root(n)
    candidates=[a for a in range(1,n) if math.gcd(a,n)==1 and multiplicative_order(a,n)==euler_phi(n)]
    assert g==(min(candidates) if candidates else None)
    assert all(is_primitive_root(a,n)==(a in candidates) for a in range(n))

def test_notebook3_invalid_contracts():
    for f,args in [(decompose,(0,)),(euler_phi,(-1,)),(multiplicative_order,(2,4)),(multiplicative_order,(2,7,2)),(primitive_root,(1,)),(is_primitive_root,(2,1))]:
        with pytest.raises(ValueError): f(*args)
    assert smallest_factor(0) is None and smallest_factor(1) is None
    assert euler_phi(1)==1 and multiplicative_order(1,2)==1
