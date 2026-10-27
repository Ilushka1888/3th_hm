"""Учебные алгоритмы теории чисел."""
from math import prod, isqrt

def crt(residues: list[int], moduli: list[int]) -> int:
    """Непустые равные списки, попарно взаимно простые модули>=2; решение в [0,prod)."""
    if len(residues) != len(moduli) or not residues:
        raise ValueError('residues and moduli must have equal non-zero length')
    if any((n < 2 for n in moduli)):
        raise ValueError('all moduli must be >= 2')
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            if gcd(moduli[i], moduli[j]) != 1:
                raise ValueError('moduli must be pairwise coprime')
    N = prod(moduli)
    x = 0
    for a_i, n_i in zip(residues, moduli):
        N_i = N // n_i
        M_i = mod_inverse(N_i, n_i)
        e_i = N_i * M_i
        x += a_i * e_i
    return x % N

def divides(a: int, b: int) -> bool:
    """Контракт: возвращает True тогда и только тогда, когда a | b (Определение 1.1).

    Крайний случай: 0 | b  <=>  b == 0.
    """
    if a == 0:
        return b == 0
    return b % a == 0

def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Не нулевая пара; возвращает (g,x,y), g=НОД(a,b)>0, ax+by=g."""
    if a == 0 and b == 0:
        raise ValueError('gcd(0, 0) is undefined by the course convention')
    A, B = (abs(a), abs(b))
    old_r, r = (A, B)
    old_s, s = (1, 0)
    old_t, t = (0, 1)
    while r != 0:
        assert old_r == old_s * A + old_t * B
        assert r == s * A + t * B
        q = old_r // r
        old_r, r = (r, old_r - q * r)
        old_s, s = (s, old_s - q * s)
        old_t, t = (t, old_t - q * t)
    assert old_r == old_s * A + old_t * B
    x = old_s if a >= 0 else -old_s
    y = old_t if b >= 0 else -old_t
    return (old_r, x, y)

def garner_crt(residues: list[int], moduli: list[int]) -> int:
    """Reconstruct the least non-negative CRT solution sequentially."""
    if len(residues) != len(moduli) or not residues:
        raise ValueError('residues and moduli must have equal non-zero length')
    if any((n < 2 for n in moduli)):
        raise ValueError('all moduli must be >= 2')
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            if gcd(moduli[i], moduli[j]) != 1:
                raise ValueError('moduli must be pairwise coprime')
    x = residues[0] % moduli[0]
    M = moduli[0]
    for a_i, m_i in zip(residues[1:], moduli[1:]):
        t = (a_i - x) * mod_inverse(M, m_i) % m_i
        x += t * M
        M *= m_i
    return x

def gcd(a: int, b: int) -> int:
    """Целые a,b, не оба нулевые; возвращает положительный НОД."""
    if a == 0 and b == 0:
        raise ValueError('gcd(0, 0) is undefined by the course convention')
    a, b = (abs(a), abs(b))
    while b:
        a, b = (b, a % b)
    return a

def mod_inverse(a: int, modulus: int) -> int:
    """modulus>=2 и НОД(a,modulus)=1; канонический обратный, иначе ValueError."""
    if modulus < 2:
        raise ValueError('modulus must be >= 2')
    g, x, _ = extended_gcd(a, modulus)
    if g != 1:
        raise ValueError('inverse does not exist')
    return x % modulus

def mod_pow(base: int, exponent: int, modulus: int) -> int:
    """exponent>=0, modulus>0; канонический остаток base**exponent."""
    if modulus <= 0:
        raise ValueError('modulus must be positive')
    if exponent < 0:
        raise ValueError('this implementation expects exponent >= 0')
    result = 1 % modulus
    base %= modulus
    e = exponent
    while e:
        if e & 1:
            result = result * base % modulus
        base = base * base % modulus
        e >>= 1
    return result

def solve_linear(a, b, n):
    """n>=2; все решения ax=b (mod n) в [0,n), по возрастанию."""
    if n < 2:
        raise ValueError('n должно быть >= 2')
    g, u, _ = extended_gcd(a, n)
    if b % g:
        return []
    step = n // g
    x0 = u * (b // g) % step
    return [x0 + k * step for k in range(g)]

def solve_linear_diophantine(a, b, c):
    """При (a,b) != (0,0): None либо (x0,y0,b/g,-a/g), ax0+by0=c."""
    g, x, y = extended_gcd(a, b)
    if c % g:
        return None
    return (x * (c // g), y * (c // g), b // g, -a // g)


def random_exact_bits(bits: int) -> int:
    """Положительное целое ровно bits бит; bits>=1, иначе ValueError."""
    import secrets
    if bits<1: raise ValueError('bits должно быть положительным')
    return secrets.randbits(bits) | (1 << (bits-1))


# Совместимый интерфейс Notebook 3; алгоритмы хранятся в primes.py.
from .primes import is_prime_fast, primes, next_prime


def decompose(N):
    """Целое N>=1; каноническое разложение {простое: положительный показатель}."""
    if N < 1:
        raise ValueError('Разложение определено для N>=1')
    n = N
    result = {}
    while n > 1:
        p = smallest_factor(n)
        result[p] = result.get(p, 0) + 1
        n //= p
    return result

def is_primitive_root(g, n, t=None):
    """n>=2; t, если передан, должен быть точным порядком группы phi(n)."""
    if n < 2:
        raise ValueError('n>=2')
    if gcd(g, n) != 1:
        return False
    if t is None:
        t = euler_phi(n)
    if t < 1:
        raise ValueError('Порядок группы положителен')
    return mod_pow(g, t, n) == 1 and all((mod_pow(g, t // q, n) != 1 for q in decompose(t)))

def multiplicative_order(a, n, t=None):
    """n>=2, gcd(a,n)=1; t>=1 — известное кратное порядка, иначе используется phi(n)."""
    if n < 2 or gcd(a, n) != 1:
        raise ValueError('Нужны n>=2 и обратимый a')
    if t is None:
        t = euler_phi(n)
    if t < 1 or mod_pow(a, t, n) != 1:
        raise ValueError('t должно быть положительным кратным порядка')
    order = t
    for q in decompose(t):
        while order % q == 0 and mod_pow(a, order // q, n) == 1:
            order //= q
    return order

def primitive_root(n):
    """Наименьший положительный первообразный корень для n>=2, либо None."""
    if n < 2:
        raise ValueError('n>=2')
    if n == 2:
        return 1
    if n == 4:
        return 3
    odd = n // 2 if n % 2 == 0 else n
    if odd % 2 == 0:
        return None
    D = decompose(odd)
    if len(D) != 1:
        return None
    t = euler_phi(n)
    qs = tuple(decompose(t))
    for g in range(2, n):
        if gcd(g, n) == 1 and all((mod_pow(g, t // q, n) != 1 for q in qs)):
            return g
    return None

def recompose(D):
    """Произведение p**e по каноническому разложению; пустой словарь даёт 1."""
    result = 1
    for p, e in D.items():
        if p < 2 or e < 1:
            raise ValueError('Нужны p>=2 и e>=1')
        result *= p ** e
    return result

def smallest_factor(n):
    """Целое n; минимальный простой делитель при n>=2, иначе None."""
    if n < 2:
        return None
    if n % 2 == 0:
        return 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return d
    return n

def euler_phi(n):
    """Целое n>=1; phi(1)=1, далее число обратимых классов по модулю n."""
    result=n
    for p in decompose(n): result=result//p*(p-1)
    return result

