"""Проверенные алгоритмы Notebook 2."""
import math

def is_prime_fast(n):
    """Целое n; True ровно для простых n. Перебор нечётных делителей до isqrt."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, math.isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True

def isprime_bytes(n):
    """Компактное решето, один байт на индекс; граница n включена."""
    if n < 0:
        raise ValueError('n должно быть неотрицательным')
    flags = bytearray(b'\x01') * (n + 1)
    flags[0] = 0
    if n >= 1:
        flags[1] = 0
    if n >= 4:
        flags[4::2] = b'\x00' * ((n - 4) // 2 + 1)
    for p in range(3, math.isqrt(n) + 1, 2):
        if flags[p]:
            flags[p * p::2 * p] = b'\x00' * ((n - p * p) // (2 * p) + 1)
    return flags

def next_prime(n):
    """Для целого n возвращает наименьшее простое строго больше n."""
    if n < 2:
        return 2
    candidate = n + 1
    if candidate % 2 == 0:
        candidate += 1
    while not is_prime_fast(candidate):
        candidate += 2
    return candidate

def previous_prime(n):
    """Целое n>2; наибольшее простое строго меньше n. При n<=2 ValueError."""
    if n <= 2:
        raise ValueError('Простого числа строго меньше n не существует')
    if n <= 3:
        return 2
    candidate = n - 1
    if candidate % 2 == 0:
        candidate -= 1
    while candidate >= 3:
        if is_prime_fast(candidate):
            return candidate
        candidate -= 2
    return 2

def prime_pi_exact(n):
    """Точный комбинаторный подсчёт π(n) с O(sqrt(n)) хранимых значений."""
    if n < 2:
        return 0
    r = math.isqrt(n)
    values = sorted(set(range(1, r + 1)) | {n // i for i in range(1, r + 1)}, reverse=True)
    counts = {v: v - 1 for v in values}
    for p in range(2, r + 1):
        if counts[p] == counts[p - 1]:
            continue
        smaller = counts[p - 1]
        square = p * p
        for v in values:
            if v < square:
                break
            counts[v] -= counts[v // p] - smaller
    return counts[n]

def prime_range(a, b):
    """Целые a,b; все простые в полуинтервале [a,b)."""
    if a >= b or b <= 2:
        return []
    flags = isprime_bytes(b - 1)
    return [p for p in range(max(2, a), b) if flags[p]]

def primes(n):
    """Все простые <=n; целое n>=0, иначе ValueError. Учебное решето с циклами."""
    if n < 0:
        raise ValueError('n должно быть неотрицательным')
    flags = bytearray(b'\x01') * (n + 1)
    flags[0] = 0
    if n >= 1:
        flags[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if flags[p]:
            for multiple in range(p * p, n + 1, p):
                flags[multiple] = 0
    return [p for p in range(n + 1) if flags[p]]
