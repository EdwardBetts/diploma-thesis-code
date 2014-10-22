cdef double cfactorial(double n):
    if n == 0 or n == 1:
        return 1
    else:
        return cfactorial(n - 1) * n

def factorial(double n):
    return cfactorial(n)
