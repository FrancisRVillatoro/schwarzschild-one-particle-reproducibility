#!/usr/bin/env python3
"""Closed-form leading ultraviolet tails used in manuscript v7.1."""
import math
c7=-1/2835
d7=4/2835

def C_tail(L):
    # (1/2pi^2) int_L^inf lambda^2 c7 lambda^-7 dlambda
    return c7/(8*math.pi**2*L**4)

def Delta_tail(L):
    # (1/4pi^2) int_L^inf lambda^4 d7 lambda^-7 dlambda
    return d7/(8*math.pi**2*L**2)

print("c7 =",c7," = -1/2835")
print("d7 =",d7," =  4/2835")
print("C_dyn leading tail beyond 8   =",C_tail(8.0))
print("Delta_dyn leading tail beyond 16 =",Delta_tail(16.0))
