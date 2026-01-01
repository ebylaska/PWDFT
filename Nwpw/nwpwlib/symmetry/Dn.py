import math

def Rz(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return [
         c,  s, 0.0,
        -s,  c, 0.0,
         0.0, 0.0, 1.0
    ]

def Rz_Rx_pi(theta):
    # Rz(theta) * Rx(pi)
    c = math.cos(theta)
    s = math.sin(theta)
    return [
         c, -s, 0.0,
        -s, -c, 0.0,
         0.0, 0.0,-1.0
    ]

def print_op(op):
    print(" ".join(f"{x: .12f}" for x in op))

if __name__ == "__main__":
    for n in range(2, 11):   # D1 is not used in Cotton chemistry tables; start at 2
        print(f"D{n}")
        print(2 * n)

        # n rotations about z
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz(theta))

        # n C2 rotations perpendicular to z
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz_Rx_pi(theta))

        print()

