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
    c = math.cos(theta)
    s = math.sin(theta)
    return [
         c, -s, 0.0,
        -s, -c, 0.0,
         0.0, 0.0,-1.0
    ]

def Rz_sigma_h(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return [
         c,  s, 0.0,
        -s,  c, 0.0,
         0.0, 0.0,-1.0
    ]

def sigma_d(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return [
         c*c - s*s,  2*c*s, 0.0,
         2*c*s,      s*s - c*c, 0.0,
         0.0,        0.0, 1.0
    ]

def print_op(op):
    print(" ".join(f"{x: .12f}" for x in op))

if __name__ == "__main__":
    for n in range(2, 11):
        print(f"D{n}d")
        print(4 * n)

        # C_n rotations
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz(theta))

        # C2' rotations
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz_Rx_pi(theta))

        # S_2n improper rotations
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz_sigma_h(theta))

        # sigma_d planes
        for k in range(n):
            theta = math.pi / (2*n) + math.pi * k / n
            print_op(sigma_d(theta))

        print()

