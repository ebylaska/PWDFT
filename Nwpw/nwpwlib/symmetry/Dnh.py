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

def sigma_h():
    return [
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
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

def sigma_v(theta):
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
        print(f"D{n}h")
        print(4 * n)

        # D_n rotations
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz(theta))

        # C2' axes
        for k in range(n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz_Rx_pi(theta))

        # sigma_h
        print_op(sigma_h())

        # improper rotations S_n
        for k in range(1, n):
            theta = 2.0 * math.pi * k / n
            print_op(Rz_sigma_h(theta))

        # vertical mirrors
        for k in range(n):
            theta = math.pi * k / n
            print_op(sigma_v(theta))

        print()

