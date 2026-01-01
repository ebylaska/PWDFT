import math

def rotation(theta):
    ct = math.cos(theta)
    st = math.sin(theta)
    return [
        ct,  st, 0.0,
       -st,  ct, 0.0,
        0.0, 0.0, 1.0
    ]

def mirror(theta):
    # σ_v = R(theta) * diag(1,-1,1)
    ct = math.cos(theta)
    st = math.sin(theta)
    return [
         ct, -st, 0.0,
        -st, -ct, 0.0,
         0.0, 0.0, 1.0
    ]

def print_op(op):
    print(" ".join(f"{x: .12f}" for x in op))

if __name__ == "__main__":

    for order in range(1, 11):
        print(f"C{order}v")
        print(2 * order)

        # Rotations
        for i in range(order):
            theta = 2.0 * math.pi * i / order
            print_op(rotation(theta))

        # Vertical mirrors
        for i in range(order):
            theta = 2.0 * math.pi * i / order
            print_op(mirror(theta))

        print()

