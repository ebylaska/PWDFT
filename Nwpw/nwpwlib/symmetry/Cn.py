import math

def rotation(theta):
    ct = math.cos(theta)
    st = math.sin(theta)
    # 3x3 rotation about z, flattened
    print(f"{ct: .12f} {st: .12f} 0.0  "
          f"{-st: .12f} {ct: .12f} 0.0  "
          f"0.0 0.0 1.0")

if __name__ == "__main__":
    for order in range(1, 11):
        print(f"C{order}")
        print(order)
        for i in range(order):
            theta = 2.0 * math.pi * i / order
            rotation(theta)
        print()

