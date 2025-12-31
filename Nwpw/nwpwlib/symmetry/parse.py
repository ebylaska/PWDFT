def parse_spacegroup_file(filename):
    groups = []

    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    nlines = len(lines)

    while i < nlines:
        # Group name
        group_name = lines[i]
        i += 1

        # Group number (rank)
        group_number = int(lines[i])
        i += 1

        # Number of symmetry operations
        num_ops = int(lines[i])
        i += 1

        # Skip symmetry operation matrices
        # Each operation has 3 rows
        i += num_ops * 3

        groups.append({
            "group_name": group_name,
            "group_number": group_number,
            "num_operations": num_ops
        })

    return groups


if __name__ == "__main__":
    filename = "spacegroups.dat"  # change as needed
    groups = parse_spacegroup_file(filename)

    for g in groups:
        print(
            f"{g['group_name']:10s}  "
            f"{g['group_number']:3d}  "
            f"{g['num_operations']:3d}"
        )
        #print(
        #    f"{g['group_name']:10s}  "
        #    f"number={g['group_number']:3d}  "
        #    f"ops={g['num_operations']:3d}"
        #)

