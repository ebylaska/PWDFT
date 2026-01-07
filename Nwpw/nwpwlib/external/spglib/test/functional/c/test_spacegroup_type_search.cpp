#include <gtest/gtest.h>

extern "C" {
#include "spglib.h"
}

TEST(SpacegroupTypeSearch, test_spg_get_international) {
    double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
    double position[][3] = {
        {0, 0, 0},     {0.5, 0.5, 0.5}, {0.3, 0.3, 0},
        {0.7, 0.7, 0}, {0.2, 0.8, 0.5}, {0.8, 0.2, 0.5},
    };
    int types[] = {1, 1, 2, 2, 2, 2};
    int num_spg;
    int num_atom = 6;
    char symbol[21];
    char symbol_expect[] = "P4_2/mnm";

    num_spg =
        spg_get_international(symbol, lattice, position, types, num_atom, 1e-5);
    ASSERT_EQ(num_spg, 136);
    ASSERT_STREQ(symbol, symbol_expect);
}

TEST(SpacegroupTypeSearch, test_spg_get_schoenflies) {
    double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
    double position[][3] = {
        {0, 0, 0},     {0.5, 0.5, 0.5}, {0.3, 0.3, 0},
        {0.7, 0.7, 0}, {0.2, 0.8, 0.5}, {0.8, 0.2, 0.5},
    };
    int types[] = {1, 1, 2, 2, 2, 2};
    int num_atom = 6;
    int num_spg;
    char symbol[7];
    char symbol_expect[] = "D4h^14";

    num_spg =
        spg_get_schoenflies(symbol, lattice, position, types, num_atom, 1e-5);
    ASSERT_EQ(num_spg, 136);
    ASSERT_STREQ(symbol, symbol_expect);
}

TEST(SpacegroupTypeSearch, test_spg_get_spacegroup_type_from_symmetry) {
    SpglibDataset *dataset;
    SpglibSpacegroupType spg_type;

    // Structure derived from 1509692 in the crystallography open database.
    double lattice[3][3] = {{0.0000, 5.245, 0.0000},
                            {6.4345, -6.4345, 0.0000},
                            {1.8319, 0.0000, -3.6638}};
    double position[][3] = {{0.6416, 0.7350, 0.0708}, {0.8284, 0.2350, 0.4142},
                            {0.3584, 0.2650, 0.4292}, {0.1716, 0.7650, 0.5858},
                            {0.6994, 0.5334, 0.2150}, {0.3674, 0.0334, 0.2990},
                            {0.3006, 0.4666, 0.5156}, {0.0000, 0.0000, 0.5902},
                            {0.0000, 0.5000, 0.8402}, {0.6326, 0.9666, 0.9316}};
    int types[] = {1, 1, 1, 1, 6, 6, 6, 6, 6, 6};
    int num_atom = 10;

    dataset = spg_get_dataset(lattice, position, types, num_atom, 1e-5);
    ASSERT_NE(dataset, nullptr);
    EXPECT_EQ(dataset->hall_number, 212);

    spg_type = spg_get_spacegroup_type_from_symmetry(
        dataset->rotations, dataset->translations, dataset->n_operations,
        lattice, 1e-5);

    EXPECT_EQ(spg_type.hall_number, 212);

    if (dataset) {
        spg_free_dataset(dataset);
    }
}
