import pytest
from spglib.error import SpglibError
from spglib.utils import _expand_cell


def test_expand_cell():
    with pytest.raises(SpglibError, match="cell has to be a tuple of 3 or 4 elements."):
        _expand_cell(([], [], [], [], []))
