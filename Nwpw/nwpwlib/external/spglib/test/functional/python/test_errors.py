from __future__ import annotations

import warnings

import pytest
import spglib
from spglib import niggli_reduce


@pytest.mark.parametrize("value", ["0", "1", "true", "false", "True", "False"])
def test_deprecation_warning_env(monkeypatch, value: str):
    monkeypatch.setenv("SPGLIB_OLD_ERROR_HANDLING", value)

    if value in ("1", "true", "True"):
        with pytest.warns(DeprecationWarning):
            warnings.simplefilter("always", DeprecationWarning)
            niggli_reduce([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            niggli_reduce([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


@pytest.mark.parametrize("value", [True, False])
def test_deprecation_warning_OLD_ERROR_HANDLING(monkeypatch, value: bool):
    monkeypatch.setattr(spglib.error, "OLD_ERROR_HANDLING", value)

    if value:
        with pytest.warns(DeprecationWarning):
            warnings.simplefilter("always", DeprecationWarning)
            niggli_reduce([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            niggli_reduce([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
