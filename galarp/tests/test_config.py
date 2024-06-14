from .. import config
from configobj.validate import Validator


def test_default_config():
    """
    Test case for the default_config function.
    """
    c = config.default_config()
    # Do type assertions
    assert isinstance(c, config.configobj.ConfigObj)


def test_default_cspec():
    """
    Test case for the default_cspec function.
    """
    cspec = config.default_cspec()
    # Do type assertions
    assert isinstance(cspec, config.configobj.ConfigObj)


def test_config_validation():
    c = config.default_config()
    cspec = config.default_cspec()

    vtor = Validator()
    c.configspec = cspec

    validated = c.validate(vtor)
    assert validated


def test_dump_config(tmp_path):
    
    config.dump_default_config(f"{tmp_path}test_config.ini")


def test_load_config(tmp_path):
    config.dump_default_config(f"{tmp_path}test_config.ini")

    c = config.load_config(f"{tmp_path}test_config.ini")
    assert isinstance(c, config.configobj.ConfigObj)