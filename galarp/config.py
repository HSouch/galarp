import configobj as configobj


def default_config():

    config = configobj.ConfigObj()

    config["OUT_DIR"] = "output/"

    # Wind parameters
    config["WindType"] = "constant"
    config["WindSpeed"] = 500.0
    config[""]
