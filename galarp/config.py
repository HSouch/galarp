import configobj as configobj
from configobj.validate import Validator


def default_config():

    config = configobj.ConfigObj()


    config["OUT_DIR"] = "output/"

    # Potential Parameters
    config["POTENTIAL"] = {}
    config["POTENTIAL"]["DMPotentialName"] = "burkert"
    config["POTENTIAL"]["DMPotentialParams"] = [15.]
    config["POTENTIAL"]

    # Wind parameters
    config["WIND"] = {}
    config["WIND"]["WindType"] = "constant"
    config["WIND"]["WindSpeed"] = 500.0
    config["WIND"]["WindInclinationDeg"] = 45.0

    config["WIND"]["LorentzianWindT0Myr"] = 500.0
    config["WIND"]["LorentzianWindWidthMyr"] = 100
    config["WIND"]["InterpolationFilename"] = "None"
    config["WIND"]["InterpolatedWindTKey"] = "col1"
    config["WIND"]["InterpolatedWindVKeys"] = ["col4", "col5", "col6"]

    config["WIND"]["DensityType"] = "constant"
    config["WIND"]["DensityGasDensityGCM3"] = 1e-26 
    return config


def default_cspec():
    cspec = configobj.ConfigObj()
    
    cspec["POTENTIAL"] = {}
    cspec["POTENTIAL"]["DMPotentialName"] = "option('burkert', 'nfw', 'mn', 'hq', 'None', default='burkert')"
    
    cspec["WIND"] = {}
    cspec["WIND"]["WindType"] = "option('constant', 'lorentzian', 'builtin', 'interpolated', 'None', default='constant')"
    cspec["WIND"]["WindSpeed"] = "float(default=500.0)"
    cspec["WIND"]["WindInclinationDeg"] = "float(default=45.0)"

    cspec["WIND"]["LorentzianWindT0Myr"] = "float(default=500.0)"
    cspec["WIND"]["LorentzianWindWidthMyr"] = "float(default=100)"
    cspec["WIND"]["InterpolationFilename"] = "string(default='None')"
    cspec["WIND"]["InterpolatedWindTKeys"] = "string(default='col1')"
    cspec["WIND"]["InterpolatedWindVKeys"] = "list(default=['col4', 'col5', 'col6'])"

    cspec["WIND"]["DensityType"] = "option('constant', 'builtin', 'interpolated', 'None', default='constant')"
    cspec["WIND"]["DensityGasDensityGCM3"] = "float(default=1e-26)"

    return cspec



def validate_config(config):
   
    vtor = Validator()
    config.configspec = default_cspec()

    validated = config.validate(vtor)


    print(config["WIND"]["WindSpeed"], type(config["WIND"]["WindSpeed"]))   


    for entry in configobj.flatten_errors(config, validated):
        # each entry is a tuple
        section_list, key, error = entry
        if key is not None:
            section_list.append(key)
        else:
            section_list.append('[missing section]')
        section_string = ', '.join(section_list)
        if error == False:
            error = 'Incorrect value or section.'
        print (section_string, ' = ', error)


def dump_default_config(path):
    config = default_config()
    config.filename = path
    config.write()


def load_config(path):
    config = configobj.ConfigObj(path)
    validate_config(config)
    return config


if __name__ == "__main__":

    load_config("config.ini")

