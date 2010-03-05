"""Base for a bunch of a resequencing model stuff.
"""

class ParamDefinition:
    def __init__(self, name, default_string, ftype, description):
        """
        @param name: the name of the parameter
        @param default_string: a string representing the default value
        @param ftype: a function that deserializes from a string
        @param description: a short description of the parameter
        """
        self.name = name
        self.default_string = default_string
        self.ftype = ftype
        self.description = description

def read_params(param_definitions, name_value_pairs, conf):
    """Read parameters in a given format.
    @param model_definition_rows: parameter definitions
    @param name_value_pairs: pairs of strings of user input
    @param conf: a config object to which attributes are added
    """
    # assert that each expected parameter has a unique name
    expected_names = [p.name for p in param_definitions]
    if len(set(expected_names)) < len(expected_names):
        msg = 'expected unique names in the param definitions'
        raise ValueError(msg)
    # map expected parameter names to definitions
    name_to_defn = dict((p.name, p) for p in param_definitions)
    # assert that each parameter provided by the user has a unique name
    observed_names = [name for name, value in name_value_pairs]
    if len(set(expected_names)) < len(expected_names):
        msg = 'expected unique names in the name value pairs'
        raise ValueError(msg)
    # map parameter names to value strings
    name_to_value = dict(name_value_pairs)
    # assert that no parameter name is invalid or missing
    expected = set(expected_names)
    observed = set(observed_names)
    invalid_params = observed - expected
    missing_params = expected - observed
    if invalid_params:
        raise Exception('invalid params: ' + str(list(invalid_params)))
    if missing_params:
        raise Exception('missing params: ' + str(list(missing_params)))
    # set attributes of the config object
    for name, value in name_value_pairs:
        defn = name_to_defn[name]
        typed_value = defn.ftype(value)
        setattr(conf, name, typed_value)
    return conf
