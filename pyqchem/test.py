import numpy as np

def trunc_dictionary_list(dic_data, decimal):
    # w : dictionary
    # decimal: number of decimal places to leave after truncation
    def round_float(num):
        if isinstance(num, float):
            return int(num * 10**decimal)
        else:
            return num

    def iterlist(d):
        lf = []
        for i, v in enumerate(d):
            if isinstance(v, dict):
                iterdict(v)
            elif isinstance(v, list):
                iterlist(v)
            else:
                d[i] = round_float(v)

    def iterdict(d):
        lf = []
        for k,v in d.items():
            if isinstance(v, dict):
                iterdict(v)
            elif isinstance(v, list):
                iterlist(v)
            else:
                d[k] = round_float(v)

    if isinstance(dic_data, dict):
        iterdict(dic_data)
    elif isinstance(dic_data, list):
        iterlist(dic_data)
    else:
        round_float(dic_data)


def standardize_dictionary(dic_data, decimal=5):

    # set all amplitudes to absolute value
    for state in dic_data['excited_states']:
        for configuration in state['configurations']:
            configuration['amplitude'] = abs(configuration['amplitude'])

    trunc_dictionary_list(dic_data, decimal)

    # delete zero amplitude configurations
    for state in dic_data['excited_states']:
        delete_list = []
        for i, configuration in enumerate(state['configurations']):
            if configuration['amplitude'] == 0:
                delete_list.append(i)
        for i in reversed(delete_list):
            del state['configurations'][i]

    return dic_data


def assertDeepAlmostEqual(test_case, expected, actual, *args, **kwargs):
    """
    Assert that two complex structures have almost equal contents.

    Compares lists, dicts and tuples recursively. Checks numeric values
    using test_case's :py:meth:`unittest.TestCase.assertAlmostEqual` and
    checks all other values with :py:meth:`unittest.TestCase.assertEqual`.
    Accepts additional positional and keyword arguments and pass those
    intact to assertAlmostEqual() (that's how you specify comparison
    precision).

    :param test_case: TestCase object on which we can call all of the basic
    'assert' methods.
    :type test_case: :py:class:`unittest.TestCase` object
    """

    import numpy
    from numpy import long
    is_root = not '__trace' in kwargs
    trace = kwargs.pop('__trace', 'ROOT')
    try:
        if isinstance(expected, (int, float, long, complex)):
            test_case.assertAlmostEqual(expected, actual, *args, **kwargs)
        elif isinstance(expected, (list, tuple, numpy.ndarray)):
            test_case.assertEqual(len(expected), len(actual))
            for index in range(len(expected)):
                v1, v2 = expected[index], actual[index]
                assertDeepAlmostEqual(test_case, v1, v2,
                                      __trace=repr(index), *args, **kwargs)
        elif isinstance(expected, dict):
            test_case.assertEqual(set(expected), set(actual))
            for key in expected:
                assertDeepAlmostEqual(test_case, expected[key], actual[key],
                                      __trace=repr(key), *args, **kwargs)
        else:
            test_case.assertEqual(expected, actual)
    except AssertionError as exc:
        exc.__dict__.setdefault('traces', []).append(trace)
        if is_root:
            trace = ' -> '.join(reversed(exc.traces))
            exc = AssertionError("%s\nTRACE: %s" % (exc.message, trace))
        raise exc
