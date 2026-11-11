FIELDS = {'observation_id': 'str', 'sample': 'str', 'line': 'str', 'day': 'int', 'new_value': 'float', 'old_value': 'float', 'operator': 'str', 'lot': 'str', 'valid': 'int'}
REQUIRED = ['observation_id', 'sample', 'line', 'new_value']
ENTITY = 'sample'
GROUP = 'line'
MEASURE = 'new_value'
PRODUCT = 'measurement-equivalence'
UNIT = '微米'
