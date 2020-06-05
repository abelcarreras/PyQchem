class ParserError(Exception):
    def __init__(self, parser_name, message):
        self.parser_name = parser_name
        self.message = message

    def __str__(self):
        return 'Error found while parsing output using "{}" parser: {}'.format(self.parser_name, self.message)


class OutputError(Exception):
    def __init__(self, output, error_output):
        self.full_output = output
        self.error_lines = error_output + '\n'.join(output.split('\n')[-20:])

    def __str__(self):
        return 'Error in Q-Chem calculation:\n{}'.format(self.error_lines)


class StructureError(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return 'Error in Structure:\n{}'.format(self._message)

class QchemInputWarning(UserWarning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class QchemInputError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message
