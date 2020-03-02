class ParserError(Exception):
    def __init__(self, parser_name, message):
        self.parser_name = parser_name
        self.message = message

    def __str__(self):
        return 'Error found while parsing output using parser {}: {}'.format(self.parser_name, self.message)


class OutputError(Exception):
    def __init__(self, output, error_lines):
        self.full_output = output
        self.error_lines = error_lines

    def __str__(self):
        return 'Error in Q-Chem calculation:\n{}'.format(self.error_lines)

class QchemInputWarning(UserWarning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message