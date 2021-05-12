import pickle
import time
import fcntl
import sys
import sqlite3
import six


# Singleton class to handle cache
class SimpleCache:
    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance__ is not None:
            return cls.__instance__

        # Py2 compatibility
        if sys.version_info[0] < 3:
            BlockingIOError = IOError

        cls._calculation_data_filename = 'calculation_data.pkl'
        cls._pickle_protocol = pickle.HIGHEST_PROTOCOL

        cls.__instance__ = super(SimpleCache, cls, ).__new__(cls)
        return cls.__instance__

    def __init__(self, filename=None):
        """
        Constructor
        """

        if filename is not None:
            self._calculation_data_filename = filename

        # python 2 compatibility
        if not '_calculation_data_filename' in dir(self):
            self._calculation_data_filename = 'calculation_data.db'

        try:
            with open(self._calculation_data_filename, 'rb') as input:
                self._calculation_data = pickle.load(input)
                print('Loaded data from {}'.format(self._calculation_data_filename))
        except (IOError, EOFError, BlockingIOError):
            print('Creating new calculation data file {}'.format(self._calculation_data_filename))
            self._calculation_data = {}
        except (UnicodeDecodeError):
            print('Warning: Calculation data file is corrupted and will be overwritten')
            self._calculation_data = {}

    def redefine_calculation_data_filename(self, filename):

        self._calculation_data_filename = filename
        print('Set data file to {}'.format(self._calculation_data_filename))

        try:
            with open(self._calculation_data_filename, 'rb') as input:
                self._calculation_data = pickle.load(input)
                print('Loaded data from {}'.format(self._calculation_data_filename))
        except (IOError, EOFError):
            print('Creating new calculation data file {}'.format(self._calculation_data_filename))
            self._calculation_data = {}

    def store_calculation_data(self, input_qchem, keyword, data, timeout=60):

        for iter in range(100):
            try:
                with open(self._calculation_data_filename, 'rb') as input:
                    self._calculation_data = pickle.load(input)
            except FileNotFoundError:
                self._calculation_data = {}
                continue
            except (UnicodeDecodeError):
                print('Warning: {} file is corrupted and will be overwritten'.format(self._calculation_data_filename))
                self._calculation_data = {}
            except (BlockingIOError, IOError, EOFError):
                # print('read_try: {}'.format(iter))
                time.sleep(timeout/100)
                continue
            break

        self._calculation_data[(hash(input_qchem), keyword)] = data

        for iter in range(100):
            try:
                with open(self._calculation_data_filename, 'wb') as f:
                    fcntl.lockf(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    pickle.dump(self._calculation_data, f, self._pickle_protocol)
            except BlockingIOError:
                # print('read_try: {}'.format(iter))
                time.sleep(timeout/100)
                continue
            break

    def retrieve_calculation_data(self, input_qchem, keyword):
        return self._calculation_data[(hash(input_qchem), keyword)] if (hash(input_qchem), keyword) in self._calculation_data else None

    @property
    def calculation_data(self):
        return self._calculation_data

    @calculation_data.setter
    def calculation_data(self, calculation_data):
        self._calculation_data = calculation_data


class SqlCache:
    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance__ is not None:
            return cls.__instance__

        cls._calculation_data_filename = 'calculation_data.db'

        cls.__instance__ = super(SqlCache, cls, ).__new__(cls)
        return cls.__instance__

    def __init__(self, filename=None):
        """
        Constructor
        """

        if filename is not None:
            self._calculation_data_filename = filename

        # python 2 compatibility
        if not '_calculation_data_filename' in dir(self):
            print('python 2')
            self._calculation_data_filename = 'calculation_data.db'

        self._conn = sqlite3.connect(self._calculation_data_filename)

        try:
            self._conn.execute('''CREATE TABLE DATA_TABLE
                              (input_hash  INT ,
                               parser      TEXT,
                               qcdata      TEXT);''')
            self._conn.commit()
            # print('Initialized database')

        except sqlite3.OperationalError as e:
            if str(e) != 'table DATA_TABLE already exists':
                raise e

        self._conn.close()

    def __del__(self):
        self._conn.close()

    def redefine_calculation_data_filename(self, filename):
        self._calculation_data_filename = filename
        self.__init__()

    def store_calculation_data(self, input_qchem, keyword, data):

        self._conn = sqlite3.connect(self._calculation_data_filename)


        serialized_data = pickle.dumps(data, protocol=2)

        # python 2 compatibility
        if sys.version_info[0] < 3:
            serialized_data = buffer(serialized_data)

        self._conn.execute("INSERT or REPLACE into DATA_TABLE (input_hash, parser, qcdata)  VALUES (?, ?, ?)",
                           (hash(input_qchem),  keyword, serialized_data))
        self._conn.commit()
        self._conn.close()

    def retrieve_calculation_data(self, input_qchem, keyword):

        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("SELECT qcdata FROM DATA_TABLE WHERE input_hash=? AND parser=?",
                                    (hash(input_qchem), keyword))
        rows = cursor.fetchall()

        self._conn.close()

        return pickle.loads(rows[0][0]) if len(rows) > 0 else None

    @property
    def calculation_data(self):

        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("SELECT input_hash, parser, qcdata from DATA_TABLE")

        self._calculation_data = {}
        for row in cursor:
            self._calculation_data[(row[0], row[1])] = pickle.loads(row[2])

        self._conn.close()

        return self._calculation_data

    @calculation_data.setter
    def calculation_data(self, calculation_data):

        self._conn = sqlite3.connect(self._calculation_data_filename)

        for key, value in calculation_data.items():
            self._conn.execute("INSERT or REPLACE into DATA_TABLE (input_hash, parser, qcdata)  VALUES (?, ?, ?)",
                               (key[0], key[1], pickle.dumps(value, protocol=2)))

        self._conn.commit()
        self._conn.close()


if __name__ == '__main__':
    a = SqlCache()
    b = SqlCache()

    #b.redefine_calculation_data_filename('calculation_data2.db')

    from pyqchem import QchemInput, Structure

    input = QchemInput(Structure(coordinates=[[0, 0, 0]], symbols=['X']))

    b.store_calculation_data(input, 'key1', {'entry1': 454, 'entry2': 2323})
    data = b.retrieve_calculation_data(input, 'key1')
    print(data)