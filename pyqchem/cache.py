import pickle
import time
import sys
import sqlite3
import numpy as np
from datetime import datetime
import zlib


# Singleton classes to handle cache


class SimpleCache(object):
    """
    Class that makes use of pickle module to hadle data cache

    :param filename: name of the file storing the cache data
    """

    class __SimpleCache:
        def __init__(self, filename='calculation_data.pkl', compress=False):
            self._calculation_data_filename = filename
            self._pickle_protocol = pickle.HIGHEST_PROTOCOL
            self._compress = compress

            # Py2 compatibility
            if 'BlockingIOError' not in vars():
                BlockingIOError = IOError

            try:
                with open(self._calculation_data_filename, 'rb') as input:
                    self.calculation_data = pickle.load(input)
                    print('Loaded data from {}'.format(self._calculation_data_filename))
            except (IOError, EOFError, BlockingIOError):
                print('Creating new calculation data file {}'.format(self._calculation_data_filename))
                self.calculation_data = {}
            except (UnicodeDecodeError):
                print('Warning: Calculation data file is corrupted and will be overwritten')
                self.calculation_data = {}

        def redefine_calculation_data_filename(self, filename):
            """
            redefine the file that contains the cache

            :param filename: cache file name
            :return:
            """

            self._calculation_data_filename = filename
            print('Set data file to {}'.format(self._calculation_data_filename))

            try:
                with open(self._calculation_data_filename, 'rb') as input:
                    self.calculation_data = pickle.load(input)
                    print('Loaded data from {}'.format(self._calculation_data_filename))
            except (IOError, EOFError):
                print('Creating new calculation data file {}'.format(self._calculation_data_filename))
                self.calculation_data = {}

        def store_calculation_data(self, input_qchem, keyword, data, timeout=60):
            """
            Store data in the cache file

            :param input_qchem: QchemInput instance
            :param keyword: string that will be used as a key to store the data
            :param data: the data dictionary to be stored
            :param timeout: time out in seconds to wait if other precesses are accesing the file at the same time
            :return:
            """

            # Py2 compatibility
            if 'FileNotFoundError' not in vars():
                FileNotFoundError = IOError

            for iter in range(100):
                try:
                    with open(self._calculation_data_filename, 'rb') as input:
                        self.calculation_data = pickle.load(input)
                except FileNotFoundError:
                    self.calculation_data = {}
                    continue
                except (UnicodeDecodeError):
                    print(
                        'Warning: {} file is corrupted and will be overwritten'.format(self._calculation_data_filename))
                    self.calculation_data = {}
                except (BlockingIOError, IOError, EOFError):
                    # print('read_try: {}'.format(iter))
                    time.sleep(timeout / 100)
                    continue
                break

            if self._compress:
                data = zlib.compress(pickle.dumps(data, protocol=2))

            self.calculation_data[(hash(input_qchem), keyword)] = data

            for iter in range(100):
                try:
                    with open(self._calculation_data_filename, 'wb') as f:
                        if sys.platform not in ["win32", "cygwin"]:
                            import fcntl
                            fcntl.lockf(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        pickle.dump(self.calculation_data, f, self._pickle_protocol)
                except BlockingIOError:
                    # print('read_try: {}'.format(iter))
                    time.sleep(timeout / 100)
                    continue
                break

        def retrieve_calculation_data(self, input_qchem, keyword):
            """
            retrieve calculation data from cache file

            :param input_qchem: QchemInput instance
            :param keyword: string that was used as a key to store the data
            :return:
            """
            data = self.calculation_data[(hash(input_qchem), keyword)] if (hash(input_qchem),
                                                                           keyword) in self.calculation_data else None

            if self._compress and data is not None:
                data = pickle.loads(zlib.decompress(data))

            return data

        def get_all_data(self):
            """
            return a list of all data stored in the cache file

            :return: list of data
            """

            calc_id_list = np.unique([r[0] for r in self.calculation_data.keys()])
            calc_list = []
            for id in calc_id_list:
                data_dict = {}
                for r in self.calculation_data:
                    if r[0] == id:
                        value = self.calculation_data[(id, r[1])]
                        if self._compress:
                            data_dict.update({r[1]: pickle.loads(zlib.decompress(value))})
                        else:
                            data_dict.update({r[1]: value})
                calc_list.append(data_dict)

            return calc_list

    instance = None

    def __new__(cls, **arguments):
        if not SimpleCache.instance:
            SimpleCache.instance = SimpleCache.__SimpleCache(**arguments)
        return SimpleCache.instance

    def __getattr__(self, nombre):
        return getattr(self.instance, nombre)

    def __setattr__(self, nombre, valor):
        return setattr(self.instance, nombre, valor)


class SqlCache:
    """
    Class that makes use a SQLlite database to hadle data cache

    :param filename: name of the database file storing the cache data
    """
    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance__ is not None:
            return cls.__instance__

        cls._calculation_data_filename = 'calculation_data.db'

        cls.__instance__ = super(SqlCache, cls, ).__new__(cls)
        return cls.__instance__

    def __init__(self, filename=None, compress=False):
        """
        Constructor
        """
        self._compress = compress

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
                               parser      LONGTEXT,
                               qcdata      LONGTEXT,
                               date        LONGTEXT);''')
            self._conn.commit()
            # print('Initialized database')

        except sqlite3.OperationalError as e:
            if str(e) != 'table DATA_TABLE already exists':
                raise e

        self._conn.close()

    def __del__(self):
        self._conn.close()

    def redefine_calculation_data_filename(self, filename):
        """
        redefine the file that contains the cache

        :param filename: cache file name
        :return:
        """
        self._calculation_data_filename = filename
        self.__init__()

    def store_calculation_data(self, input_qchem, keyword, data):
        """
        Store data in the cache file

        :param input_qchem: QchemInput instance
        :param keyword: string that will be used as a key to store the data
        :param data: the data dictionary to be stored
        :param timeout: time out in seconds to wait if other precesses are accesing the file at the same time
        :return:
        """

        date_time = datetime.now()

        self._conn = sqlite3.connect(self._calculation_data_filename)

        serialized_data = pickle.dumps(data, protocol=2)

        if self._compress:
            serialized_data = zlib.compress(serialized_data)

        # python 2 compatibility
        if sys.version_info[0] < 3:
            serialized_data = buffer(serialized_data) # noqa

        self._conn.execute("DELETE FROM DATA_TABLE WHERE input_hash=? AND parser=?",
                           (hash(input_qchem), keyword))

        try:
            self._conn.execute("INSERT into DATA_TABLE (input_hash, parser, qcdata, date)  VALUES (?, ?, ?, ?)",
                               (hash(input_qchem),  keyword, serialized_data, date_time))
        except sqlite3.OperationalError:
            self._conn.execute("INSERT into DATA_TABLE (input_hash, parser, qcdata)  VALUES (?, ?, ?)",
                               (hash(input_qchem),  keyword, serialized_data))

        self._conn.commit()
        self._conn.close()

    def retrieve_calculation_data(self, input_qchem, keyword):
        """
        retrieve calculation data from cache file

        :param input_qchem: QchemInput instance
        :param keyword: string that was used as a key to store the data
        :return:
        """
        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("SELECT qcdata FROM DATA_TABLE WHERE input_hash=? AND parser=?",
                                    (hash(input_qchem), keyword))
        rows = cursor.fetchall()

        self._conn.close()

        if self._compress:
            return pickle.loads(zlib.decompress(rows[0][0])) if len(rows) > 0 else None

        return pickle.loads(rows[0][0]) if len(rows) > 0 else None

    def retrieve_calculation_data_from_id(self, id, keyword=None):
        """
        return data using database entry ID
        [Only for SQL database cache]

        :param id: databse entry ID
        :param keyword: string that was used as a key to store the data
        :return:
        """

        self._conn = sqlite3.connect(self._calculation_data_filename)

        if keyword is None:
            cursor = self._conn.execute("SELECT qcdata FROM DATA_TABLE WHERE input_hash=?", (id,))
            rows = cursor.fetchall()
        else:
            cursor = self._conn.execute("SELECT qcdata FROM DATA_TABLE WHERE input_hash=? AND parser=?",
                                        (id, keyword))
            rows = cursor.fetchall()

        self._conn.close()

        if len(rows) <= 0:
            return None
        elif len(rows) == 1:
            if self._compress:
                return pickle.loads(zlib.decompress(rows[0][0])) if len(rows) > 0 else None
            else:
                return pickle.loads(rows[0][0]) if len(rows) > 0 else None
        else:
            if self._compress:
                return [pickle.loads(zlib.decompress(r[0])) for r in rows]
            else:
                return [pickle.loads(r[0]) for r in rows]

    def list_database(self):
        """
        prints data inside database
        [Only for SQL databse cache]

        :return: None
        """
        self._conn = sqlite3.connect(self._calculation_data_filename)

        try:
            cursor = self._conn.execute("SELECT input_hash, parser, date from DATA_TABLE")
        except sqlite3.OperationalError:
            cursor = self._conn.execute("SELECT input_hash, parser from DATA_TABLE")


        print('{:^25} {:^25} {:^25}'.format('ID', 'KEYWORD', 'DATE'))
        print('--'*40)
        for row in cursor:
            try:
                print('{:<25} {:<25} {}'.format(*row))
            except IndexError:
                print('{:<25} {:<25}'.format(*row))

        self._conn.close()

    def integrity_check(self):
        """
        Check integrity of the database
        [Only for SQL databse cache]

        :return:
        """
        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("PRAGMA integrity_check")
        rows = ''.join(*cursor.fetchall()[0])
        print(rows)
        self._conn.close()

        pass

    def fix_database(self, filename):
        """
        fix correupted database and store the recovered data in a new recovered database file
        [Only for SQL databse cache]

        :param filename: recovered database filename
        :return:
        """

        import subprocess, os

        dump_file = self._calculation_data_filename + '.dump'
        # dump_file = '_recovery.test'

        schema = subprocess.run(
            ['sqlite3',
             self._calculation_data_filename,
             '.output {}'.format(dump_file),
             '.dump',
             ],
            capture_output=True
        )

        with open(dump_file, 'r') as f:
            data = f.read().replace("ROLLBACK", "COMMIT")

        with open(dump_file, 'w') as f:
            f.write(data)

        try:
            os.remove(filename)
        except FileNotFoundError:
            pass

        schema = subprocess.run(
            ['sqlite3',
             filename,
             '.read {}'.format(dump_file)
             ],
            capture_output=True
        )

        os.remove(dump_file)

    def _recovery_kill(self, file):

        # cursor = self._conn.execute(".save ?", (file,))

        import subprocess
        schema = subprocess.run(
            ['sqlite3',
             self._calculation_data_filename,
             '.recover'.format(file)
             ],
            capture_output=True
        )

        print(schema.stdout)

    def get_all_data(self):
        """
        return a list of all data stored in the cache file

        :return: list of data
        """

        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("SELECT * FROM DATA_TABLE")
        rows = cursor.fetchall()

        self._conn.close()

        calc_id_list = np.unique([r[0] for r in rows])

        calc_list = []
        for id in calc_id_list:
            data_dict = {}
            for r in rows:
                if r[0] == id:
                    if self._compress:
                        data_dict.update({r[1]: pickle.loads(zlib.decompress(r[2]))})
                    else:
                        data_dict.update({r[1]: pickle.loads(r[2])})
            calc_list.append(data_dict)

        return calc_list

    @property
    def calculation_data(self):

        self._conn = sqlite3.connect(self._calculation_data_filename)

        cursor = self._conn.execute("SELECT input_hash, parser, qcdata from DATA_TABLE")

        self._calculation_data = {}
        for row in cursor:

            if self._compress:
                self._calculation_data[(row[0], row[1])] = pickle.loads(zlib.decompress(row[2]))
            else:
                self._calculation_data[(row[0], row[1])] = pickle.loads(row[2])

        self._conn.close()

        return self._calculation_data

    @calculation_data.setter
    def calculation_data(self, calculation_data):

        self._conn = sqlite3.connect(self._calculation_data_filename)

        for key, value in calculation_data.items():

            serialized_data = pickle.dumps(value, protocol=2)
            if self._compress:
                serialized_data = zlib.compress(serialized_data)

            self._conn.execute("INSERT or REPLACE into DATA_TABLE (input_hash, parser, qcdata)  VALUES (?, ?, ?)",
                               (key[0], key[1], serialized_data))

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