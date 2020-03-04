from setuptools import setup

def get_version_number():
    for l in open('pyqchem/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


# Make python package
setup(name='pyqchem',
      version=get_version_number(),
      description='qchem python module',
      install_requires=['numpy', 'scipy', 'lxml', 'requests', 'matplotlib', 'PyYAML'],
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['pyqchem', 'pyqchem.parsers'],
      url='https://github.com/abelcarreras/PyQchem',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
