from distutils.core import setup

# Make python package
setup(name='pyqchem',
      version=0.1,
      description='qchem python module',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['pyqchem', 'pyqchem.parsers'],
      scripts=['scripts/diabatic.py',
               'scripts/plot_diabatic_1d.py',
               'scripts/plot_diabatic_3d.py'])