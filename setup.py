from setuptools import setup

# Make python package
setup(name='pyqchem',
      version=0.1,
      description='qchem python module',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['pyqchem', 'pyqchem.parsers'],
      url='https://github.com/abelcarreras/PyQchem',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License",
          ],
      scripts=['scripts/diabatic.py',
               'scripts/plot_diabatic_1d.py',
               'scripts/plot_diabatic_3d.py'])
