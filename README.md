# pyreloc
A simple project to perform earthquake relocation (for student projets).

## Dependencies
- python2 or python3
- numpy
- scipy.signal
- matplotlib.pyplot

## Some instructions
To get pyreloc, use the following git command:
```
git clone https://github.com/eost/pyreloc.git
```
or download the module directly from [https://github.com/eost/pyreloc.git].

The pyreloc directory must be placed in a directory pointed by the PYTHONPATH environment variable. You can modify the environment variable in your .bashrc using
```
export PYTHONPATH=path_to_the_directory_including_pyreloc:$PYTHONPATH
```

To use pyreloc, simply import the pyreloc module
```
import pyreloc
```
or if you want to import the seismogram class:
```
from pyreloc import seismogram
```

For now, the module only includes the seismogram class (the class used to deal with seismograms).
