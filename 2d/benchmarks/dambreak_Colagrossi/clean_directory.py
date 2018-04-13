import os

def clean_directory():
    """A preprocessing command to clean the current directory. """
    rdomain_ext = ('edge','ele','neig','node','poly','prof0','info',
                   'm','log','h5','xmf')
    for currentFile in os.listdir('.'):
        if any(currentFile.endswith(ext) for ext in rdomain_ext):
            os.remove(currentFile)

clean_directory()
